// vim: sw=4 ts=4 sts=4 expandtab :

#include "ntc.hh"
#include "util.hh"
#include "ui/util.hh"

#include "coot-utils/coot-coord-utils.hh"

#include <mmdb2/mmdb_manager.h>

#include <LLKA/llka_nucleotide.h>
#include <LLKA/llka_resource_loaders.h>
#include <LLKA/llka_superposition.h>
#include <LLKA/llka_util.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <mutex>

enum class RelatedSteps {
    Previous,
    Next
};

static const std::string CLUSTERS_FILE{"clusters.csv"};
static const std::string CONFALS_FILE{"confals.csv"};
static const std::string GOLDEN_STEPS_FILE{"golden_steps.csv"};
static const std::string NU_ANGLES_FILE{"nu_angles.csv"};

static LLKA_ClassificationContext *classificationContext{nullptr};
static std::mutex initializationLock;

static const std::vector<LLKA_NtC> AllNtCs = []() {
    auto r = make_ntc_range(LLKA_AA00, LLKA_LAST_NTC);
    r.push_back(LLKA_INVALID_NTC);
    return r;
}();

static
bool is_nucleotide(const char *compId) {
    LLKA_Bool res = LLKA_isNucleotideCompound(compId);
    return res == LLKA_TRUE;
}

static
void destroy_connectivities(LLKA_Connectivities &conns) {
    delete [] conns.conns;
}

static
void destroy_similarities(LLKA_Similarities &simils) {
    delete [] simils.similars;
}

static
NtCSteps expand_residue_to_steps(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue, const std::string &firstAltconf = "", const std::string &secondAltconf = "") {
    if (!is_nucleotide(residue->GetLabelCompID())) {
        return {};
    }

    mmdb::Residue *residue2 = coot::util::get_following_residue(coot::residue_spec_t(residue), srcMmdbStru);
    if (!residue2 || !is_nucleotide(residue2->GetLabelCompID())) {
        return {};
    }

    std::vector<std::string> altconfs1 = firstAltconf.empty() ? all_altconfs(residue) : std::vector<std::string>{ firstAltconf };
    std::vector<std::string> altconfs2 = secondAltconf.empty() ? all_altconfs(residue2) : std::vector<std::string>{ secondAltconf };

    if (altconfs1.empty()) {
        altconfs1.push_back("");
    }
    if (altconfs2.empty()) {
        altconfs2.push_back("");
    }

    NtCSteps steps;
    for (const auto &ac1 : altconfs1) {
        mmdb::Residue *filteredResidue1 = clone_mmdb_residue(residue, ac1);

        for (const auto &ac2 : altconfs2) {
            mmdb::Residue *filteredResidue2 = clone_mmdb_residue(residue2, ac2);

            mmdb::Manager *mmdbStru = new mmdb::Manager;
            mmdb::Model *model = new mmdb::Model;
            mmdb::Chain *chain = new mmdb::Chain;

            chain->AddResidue(clone_mmdb_residue(filteredResidue1)); // We need a new residue for each variant
            chain->AddResidue(filteredResidue2);
            chain->SetChainID(residue->GetChainID());
            model->AddChain(chain);
            mmdbStru->AddModel(model);
            mmdbStru->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);
            mmdbStru->FinishStructEdit();

            LLKA_Structure llkaStru = mmdb_structure_to_LLKA_structure(mmdbStru);
            assert(llkaStru.nAtoms == mmdbStru->GetNumberOfAtoms());

            steps.emplace_back(NtCStepAltConf{ac1, ac2}, NtCStructure{mmdbStru, llkaStru});
        }

        delete filteredResidue1;
    }

    // WARNING: We create a bunch of raw pointers here, freeing just the "mol" object
    // hopefully should be enough to reclaim the resources.

    return steps;
}

static
NtCSteps get_related_steps(RelatedSteps which, const NtCStep &step, mmdb::Manager *srcMmdbStru) {
    assert(step.isValid());

    mmdb::Residue **mmdbRes = nullptr;
    int numResidues = 0;
    step.stru.mmdbStru->GetResidueTable(mmdbRes, numResidues);
    if (numResidues != 2 || !mmdbRes) {
        return {};
    }

    NtCSteps steps;
    if (which == RelatedSteps::Next) {
        steps = expand_residue_to_steps(srcMmdbStru, mmdbRes[1], step.altconf.second);
    } else if (which == RelatedSteps::Previous) {
        mmdb::Residue *prevRes = coot::util::get_previous_residue(coot::residue_spec_t(mmdbRes[0]), srcMmdbStru);
        if (!prevRes) {
            return {};
        } else {
            steps = expand_residue_to_steps(srcMmdbStru, prevRes, "", step.altconf.first);
        }
    } else {
        assert(false); // Should never happen
    }

    return steps;
}

static
LLKA_Connectivities init_connectivities(size_t count) {
    return {
        new LLKA_Connectivity[count],
        count
    };
}

static
LLKA_Similarities init_similarities(size_t count) {
    return {
        new LLKA_Similarity[count],
        count
    };
}

static
std::vector<NtCConnectivity> map_connectivities(const LLKA_Connectivities &llkaConns, const std::vector<LLKA_NtC> &NtCs) {
    assert(llkaConns.nConns == NtCs.size() - 1); // -1 because AllNtCs are terminated by LLKA_INVALID_NTC

    std::vector<NtCConnectivity> conns(llkaConns.nConns);
    for (size_t idx = 0; idx < llkaConns.nConns; idx++) {
        conns[idx].connectivity = llkaConns.conns[idx];
        conns[idx].NtC = LLKA_NtCToName(NtCs[idx]);
    }

    return conns;
}

NtCResult<LLKA_ClassifiedStep, LLKA_RetCode> ntc_classify(const NtCStructure &stru) {
    using RT = NtCResult<LLKA_ClassifiedStep, LLKA_RetCode>;

    LLKA_ClassifiedStep classified{};

    auto tRet = LLKA_classifyStep(&stru.llkaStru, classificationContext, &classified);
    if (tRet != LLKA_OK) {
        return RT::fail(tRet); // We need to return some LLKA_RetCode here
    }

    return RT::succeed(classified);
}

NtCStructure ntc_get_reference_structure(LLKA_NtC ntc) {
    assert(ntc != LLKA_INVALID_NTC);

    LLKA_Structure llkaStru = LLKA_NtCStructure(ntc);
    mmdb::Manager *mmdbStru = LLKA_structure_to_mmdb_structure(llkaStru);

    return { mmdbStru, llkaStru };
}

bool ntc_initialize_classification_context(const std::string &path, std::string &error) {
    std::lock_guard<std::mutex> lk{initializationLock};
    bool res = false;

    if (classificationContext)
        LLKA_destroyClassificationContext(classificationContext);

    LLKA_ClassificationLimits limits;
    LLKA_Resource goldenSteps;
    LLKA_Resource clusters;
    LLKA_Resource confals;
    LLKA_Resource nuAngles;

    // TODO: This should be made configurable
    limits.averageNeighborsTorsionCutoff = LLKA_deg2rad(28.0);
    limits.nearestNeighborTorsionsCutoff = LLKA_deg2rad(28.0);
    limits.totalDistanceCutoff = LLKA_deg2rad(60.0);
    limits.pseudorotationCutoff = LLKA_deg2rad(72.0);
    limits.minimumClusterVotes = 0.001111;
    limits.minimumNearestNeighbors = 7;
    limits.numberOfUsedNearestNeighbors = 11;

    goldenSteps.type = LLKA_RES_GOLDEN_STEPS;
    clusters.type = LLKA_RES_CLUSTERS;
    confals.type = LLKA_RES_CONFALS;
    nuAngles.type = LLKA_RES_AVERAGE_NU_ANGLES;

    /* Load the data necessary to initialize the classification context */
    LLKA_RetCode tRet = LLKA_loadResourceFile(LLKAPathConverter<LLKA_PathChar>::convert(path + "/" + GOLDEN_STEPS_FILE).c_str(), &goldenSteps);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load golden steps data: "} + LLKA_errorToString(tRet);
        goto out;
    }
    tRet = LLKA_loadResourceFile(LLKAPathConverter<LLKA_PathChar>::convert(path + "/" + CLUSTERS_FILE).c_str(), &clusters);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load clusters data: "} + LLKA_errorToString(tRet);
        goto out_golden_steps;
    }
    tRet = LLKA_loadResourceFile(LLKAPathConverter<LLKA_PathChar>::convert(path + "/" + CONFALS_FILE).c_str(), &confals);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load confals data: "} + LLKA_errorToString(tRet);
        goto out_clusters;
    }
    tRet = LLKA_loadResourceFile(LLKAPathConverter<LLKA_PathChar>::convert(path + "/" + NU_ANGLES_FILE).c_str(), &nuAngles);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load Nu angles data: "} + LLKA_errorToString(tRet);
        goto out_confals;
    }

    /* Try to initialize the context */
    tRet = LLKA_initializeClassificationContext(
        clusters.data.clusters, clusters.count,
        goldenSteps.data.goldenSteps, goldenSteps.count,
        confals.data.confals, confals.count,
        nuAngles.data.clusterNuAngles, nuAngles.count,
        &limits,
        &classificationContext
    );

    if (tRet != LLKA_OK) {
        error = std::string{"Failed to initialize classification context: "} + LLKA_errorToString(tRet);
    } else {
        res = true;
    }

    LLKA_destroyResource(&nuAngles);
out_confals:
    LLKA_destroyResource(&confals);
out_clusters:
    LLKA_destroyResource(&clusters);
out_golden_steps:
    LLKA_destroyResource(&goldenSteps);
out:

    return res;
}

bool ntc_initialize_classification_context_if_needed(std::string path, std::string &error) {
    if (classificationContext)
        return true;

    if (path.empty()) {
        path = pick_ntc_parameters_directory();
    }
    if (path.empty()) {
        error = "No path to NtC parameters files was specified";
        return false;
    }

    return ntc_initialize_classification_context(path, error);
}

NtCConnectivitiesResult ntc_calculate_connectivities(LLKA_NtC ntc, const NtCStep &step, mmdb::Manager *srcMmdbStru) {
    std::vector<AltConfNtCConnectivities> prevConns;
    NtCSteps prevSteps = get_related_steps(RelatedSteps::Previous, step, srcMmdbStru);
    LLKA_Connectivities cConns = init_connectivities(AllNtCs.size() - 1);
    for (auto &prevStep : prevSteps) {
        LLKA_RetCode tRet = LLKA_measureStepConnectivityNtCsMultipleFirst(&prevStep.stru.llkaStru, AllNtCs.data(), &prevStep.stru.llkaStru, ntc, &cConns);
        if (tRet != LLKA_OK) {
            destroy_connectivities(cConns);
            return NtCConnectivitiesResult::fail(tRet);
        }

        prevConns.emplace_back(prevStep.altconf.first, map_connectivities(cConns, AllNtCs));
    }

    std::vector<AltConfNtCConnectivities> nextConns;
    NtCSteps nextSteps = get_related_steps(RelatedSteps::Next, step, srcMmdbStru);
    for (auto &nextStep : nextSteps) {
        LLKA_RetCode tRet = LLKA_measureStepConnectivityNtCsMultipleSecond(&step.stru.llkaStru, ntc, &nextStep.stru.llkaStru, AllNtCs.data(), &cConns);
        if (tRet != LLKA_OK) {
            destroy_connectivities(cConns);
            return NtCConnectivitiesResult::fail(tRet);
        }

        nextConns.emplace_back(nextStep.altconf.second, map_connectivities(cConns, AllNtCs));
    }

    destroy_connectivities(cConns);

    return NtCConnectivitiesResult::succeed(prevConns, nextConns);
}

NtCSimilaritiesResult ntc_calculate_similarities(const NtCStructure &stru) {
    LLKA_Similarities cSimils = init_similarities(AllNtCs.size() - 1);

    LLKA_RetCode tRet = LLKA_measureStepSimilarityNtCMultiple(&stru.llkaStru, AllNtCs.data(), &cSimils);
    if (tRet != LLKA_OK) {
        return NtCSimilaritiesResult::fail(tRet);
    }

    std::vector<NtCSimilarity> simils{};
    for (size_t idx = 0; idx < cSimils.nSimilars; idx++) {
        simils.emplace_back(cSimils.similars[idx], LLKA_NtCToName(AllNtCs[idx]));
    }
    destroy_similarities(cSimils);

    return NtCSimilaritiesResult::succeed(std::move(simils));
}

NtCSteps ntc_dinucleotides(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue) {
    if (!residue) {
        return {};
    }

    return expand_residue_to_steps(srcMmdbStru, residue);
}

NtCSuperposition ntc_superpose_reference(const NtCStructure &stru, LLKA_NtC ntc) {
    LLKA_Structure llkaRefStru = LLKA_NtCStructure(ntc);
    assert(llkaRefStru.nAtoms > 0);

    // Extract backbones
    LLKA_Structure bkbn{};
    LLKA_Structure refBkbn{};
    auto tRet = LLKA_extractBackbone(&stru.llkaStru, &bkbn);
    if (tRet != LLKA_OK) {
        return {};
    }

    tRet = LLKA_extractBackbone(&llkaRefStru, &refBkbn);
    if (tRet != LLKA_OK) {
        LLKA_destroyStructure(&bkbn);

        return {};
    }

    // Construct transformation matrix
    LLKA_Matrix transformation{};
    tRet = LLKA_superpositionMatrixStructures(&refBkbn, &bkbn, &transformation);
    if (tRet != LLKA_OK) {
        return {};
    }

    tRet = LLKA_applyTransformationStructure(&llkaRefStru, &transformation);
    assert(tRet == LLKA_OK);

    mmdb::Manager *mmdbRefStru = LLKA_structure_to_mmdb_structure(llkaRefStru);

    // We need to relabel the chains and residue numbers to allow coot to replace the current structure
    // instead of just creating a new one
    relabel_mmdb_step(mmdbRefStru, stru.mmdbStru, true);
    // Reference NtC may have different base than the step in the structure. Replace the bases on the reference
    replace_bases(mmdbRefStru, stru.mmdbStru);

    // Superpose the reference backbone onto the actual backbone to allow us to calculate RMSD
    double rmsd;
    tRet = LLKA_applyTransformationStructure(&refBkbn, &transformation);
    assert(tRet == LLKA_OK);
    tRet = LLKA_rmsdStructures(&refBkbn, &bkbn, &rmsd);
    assert(tRet == LLKA_OK);

    LLKA_destroyStructure(&bkbn);
    LLKA_destroyStructure(&refBkbn);

    LLKA_destroyMatrix(&transformation);
    LLKA_destroyStructure(&llkaRefStru);

    return { mmdbRefStru, rmsd };
}
