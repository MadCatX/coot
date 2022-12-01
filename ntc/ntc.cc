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

#include <cassert>
#include <cmath>
#include <mutex>

enum class RelatedStep {
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

NtCStructure::NtCStructure() :
    isValid{false},
    m_released{false}
{}

NtCStructure::NtCStructure(mmdb::Manager *mmdbStru, LLKA_Structure llkaStru) :
    mmdbStru{mmdbStru},
    llkaStru{llkaStru},
    isValid{true},
    m_released{false}
{}

NtCStructure::NtCStructure(NtCStructure &&other) noexcept :
    mmdbStru{other.mmdbStru},
    llkaStru{other.llkaStru},
    isValid{other.isValid},
    m_released{other.m_released}
{
    other.release();
}

NtCStructure::~NtCStructure() {
    if (!m_released && isValid) {
        LLKA_destroyStructure(&llkaStru);
        delete mmdbStru;
    }
}

NtCStructure & NtCStructure::operator=(NtCStructure &&other) noexcept {
    this->mmdbStru = other.mmdbStru;
    this->llkaStru = other.llkaStru;
    this->isValid = other.isValid;
    this->m_released = other.m_released;

    other.release();

    return *this;
}

void NtCStructure::release() {
    m_released = true;
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
mmdb::Manager * expand_residue_to_step(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue, const std::string &altconf) {
    if (!is_nucleotide(residue->GetLabelCompID())) {
        return nullptr;
    }

    mmdb::Residue *residue2 = coot::util::get_following_residue(coot::residue_spec_t(residue), srcMmdbStru);
    if (!residue2 || !is_nucleotide(residue2->GetLabelCompID())) {
        return nullptr;
    }

    mmdb::Residue *filteredResidue = clone_mmdb_residue(residue, altconf);
    mmdb::Residue *filteredResidue2 = clone_mmdb_residue(residue2, altconf);

    mmdb::Manager *mmdbStru = new mmdb::Manager;
    mmdb::Model *model = new mmdb::Model;
    mmdb::Chain *chain = new mmdb::Chain;

    chain->AddResidue(filteredResidue);
    chain->AddResidue(filteredResidue2);
    chain->SetChainID(residue->GetChainID());
    model->AddChain(chain);
    mmdbStru->AddModel(model);
    mmdbStru->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);
    mmdbStru->FinishStructEdit();

    // WARNING: We create a bunch of raw pointers here, freeing just the "mol" object
    // hopefully should be enough to reclaim the resources

    return mmdbStru;
}

static
NtCStructure get_related_step(RelatedStep which, const NtCStructure &stru, mmdb::Manager *srcMmdbStru) {
    assert(stru.isValid);

    mmdb::Residue **mmdbRes = nullptr;
    int numResidues = 0;
    stru.mmdbStru->GetResidueTable(mmdbRes, numResidues);
    if (numResidues != 2 || !mmdbRes) {
        return {};
    }

    const mmdb::Residue *res = mmdbRes[0];
    coot::residue_spec_t rs(stru.mmdbStru->GetFirstModelNum(), res->chain->GetChainID(), res->seqNum, res->insCode);

    mmdb::Residue *related = nullptr;
    switch (which) {
    case RelatedStep::Next:
        related = coot::util::get_following_residue(rs, srcMmdbStru);
        break;
    case RelatedStep::Previous:
        related = coot::util::get_previous_residue(rs, srcMmdbStru);
        break;
    default:
        return {};
    }

    if (!related) {
        return {};
    }

    // TODO: Handle altconfs
    mmdb::Manager *mmdbStru = expand_residue_to_step(srcMmdbStru, related, "");
    if (!mmdbStru) {
        return {};
    }

    LLKA_Structure llkaStru = mmdb_structure_to_LLKA_structure(mmdbStru);
    if (llkaStru.nAtoms == 0) {
        delete mmdbStru;
        return {};
    }

    return { mmdbStru, llkaStru };
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

NtCStructure ntc_dinucleotide(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue, const std::string &altconf) {
    if (!residue) {
        return {};
    }

    mmdb::Manager *mmdbStru = expand_residue_to_step(srcMmdbStru, residue, altconf);
    if (!mmdbStru) {
        return {};
    }

    LLKA_Structure llkaStru = mmdb_structure_to_LLKA_structure(mmdbStru);
    if (llkaStru.nAtoms == 0) {
        delete mmdbStru;
        return {};
    }

    return { mmdbStru, llkaStru };
}

NtCConnectivityResult ntc_calculate_connectivity(LLKA_NtC ntc, const NtCStructure &stru, mmdb::Manager *srcMmdbStru) {
    LLKA_Connectivities cPrevConns = init_connectivities(AllNtCs.size() - 1);
    NtCStructure prevStep = get_related_step(RelatedStep::Previous, stru, srcMmdbStru);
    if (prevStep.isValid) {
        LLKA_RetCode tRet = LLKA_measureStepConnectivityNtCsMultipleFirst(&prevStep.llkaStru, AllNtCs.data(), &stru.llkaStru, ntc, &cPrevConns);
        if (tRet != LLKA_OK) {
            NtCConnectivityResult::fail(tRet);
        }
    }

    LLKA_Connectivities cNextConns = init_connectivities(AllNtCs.size() - 1);
    NtCStructure nextStep = get_related_step(RelatedStep::Next, stru, srcMmdbStru);
    if (nextStep.isValid) {
        LLKA_RetCode tRet = LLKA_measureStepConnectivityNtCsMultipleSecond(&stru.llkaStru, ntc, &nextStep.llkaStru, AllNtCs.data(), &cNextConns);
        if (tRet != LLKA_OK) {
            NtCConnectivityResult::fail(tRet);
        }
    }

    NtCConnectivityResult ret = NtCConnectivityResult::succeed(
        prevStep.isValid ? map_connectivities(cPrevConns, AllNtCs) : std::vector<NtCConnectivity>(),
        nextStep.isValid ? map_connectivities(cNextConns, AllNtCs) : std::vector<NtCConnectivity>()
    );

    destroy_connectivities(cPrevConns);
    destroy_connectivities(cNextConns);;

    return ret;
}

NtCSimilarityResult ntc_calculate_similarities(const NtCStructure &stru) {
    LLKA_Similarities cSimils = init_similarities(AllNtCs.size() - 1);

    LLKA_RetCode tRet = LLKA_measureStepSimilarityNtCMultiple(&stru.llkaStru, AllNtCs.data(), &cSimils);
    if (tRet != LLKA_OK) {
        return NtCSimilarityResult::fail(tRet);
    }

    std::vector<NtCSimilarity> simils{};
    for (size_t idx = 0; idx < cSimils.nSimilars; idx++) {
        simils.emplace_back(cSimils.similars[idx], LLKA_NtCToName(AllNtCs[idx]));
    }
    destroy_similarities(cSimils);

    return NtCSimilarityResult::succeed(std::move(simils));
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
