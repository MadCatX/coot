// vim: sw=4 ts=4 sts=4 expandtab :

#include "ntc.hh"
#include "util.hh"
#include "ui/util.hh"

#include "src/molecule-class-info.h"

#include <LLKA/llka_nucleotide.h>
#include <LLKA/llka_resource_loaders.h>
#include <LLKA/llka_superposition.h>
#include <LLKA/llka_util.h>

#include <cassert>
#include <cmath>
#include <mutex>

static const std::string CLUSTERS_FILE{"clusters.csv"};
static const std::string CONFALS_FILE{"confals.csv"};
static const std::string GOLDEN_STEPS_FILE{"golden_steps.csv"};
static const std::string NU_ANGLES_FILE{"nu_angles.csv"};

static LLKA_ClassificationContext *classificationContext{nullptr};
static std::mutex initializationLock;

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
    limits.manhattanDistanceCutoff = LLKA_deg2rad(60.0);
    limits.pseudorotationCutoff = LLKA_deg2rad(72.0);
    limits.minimumClusterVotes = 0.001111;
    limits.minimumNearestNeighbors = 7;
    limits.numberOfUsedNearestNeighbors = 11;

    goldenSteps.type = LLKA_RES_GOLDEN_STEPS;
    clusters.type = LLKA_RES_CLUSTERS;
    confals.type = LLKA_RES_CONFALS;
    nuAngles.type = LLKA_RES_AVERAGE_NU_ANGLES;

    /* Load the data necessary to initialize the classification context */
    LLKA_RetCode tRet = LLKA_loadResource((path + "/" + GOLDEN_STEPS_FILE).c_str(), LLKA_RES_AS_FILE, &goldenSteps);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load golden steps data: "} + LLKA_errorToString(tRet);
        goto out;
    }
    tRet = LLKA_loadResource((path + "/" + CLUSTERS_FILE).c_str(), LLKA_RES_AS_FILE, &clusters);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load clusters data: "} + LLKA_errorToString(tRet);
        goto out_golden_steps;
    }
    tRet = LLKA_loadResource((path + "/" + CONFALS_FILE).c_str(), LLKA_RES_AS_FILE, &confals);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load confals data: "} + LLKA_errorToString(tRet);
        goto out_clusters;
    }
    tRet = LLKA_loadResource((path + "/" + NU_ANGLES_FILE).c_str(), LLKA_RES_AS_FILE, &nuAngles);
    if (tRet != LLKA_OK) {
        error = std::string{"Failed to load Nu angles data: "} + LLKA_errorToString(tRet);
        goto out_confals;
    }

    /* Try to initialize the context */
    tRet = LLKA_initializeClassificationContext(
        clusters.data.clusters, clusters.count,
        goldenSteps.data.goldenSteps, goldenSteps.count,
        confals.data.confals, confals.count,
        nuAngles.data.averageNuAngles, nuAngles.count,
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

NtCStructure ntc_dinucleotide_from_atom(int atom_index, int imol, const std::vector<molecule_class_info_t> &molecules) {
    auto &molecule = molecules[imol];

    mmdb::Atom *atom = molecule.atom_sel.atom_selection[atom_index];
    mmdb::Residue *residue = atom->residue;
    const std::string &altconf = atom->altLoc;
    mmdb::Residue *residue2 = molecule.get_following_residue(coot::residue_spec_t(residue));
    if (!residue2)
        return {};

    if (!is_nucleotide(residue->GetLabelCompID()) || !is_nucleotide(residue2->GetLabelCompID()))
        return {};

    mmdb::Residue *filteredResidue = clone_mmdb_residue(residue, altconf);
    mmdb::Residue *filteredResidue2 = clone_mmdb_residue(residue2, altconf);

    std::cout << filteredResidue->GetLabelCompID() << ", " << filteredResidue2->GetLabelCompID() << "\n";

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

    LLKA_Structure llkaStru = mmdb_structure_to_LLKA_structure(mmdbStru);
    if (llkaStru.nAtoms == 0)
        return {};

    return { mmdbStru, llkaStru };
}

NtCSimilarityResult ntc_calculate_similarities(const NtCStructure &stru) {
    static std::vector<LLKA_NtC> AllNtCs = []() {
        auto r = make_ntc_range(LLKA_AA00, LLKA_LAST_NTC);
        r.push_back(LLKA_INVALID_NTC);
        return r;
    }();

    LLKA_Similarities cSimils{};
    cSimils.similars = new LLKA_Similarity[AllNtCs.size() - 1];
    cSimils.nSimilars = AllNtCs.size() - 1;

    LLKA_RetCode tRet = LLKA_measureStepSimilarityNtCMultiple(&stru.llkaStru, AllNtCs.data(), &cSimils);
    if (tRet != LLKA_OK) {
        return NtCSimilarityResult::fail(tRet);
    }

    std::vector<NtCSimilarity> simils{};
    for (size_t idx = 0; idx < cSimils.nSimilars; idx++) {
        simils.emplace_back(cSimils.similars[idx], LLKA_NtCToName(AllNtCs[idx]));
    }

    // WTF? Is there no deleter for cSimils???

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
