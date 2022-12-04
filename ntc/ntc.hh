// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_NTC_H
#define _NTC_NTC_H

#ifdef COOT_ENABLE_NTC

#include "types.hh"

#include <LLKA/llka_connectivity_similarity.h>
#include <LLKA/llka_classification.h>
#include <LLKA/llka_structure.h>

#include <string>
#include <utility>
#include <vector>

class NtCStep {
public:
    NtCStep() = default;
    NtCStep(NtCStepAltConf altconf, NtCStructure stru) noexcept :
        altconf{std::move(altconf)},
        stru{std::move(stru)}
    {}
    NtCStep(const NtCStep &) = delete;
    NtCStep(NtCStep &&other) noexcept :
        altconf{std::move(other.altconf)},
        stru{std::move(other.stru)}
    {}

    NtCStepAltConf altconf;
    NtCStructure stru;

    bool isValid() const { return stru.isValid; }
};
using NtCSteps = std::vector<NtCStep>;

class NtCSuperposition {
public:
    mmdb::Manager *mmdbStru;
    double rmsd;
};

NtCResult<LLKA_ClassifiedStep, LLKA_RetCode> ntc_classify(const NtCStructure &stru);
NtCSteps ntc_dinucleotides(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue);
NtCConnectivitiesResult ntc_calculate_connectivities(LLKA_NtC ntc, const NtCStep &step, mmdb::Manager *srcMmdbStru);
NtCSimilaritiesResult ntc_calculate_similarities(const NtCStructure &stru);
NtCStructure ntc_get_reference_structure(LLKA_NtC ntc);
bool ntc_initialize_classification_context(const std::string &path, std::string &error);
bool ntc_initialize_classification_context_if_needed(std::string path, std::string &error);
NtCSuperposition ntc_superpose_reference(const NtCStructure &stru, LLKA_NtC ntc);

#endif // COOT_ENABLE_NTC

#endif // _NTC_NTC_H
