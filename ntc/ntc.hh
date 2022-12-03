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

class AltConfNtCStep {
public:
    AltConfNtCStep() = default;
    AltConfNtCStep(std::string altconf1, std::string altconf2, NtCStructure stru) noexcept :
        altconf1{std::move(altconf1)},
        altconf2{std::move(altconf2)},
        stru{std::move(stru)}
    {}
    AltConfNtCStep(const AltConfNtCStep &) = delete;
    AltConfNtCStep(AltConfNtCStep &&other) noexcept :
        altconf1{std::move(other.altconf1)},
        altconf2{std::move(other.altconf2)},
        stru{std::move(other.stru)}
    {}

    std::string altconf1;
    std::string altconf2;
    NtCStructure stru;

    bool isValid() const { return stru.isValid; }
};
using AltConfNtCSteps = std::vector<AltConfNtCStep>;

class NtCSuperposition {
public:
    mmdb::Manager *mmdbStru;
    double rmsd;
};

NtCResult<LLKA_ClassifiedStep, LLKA_RetCode> ntc_classify(const NtCStructure &stru);
AltConfNtCSteps ntc_dinucleotides(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue);
NtCConnectivitiesResult ntc_calculate_connectivities(LLKA_NtC ntc, const AltConfNtCStep &step, mmdb::Manager *srcMmdbStru);
NtCSimilaritiesResult ntc_calculate_similarities(const NtCStructure &stru);
NtCStructure ntc_get_reference_structure(LLKA_NtC ntc);
bool ntc_initialize_classification_context(const std::string &path, std::string &error);
bool ntc_initialize_classification_context_if_needed(std::string path, std::string &error);
NtCSuperposition ntc_superpose_reference(const NtCStructure &stru, LLKA_NtC ntc);

#endif // COOT_ENABLE_NTC

#endif // _NTC_NTC_H
