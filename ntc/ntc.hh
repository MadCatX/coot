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

namespace mmdb {
    class Manager;
    class Residue;
}

template <typename Success, typename Failure>
class NtCResult {
public:
    // TODO: We could do this with a union

    Success success;
    Failure failure;
    bool succeeded;

    static NtCResult fail(Failure f) {
        return NtCResult({}, std::move(f), false);
    }

    static NtCResult succeed(Success s) {
        return NtCResult(std::move(s), {}, true);
    }

    template <typename ...Args>
    static NtCResult succeed(Args&& ...args) noexcept {
        return NtCResult(Success{std::forward<Args>(args)...}, {}, true);
    }

private:
    NtCResult(Success s, Failure f, bool succeeded) noexcept :
        success{std::move(s)},
        failure{std::move(f)},
        succeeded{succeeded}
    {}
};

class NtCStructure {
public:
    NtCStructure();
    NtCStructure(mmdb::Manager *mmdbStru, LLKA_Structure llkaStru);
    NtCStructure(const NtCStructure &) = delete;
    NtCStructure(NtCStructure &&other) noexcept;
    ~NtCStructure();

    NtCStructure & operator=(const NtCStructure &) = delete;
    NtCStructure & operator=(NtCStructure &&other) noexcept;

    void release();

    mmdb::Manager *mmdbStru;
    LLKA_Structure llkaStru;
    bool isValid;

private:
    bool m_released;
};

class NtCSuperposition {
public:
    mmdb::Manager *mmdbStru;
    double rmsd;
};

using NtCConnectivityResult = NtCResult<NtCConnectivities, LLKA_RetCode>;
using NtCSimilarityResult = NtCResult<std::vector<NtCSimilarity>, LLKA_RetCode>;

NtCResult<LLKA_ClassifiedStep, LLKA_RetCode> ntc_classify(const NtCStructure &stru);
NtCStructure ntc_dinucleotide(mmdb::Manager *srcMmdbStru, mmdb::Residue *residue, const std::string &altconf);
NtCConnectivityResult ntc_calculate_connectivity(LLKA_NtC ntc, const NtCStructure &stru, mmdb::Manager *srcMmdbStru);
NtCSimilarityResult ntc_calculate_similarities(const NtCStructure &stru);
NtCStructure ntc_get_reference_structure(LLKA_NtC ntc);
bool ntc_initialize_classification_context(const std::string &path, std::string &error);
bool ntc_initialize_classification_context_if_needed(std::string path, std::string &error);
NtCSuperposition ntc_superpose_reference(const NtCStructure &stru, LLKA_NtC ntc);

#endif // COOT_ENABLE_NTC

#endif // _NTC_NTC_H
