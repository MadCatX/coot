// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_TYPES_HH
#define _NTC_TYPES_HH

#include <LLKA/llka_connectivity_similarity.h>

#include <string>
#include <vector>

struct NtCConnectivity {
    NtCConnectivity() noexcept {
        connectivity.C5PrimeDistance = 0.0;
        connectivity.O3PrimeDistance = 0.0;
    }

    NtCConnectivity(LLKA_Connectivity connectivity, std::string NtC) noexcept :
        connectivity{connectivity},
        NtC{std::move(NtC)}
    {}

    LLKA_Connectivity connectivity;
    std::string NtC;
};

class NtCConnectivities {
public:
    NtCConnectivities() = default;
    NtCConnectivities(std::vector<NtCConnectivity> previous, std::vector<NtCConnectivity> next) noexcept :
        previous{std::move(previous)},
        next{std::move(next)}
    {}

    std::vector<NtCConnectivity> previous;
    std::vector<NtCConnectivity> next;
};

struct NtCSimilarity {
    NtCSimilarity(LLKA_Similarity similarity, std::string NtC) noexcept :
        similarity{similarity},
        NtC{std::move(NtC)}
    {}

    LLKA_Similarity similarity;
    std::string NtC;
};

#endif // _NTC_TYPES_HH
