// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_TYPES_HH
#define _NTC_TYPES_HH

#include <LLKA/llka_connectivity_similarity.h>

#include <string>

struct NtCSimilarity {
    NtCSimilarity(LLKA_Similarity similarity, std::string NtC) noexcept :
        similarity{similarity},
        NtC{std::move(NtC)}
    {}

    LLKA_Similarity similarity;
    std::string NtC;
};

#endif // _NTC_TYPES_HH
