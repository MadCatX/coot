// vim: sw=4 ts=4 sts=4 expandtab

#ifndef COOT_UTIL_FILE_TIMES_HH
#define COOT_UTIL_FILE_TIMES_HH

#include <cstdint>

namespace coot {

struct FileTimes {
    FileTimes() :
        valid{false}
    {
    }

    FileTimes(uint_least64_t creation, uint_least64_t lastModification, uint_least64_t lastAccess) :
        creation{creation},
        lastModification{lastModification},
        lastAccess{lastAccess},
        valid{true}
    {
    }

    uint_least64_t creation;
    uint_least64_t lastModification;
    uint_least64_t lastAccess;

    bool valid;
};

}

#endif // COOT_UTIL_FILE_TIMES_HH
