// vim: sw=4 ts=4 sts=4 expandtab

#ifndef COOT_UTILS_GATHER_OPTIONS_HH
#define COOT_UTILS_GATHER_OPTIONS_HH

#include <type_traits>

namespace coot {
    enum GatherOptions {
        GATHER_FILES       = 1 << 0,
        GATHER_DIRECTORIES = 1 << 1,
        GATHER_LINKS       = 1 << 2
    };
    inline GatherOptions operator|(GatherOptions lhs, GatherOptions rhs) {
        using UT = std::underlying_type<GatherOptions>::type;
        return static_cast<GatherOptions>(static_cast<UT>(lhs) | static_cast<UT>(rhs));
    }
} // namespace coot

#endif // COOT_UTILS_GATHER_OPTIONS_HH
