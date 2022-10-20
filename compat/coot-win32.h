// vim: sw=4 ts=4 sts=4 expandtab :

#ifndef COMPAT_COOT_WIN32_H
#define COMPAT_COOT_WIN32_H

#include <string>

namespace coot {
namespace sysdep {
    std::wstring local_to_wide_string(const std::string &str, bool *ok = nullptr);
    std::string wide_string_to_local(const std::wstring &w_str, bool *ok = nullptr);
} // namespace sysdep
} // namespace coot

#endif // COMPAT_COOT_WIN32_H
