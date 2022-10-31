// vim: sw=4 ts=4 sts=4 expandtab

/* compat/coot-sysdep.h
 * 
 * Copyright 2008, The University of York
 * Author: Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COMPAT_COOT_SYSDEP_H
#define COMPAT_COOT_SYSDEP_H

#include <utils/gather-options.hh>
#include <utils/file-times.hh>

#include <cstdint>
#include <string>
#include <vector>
#include <type_traits>

namespace coot {
namespace sysdep {
	int cpu_count();
    std::string current_working_dir();
    int create_directory(const std::string &path);
    std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &pattern, GatherOptions options = (GATHER_FILES | GATHER_LINKS));
    bool file_exists(const std::string &file_path);
    FileTimes get_file_times(const std::string &file_path);
    int_least64_t get_file_size(const std::string &file_path);
    std::string get_fixed_font();
    std::string get_home_dir();
    bool is_dir(const std::string &file_path);
    bool is_file_writeable(const std::string &file_path);
    bool is_link(const std::string &file_path);
    bool is_regular_file(const std::string &file_path);
    bool rename_file(const std::string &old_file_path, const std::string &new_file_path, std::string &error_message);
    bool set_current_directory(const std::string &path);
    void set_os_error_mode();
    void sleep(unsigned int secs);
    void usleep(unsigned int usecs);
    std::string user_account_name();
    std::string user_full_name();
} // namespace sysdep
} // namespace coot

#if defined(COOT_BUILD_WINDOWS)
# include "coot-win32.h"
#elif defined(COOT_BUILD_POSIX)
# include "coot-posix.h"
#else
# error "Misdetected or unsupported platform"
#endif // COOT_BUILD_

#endif // COMPAT_COOT_SYSDEP_H
