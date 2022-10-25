// vim: sw=4 ts=4 sts=4 expandtab :

#include "coot-sysdep.h"

#include <fcntl.h>
#include <glob.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <pwd.h>

#include <array>
#include <cstring>
#include <iostream>

#define COOT_PATH_MAX 4096UL

namespace coot {
namespace sysdep {

int cpu_count() {
    return sysconf(_SC_NPROCESSORS_CONF);
}

std::string current_working_dir() {
    std::array<char, COOT_PATH_MAX> bytes{};
    if (getcwd(bytes.data(), bytes.size()))
        return std::string{bytes.data()};
    return {};
}

int create_directory(const std::string &path) {
    int istat = -1;
    struct stat s;

    int fstat = stat(path.c_str(), &s);

    // 20060411 Totally bizarre pathology!  (FC4) If I comment out the
    // following line, we fail to create the directory, presumably
    // because the S_ISDIR returns true (so we don't make the
    // directory).

    if (fstat == -1) { // file not exist
        // not exist
        std::cout << "INFO:: Creating directory " << path << std::endl;

        mode_t mode = S_IRUSR|S_IWUSR|S_IXUSR; // over-ridden
        mode = 511; // octal 777
        mode_t mode_o = umask(0);
        mode_t mkdir_mode = mode - mode_o;

        istat = mkdir(path.c_str(), mkdir_mode);
        umask(mode_o); // oh yes we do, crazy.
    } else {
        if (!S_ISDIR(s.st_mode)) {
            // exists but is not a directory
            istat = -1;
        } else {
            // was a directory already
            istat = 0; // return as if we made it
        }
    }

    return istat;
}

std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &patterns, GatherOptions options) {
    if (patterns.empty())
        return {};

    glob_t gb;
    int gb_flags = 0;
    for (const auto &p : patterns) {
        std::string path = dir_path + "/" + p;
        glob(path.c_str(), gb_flags, nullptr, &gb);
        gb_flags = GLOB_APPEND;
    }

    std::vector<std::string> found;
    size_t count = gb.gl_pathc;
    for (char **p = gb.gl_pathv; count ; p++, count--) {
        std::string sp{*p};

        if (!(options & GATHER_LINKS) && is_link(sp)) {
            continue;
        }

        if ((options & GATHER_FILES) && is_regular_file(sp)) {
            found.push_back(std::move(sp));
        } else if ((options & GATHER_DIRECTORIES) && is_dir(sp)) {
            found.push_back(std::move(sp));
        }
    }

    globfree(&gb);

    return found;
}

std::string get_fixed_font() {
    return "Sans 9";
}

std::string get_home_dir() {
    const char *s = getenv("HOME");
    if (s) {
        return std::string(s);
    } else {
        s = getenv("COOT_HOME");
        if (s)
            return std::string(s);
    }

    return {};
}

bool is_dir(const std::string &file_path) {
    struct stat buf;
    if (stat(file_path.c_str(), &buf) == -1)
        return false;

    return S_ISDIR(buf.st_mode);
}

bool is_file_writeable(const std::string &file_path) {
    struct stat buf;
    if (stat(file_path.c_str(), &buf) == -1)
        return false;

    return buf.st_mode & S_IWUSR;
}

bool is_link(const std::string &file_path) {
    struct stat buf;
    if (lstat(file_path.c_str(), &buf) == -1)
        return false;

    return S_ISLNK(buf.st_mode);
}

bool is_regular_file(const std::string &file_path) {
    struct stat buf;
    stat(file_path.c_str(), &buf);

    return S_ISREG(buf.st_mode);
}

bool rename(const char *old_file_path, const char *new_file_path, std::string &error_message) {
    int ret = ::rename(old_file_path, new_file_path);
    if (ret == 0)
        return true;

    error_message = std::string{std::strerror(ret)};

    return false;
}

bool set_current_directory(const std::string &path) {
    return chdir(path.c_str()) == 0;
}

void set_os_error_mode() {
    // NOOP
}

void sleep(unsigned int secs) {
    ::sleep(secs);
}

void usleep(unsigned int usecs) {
    ::usleep(useconds_t{usecs});
}

std::string user_account_name() {
    const char *u = getenv("USER");
    return u ? u : "";
}

std::string user_full_name() {
    std::string username = user_account_name();
    if (username.empty()) {
        return "";
    }

    passwd *pw = getpwnam(username.c_str());
    std::string full_name(pw->pw_gecos);

    return full_name.empty() ? username : full_name;
}

} // namespace sysdep
} // namespace coot
