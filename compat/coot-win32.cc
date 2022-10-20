// vim: sw=4 ts=4 sts=4 expandtab :

#include "coot-sysdep.h"

#include <windows.h>
#include <secext.h>

#include <algorithm>
#include <cwctype>
#include <array>
#include <memory>
#include <string>

#define W_PATH_BUF_SIZE 32768UL
#define W_PATH_PREFIX L"\\\\?\\"

// Let us just declare all of the internal utility functions here
static std::wstring absolute_path(const std::wstring &path);
static std::string get_error_string(DWORD error_code);
static DWORD get_file_attributes(const std::wstring &w_file_path);
static bool is_dir_attrs(DWORD attrs);
static bool is_regular_file_attrs(DWORD attrs);
static bool is_dir_internal(const std::wstring &w_file_path);
static bool is_link_internal(const std::wstring &w_file_path);
static bool is_regular_file_internal(const std::wstring &w_file_path);
static bool paths_equal(const std::wstring &first, const std::wstring &second);
static std::string user_account_info_internal(EXTENDED_NAME_FORMAT fmt);
static std::wstring windowsize_path(const std::wstring &path);
static std::wstring windowsize_path(const std::string &path, bool *ok);

// Getting Secur32.lib to link correctly has proven troublesome
// so we will just grab the library dynamically in runtime
typedef BOOLEAN SEC_ENTRY (*FuncGetUserNameExW)(EXTENDED_NAME_FORMAT NameFormat, LPWSTR lpNameBuffer, PULONG nSize);

namespace coot {
namespace sysdep {

int cpu_count() {
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);

  return sysinfo.dwNumberOfProcessors;
}

std::string current_working_dir() {
    // ret DOES include the null-terminator
    DWORD ret = GetCurrentDirectoryW(0, NULL);
    if (ret < 1) {
        return {};
    }

    std::vector<wchar_t> w_bytes(ret);
    DWORD ret2 = GetCurrentDirectoryW(ret, w_bytes.data());
    // ret2 does NOT include the null-terminator when GetCurrentDirectoryW() is called like this
    if (ret2 < 1 || ret2 != ret - 1) {
        return {};
    }

    return wide_string_to_local(w_bytes.data());
}

int create_directory(const std::string &path) {
    // Mind that we emulate POSIX mkdir() return values

    bool ok;
    std::wstring w_path = windowsize_path(path, &ok);
    if (!ok) {
        return -1;
    }

    BOOL ret = CreateDirectoryW(w_path.c_str(), NULL);
    if (ret) {
        return 0;
    }

    DWORD error = GetLastError();
    if (error == ERROR_ALREADY_EXISTS) {
        return is_dir_internal(w_path) ? 0 : -1;
    }

    return -1;
}

std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &patterns, GatherOptions options) {
    bool ok;

    std::wstring w_dir_path = windowsize_path(dir_path, &ok);
    if (!ok) {
        return {};
    }

    std::vector<std::wstring> w_patterns;
    for (const auto &p : patterns) {
        auto w_p = local_to_wide_string(p, &ok);
        if (!ok) {
            return {};
        }
        w_patterns.push_back(std::move(w_p));
    }

    std::vector<std::string> found;
    WIN32_FIND_DATAW find_data;
    for (const auto &p : w_patterns) {
        std::wstring glob = w_dir_path + L"\\" + p;

        HANDLE h = FindFirstFileW(glob.c_str(), &find_data);
        if (h == INVALID_HANDLE_VALUE)
            continue;

        do {
            std::wstring w_file_path = w_dir_path + L"\\" + std::wstring{find_data.cFileName};

            if (!(options & GATHER_LINKS) && is_link_internal(w_file_path)) {
                continue;
            }

            if (is_regular_file_attrs(find_data.dwFileAttributes) && (options & GATHER_FILES)) {
                std::string file_path = wide_string_to_local(w_file_path, &ok);
                if (ok)
                    found.push_back(std::move(file_path));
            } else if (is_dir_attrs(find_data.dwFileAttributes) && (options & GATHER_DIRECTORIES)) {
                std::string file_path = wide_string_to_local(w_file_path, &ok);
                if (ok)
                    found.push_back(std::move(file_path));
            }
        } while (FindNextFileW(h, &find_data) != 0);

        FindClose(h);
    }

    return found;
}

std::string get_fixed_font() {
    return "monospace";
}

bool is_dir(const std::string &file_path) {
    bool ok;
    std::wstring w_file_path = windowsize_path(file_path, &ok);
    if (!ok) {
        return false;
    }

    return is_dir_internal(w_file_path);
}

bool is_link(const std::string &file_path) {
    bool ok;
    std::wstring w_file_path = windowsize_path(file_path, &ok);
    if (!ok) {
        return false;
    }
    w_file_path = absolute_path(w_file_path);

    return is_link_internal(w_file_path);
}

bool is_regular_file(const std::string &file_path) {
    bool ok;
    std::wstring w_file_path = windowsize_path(file_path, &ok);
    if (!ok) {
        return false;
    }

    return is_regular_file_internal(w_file_path);
}

std::wstring local_to_wide_string(const std::string &str, bool *ok) {
    int ret = MultiByteToWideChar(
        CP_ACP,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        nullptr,
        0
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    std::vector<wchar_t> w_bytes(ret);
    ret = MultiByteToWideChar(
        CP_ACP,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        w_bytes.data(),
        ret
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    if (ok) *ok = true;
    return std::wstring(w_bytes.data());
}

bool rename(const char *old_file_path, const char *new_file_path, std::string &error_message) {
    bool ok;
    std::wstring w_old_file_path = local_to_wide_string(old_file_path, &ok);
    if (!ok) {
        error_message = "Failed to convert old_file_path to Unicode";
        return false;
    }
    std::wstring w_new_file_path = local_to_wide_string(new_file_path, &ok);
    if (!ok) {
        error_message = "Failed to convert new_file_path to Unicode";
        return false;
    }

    w_old_file_path = windowsize_path(w_old_file_path);
    w_new_file_path = windowsize_path(w_new_file_path);

    DWORD attrs = get_file_attributes(w_old_file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES) {
        error_message = "File does not exist";
        return false; // File to rename does not exist
    }

    BOOL ret = MoveFileExW(
        w_old_file_path.c_str(),
        w_new_file_path.c_str(),
        MOVEFILE_COPY_ALLOWED | MOVEFILE_REPLACE_EXISTING
    );
    if (ret) {
        return true;
    }

    error_message = get_error_string(GetLastError());
    return false;
}

void set_os_error_mode() {
    SetErrorMode(SetErrorMode(SEM_NOGPFAULTERRORBOX) | SEM_NOGPFAULTERRORBOX);
}

std::string wide_string_to_local(const std::wstring &w_str, bool *ok) {
    int ret = WideCharToMultiByte(
        CP_ACP,
        WC_NO_BEST_FIT_CHARS,
        w_str.c_str(),
        -1,
        nullptr,
        0,
        nullptr,
        nullptr
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    std::vector<char> bytes(ret);
    ret = WideCharToMultiByte(
        CP_ACP,
        WC_NO_BEST_FIT_CHARS,
        w_str.c_str(),
        -1,
        bytes.data(),
        ret,
        nullptr,
        nullptr
    );
    if (ret < 1) {
        if (ok) *ok = false;
        return {};
    }

    if (ok) *ok = true;
    return std::string(bytes.data());
}

void sleep(unsigned int secs) {
    Sleep(1000 * secs);
}

void usleep(unsigned int usecs) {
    Sleep(usecs / 1000);
}

std::string user_account_name() {
    return user_account_info_internal(NameUserPrincipal);
}

std::string user_full_name() {
    return user_account_info_internal(NameDisplay);
}

} // namespace sysdep
} // namespace coot

// Internal utility functions

static
std::wstring absolute_path(const std::wstring &path) {
    std::array<wchar_t, W_PATH_BUF_SIZE> buf{};
    LPWSTR file_name = nullptr;

    if (GetFullPathNameW(path.c_str(), buf.size(), buf.data(), &file_name) == 0)
        return {};

    return std::wstring{buf.data()};
}

static
std::string get_error_string(DWORD error_code) {
    LPSTR buf = NULL;

    if (FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        error_code,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_NEUTRAL),
        (LPSTR) &buf,
        0,
        NULL
    ))
        return buf;
    else
        return "Unknown error";
}

static
DWORD get_file_attributes(const std::wstring &w_file_path) {
    return GetFileAttributesW(w_file_path.c_str());
}

static
bool is_dir_attrs(DWORD attrs) {
    return attrs & FILE_ATTRIBUTE_DIRECTORY;
}

static
bool is_regular_file_attrs(DWORD attrs) {
    const DWORD REG_FILE_MASK = ~(FILE_ATTRIBUTE_DIRECTORY | FILE_ATTRIBUTE_HIDDEN | FILE_ATTRIBUTE_DEVICE | FILE_ATTRIBUTE_OFFLINE);
    if (attrs == FILE_ATTRIBUTE_NORMAL)
        return true;
    return attrs & REG_FILE_MASK;
}

static
bool is_dir_internal(const std::wstring &w_file_path) {
    DWORD attrs = get_file_attributes(w_file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES)
        return false;

    return attrs & FILE_ATTRIBUTE_DIRECTORY;
}

static
bool is_link_internal(const std::wstring &w_file_path) {
    HANDLE h = CreateFileW(
        w_file_path.c_str(),
        GENERIC_READ,
        FILE_SHARE_WRITE,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL
    );
    if (h == INVALID_HANDLE_VALUE) {
        return false;
    }

    std::array<wchar_t, W_PATH_BUF_SIZE> buf{};
    DWORD ret = GetFinalPathNameByHandleW(h, buf.data(), buf.size(), FILE_NAME_NORMALIZED);
    CloseHandle(h);

    if (ret == 0) {
        return false;
    }

    std::wstring w_final_path{buf.data()};
    return !paths_equal(w_file_path, w_final_path);
}

static
bool is_regular_file_internal(const std::wstring &w_file_path) {
    DWORD attrs = get_file_attributes(w_file_path);
    if (attrs == INVALID_FILE_ATTRIBUTES) {
        return false;
    }

    return is_regular_file_attrs(attrs);
}

static
bool paths_equal(const std::wstring &first, const std::wstring &second) {
    if (first.length() != second.length()) {
        return false;
    }

    for (size_t idx = 0; idx < first.length(); idx++) {
        wchar_t a = std::towlower(first[idx]);
        wchar_t b = std::towlower(second[idx]);
        if (a != b) {
            return false;
        }
    }

    return true;
}

static
std::string user_account_info_internal(EXTENDED_NAME_FORMAT fmt) {
    ULONG length = 32;
    BOOLEAN ret;

    HMODULE hSecur32 = LoadLibraryW(L"secur32.dll");
    if (hSecur32 == NULL) {
        return {};
    }
    FuncGetUserNameExW _GetUserNameExW = (FuncGetUserNameExW)GetProcAddress(hSecur32, "GetUserNameExW");
    if (_GetUserNameExW == nullptr) {
        FreeLibrary(hSecur32);
        return {};
    }

    // MSDN docn does not specify if we can pass NULL to GetUserNameExW()
    // Let us allocate a "large enough" buffer and try.
    auto buf = std::unique_ptr<wchar_t[]>(new wchar_t[length + 1]);

    ret = _GetUserNameExW(fmt, NULL, &length);
    if (!ret && GetLastError() != ERROR_MORE_DATA) {
        // We failed but not because the buffer is too small
        FreeLibrary(hSecur32);
        return {};
    }

    // Try again with larger buffer. "length" now contains the required size of the buffer.
    // None that ehavior of GetUserNameExW is inconsistent with regard to
    // whether the "length" accounts for zero-terminator.
    // Let us add 1 to be on the safe side.
    buf = std::unique_ptr<wchar_t[]>(new wchar_t[length + 1]);

    ret = _GetUserNameExW(NameUserPrincipal, NULL, &length);
    FreeLibrary(hSecur32);

    return ret ? coot::sysdep::wide_string_to_local(std::wstring{buf.get()}, nullptr) : std::string{};
}

static
std::wstring windowsize_path(const std::wstring &path) {
    std::wstring win_path{path};
    std::transform(win_path.begin(), win_path.end(), win_path.begin(),
                   [](wchar_t ch) {
                       if (ch == L'/')
                           return L'\\';
                        return ch;
                    });
    /* W_PATH_PREFIX increases the maximum path length limit to 32767 characters when used with Unicode
       variants of the respective WinAPI functions */
    if (win_path.length() > 4 && win_path.substr(0, 4) != W_PATH_PREFIX)
        return W_PATH_PREFIX + win_path;
    return win_path;
}

static
std::wstring windowsize_path(const std::string &path, bool *ok) {
    std::wstring w_path = coot::sysdep::local_to_wide_string(path, ok);
    if (!ok) {
        return L"";
    }

    return windowsize_path(w_path);
}
