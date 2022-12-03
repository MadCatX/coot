// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_TYPES_HH
#define _NTC_TYPES_HH

#include <LLKA/llka_connectivity_similarity.h>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmdb {
    class Manager;
    class Residue;
}

using NtCStepAltConf = std::pair<std::string, std::string>;
inline
bool operator==(const NtCStepAltConf &lhs, const NtCStepAltConf &rhs) noexcept {
    return (lhs.first == rhs.first) && (lhs.second == rhs.second);
}

using NtCStepAltConfs = std::vector<NtCStepAltConf>;

template <typename T>
class NtCMaybe {
public:
    bool isEmpty() const {
        return !m_hasValue;
    }

    T & value() {
        assert(m_hasValue);
        return m_storage.get();
    }

    const T & value() const {
        assert(m_hasValue);
        return m_storage.get();
    }

    operator bool() const {
        return !isEmpty();
    }

    static
    NtCMaybe empty() {
        return NtCMaybe{};
    }

    static
    NtCMaybe filled(const T &v) noexcept {
        return NtCMaybe{Storage<T, std::is_reference<T>::value>(v)};
    }

    template <typename ...Args>
    static
    NtCMaybe filled(Args &&... args) {
        return NtCMaybe{Storage<T, std::is_reference<T>::value>(std::forward<Args>(args)...)};
    }

private:
    template <typename IType, bool IsReference>
    struct Storage;

    template <typename IType>
    struct Storage<IType, false> {
        std::unique_ptr<IType> data;

        Storage() = default;
        Storage(const IType &v) noexcept
        {
            data = std::unique_ptr<IType>(new IType(v));
        }

        template <typename ...Args>
        Storage(Args &&... args)
        {
            data = std::unique_ptr<IType>(new IType(std::forward<Args>(args)...));
        }

        IType & get() { return *data.get(); }
        const IType & get() const { return *data.get(); }
    };

    template <typename IType>
    struct Storage<IType, true> {
        using CIType = const typename std::remove_cv<IType>::type;
        CIType data;

        Storage() :
            data(CIType{})
        {}

        Storage(CIType ref) noexcept :
            data{ref}
        {}

        CIType get() { return data; }
        CIType get() const { return data; }
    };

    NtCMaybe() = default;
    NtCMaybe(Storage<T, std::is_reference<T>::value> &&storage) :
        m_storage{std::move(storage)},
        m_hasValue{true}
    {}

    Storage<T, std::is_reference<T>::value> m_storage;
    bool m_hasValue;
};

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
    NtCResult(Success &&s, Failure f, bool succeeded) noexcept :
        success{std::move(s)},
        failure{std::move(f)},
        succeeded{succeeded}
    {}
};

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

struct AltConfNtCConnectivities {
    AltConfNtCConnectivities() = default;
    AltConfNtCConnectivities(std::string altconf, std::vector<NtCConnectivity> conns) noexcept :
        altconf{std::move(altconf)},
        conns{std::move(conns)}
    {}

    // Altconf of either previous or next step, depending on which connectivity was calculated
    std::string altconf;
    std::vector<NtCConnectivity> conns;
};

struct NtCConnectivities {
    NtCConnectivities() = default;
    NtCConnectivities(std::vector<AltConfNtCConnectivities> previous, std::vector<AltConfNtCConnectivities> next) noexcept :
        previous{std::move(previous)},
        next{std::move(next)}
    {}

    std::vector<AltConfNtCConnectivities> previous;
    std::vector<AltConfNtCConnectivities> next;
};
using NtCConnectivitiesResult = NtCResult<NtCConnectivities, LLKA_RetCode>;

struct NtCSimilarity {
    NtCSimilarity(LLKA_Similarity similarity, std::string NtC) noexcept :
        similarity{similarity},
        NtC{std::move(NtC)}
    {}

    LLKA_Similarity similarity;
    std::string NtC;
};
using NtCSimilarities = std::vector<NtCSimilarity>;
using NtCSimilaritiesResult = NtCResult<NtCSimilarities, LLKA_RetCode>;

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

#endif // _NTC_TYPES_HH
