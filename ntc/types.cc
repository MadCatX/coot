#include "types.hh"

#include <mmdb2/mmdb_manager.h>

NtCStructure::NtCStructure() :
    isValid{false},
    m_released{false}
{}

NtCStructure::NtCStructure(mmdb::Manager *mmdbStru, LLKA_Structure llkaStru) :
    mmdbStru{mmdbStru},
    llkaStru{llkaStru},
    isValid{true},
    m_released{false}
{}

NtCStructure::NtCStructure(NtCStructure &&other) noexcept :
    mmdbStru{other.mmdbStru},
    llkaStru{other.llkaStru},
    isValid{other.isValid},
    m_released{other.m_released}
{
    other.release();
}

NtCStructure::~NtCStructure() {
    if (!m_released && isValid) {
        LLKA_destroyStructure(&llkaStru);
        delete mmdbStru;
    }
}

NtCStructure & NtCStructure::operator=(NtCStructure &&other) noexcept {
    this->mmdbStru = other.mmdbStru;
    this->llkaStru = other.llkaStru;
    this->isValid = other.isValid;
    this->m_released = other.m_released;

    other.release();

    return *this;
}

void NtCStructure::release() {
    m_released = true;
}
