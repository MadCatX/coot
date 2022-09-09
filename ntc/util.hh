// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UTIL_HH
#define _NTC_UTIL_HH

#include <LLKA/llka_ntc.h>
#include <LLKA/llka_structure.h>

#include <string>
#include <vector>

namespace mmdb {
    class Manager;
    class Residue;
}

mmdb::Residue * clone_mmdb_residue(mmdb::Residue *original, const std::string &onlyAltConf = {});
mmdb::Manager * clone_mmdb_structure(mmdb::Manager *original);
mmdb::Manager * LLKA_structure_to_mmdb_structure(const LLKA_Structure &llkaStru);
std::vector<LLKA_NtC> make_ntc_range(LLKA_NtC first, LLKA_NtC last);
LLKA_Structure mmdb_structure_to_LLKA_structure(mmdb::Manager *mmdbStru);
void relabel_mmdb_step(mmdb::Manager *relabelee, mmdb::Manager *relabeler);

#endif // _NTC_UTIL_HH
