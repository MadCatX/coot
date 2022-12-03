// vim: set sw=4 ts=4 sts=4 expandtab :

#include "util.hh"

#include "coot-utils/coot-coord-utils.hh"
#include "compat/coot-sysdep.h"
#include "coords/mmdb-crystal.h"

#include <mmdb2/mmdb_manager.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstring>

#ifdef COOT_ENABLE_WINAPI_SUSPENSION
#undef GetAtomName
#endif // COOT_ENABLE_WINAPI_SUSPENSION

#ifdef COOT_BUILD_WINDOWS
#include <Windows.h>
#endif // COOT_BUILD_WINDOWS

static
std::string trim(const char *str);

template <size_t Size>
static
void CopyCharArray(char (&dst)[Size], const char (&src)[Size]) {
    std::memcpy(dst, src, Size);
}

static
mmdb::Residue * get_standard_residue_instance(const std::string &residueName, mmdb::Manager *standardResidues) {
    int selHnd = standardResidues->NewSelection();
    standardResidues->Select(
        selHnd,
        mmdb::STYPE_RESIDUE,
        1,           //iModel
        "*",
        mmdb::ANY_RES, "*",  // starting res
        mmdb::ANY_RES, "*",  // ending res
        residueName == "Ud" ? "Ur" : residueName.c_str(), // Replace deoxyuridine with uridine because we do not have Ud as a std. template
	    "*",  // Residue must contain this atom name?
        "*",  // Residue must contain this Element?
        "*",  // altLocs
        mmdb::SKEY_NEW // selection key
    );

    mmdb::PPResidue SelResidue;
    int nSelResidues;
    standardResidues->GetSelIndex(selHnd, SelResidue, nSelResidues);

    assert(nSelResidues == 1);

    mmdb::Residue *stdResidue = coot::util::deep_copy_this_residue(SelResidue[0]);
    standardResidues->DeleteSelection(selHnd);

    return stdResidue;
}

static
void fix_up_atom_names(mmdb::Residue *relabelee, mmdb::Residue *relabeler) {
    // Change atom names back to the padded form

    // Do not assume that the atoms are in the same order
    for (int aixR = 0; aixR < relabeler->GetNumberOfAtoms(); aixR++) {
        mmdb::Atom *atomR = relabeler->GetAtom(aixR);
        std::string nameR = trim(atomR->GetAtomName());

        for (int aixE = 0; aixE < relabelee->GetNumberOfAtoms(); aixE++) {
            mmdb::Atom *atomE = relabelee->GetAtom(aixE);
            std::string nameE = trim(atomE->GetAtomName());

            if (nameR == nameE) {
                atomE->SetAtomName(atomR->GetAtomName());
            }
        }
    }
}

static
char mmdb_altloc_to_altloc(const char *mmdbAltLoc) {
    if (std::strlen(mmdbAltLoc) == 0)
        return LLKA_NO_ALTID;
    return mmdbAltLoc[0];
}

static
LLKA_Atom mmdb_atom_to_LLKA_atom(mmdb::Atom *mmdbAtom, mmdb::Residue *mmdbResidue, const char *entityId) {
    LLKA_Point coords{mmdbAtom->x, mmdbAtom->y, mmdbAtom->z};

    return LLKA_makeAtom(
        mmdbAtom->serNum,
        trim(mmdbAtom->element).c_str(),
        trim(mmdbAtom->GetAtomName()).c_str(),
        entityId,
        trim(mmdbResidue->GetLabelCompID()).c_str(),
        trim(mmdbResidue->GetLabelAsymID()).c_str(),
        nullptr,
        nullptr,
        nullptr,
        mmdbResidue->GetLabelSeqID(),
        mmdb_altloc_to_altloc(mmdbAtom->altLoc),
        mmdbResidue->GetLabelSeqID(),
        mmdbResidue->GetInsCode(),
        0,
        &coords
    );
}

static
std::string refmac_residue_name(std::string residueName) {
    residueName = trim(residueName.c_str());

    if (residueName == "A")  return "Ar";
    if (residueName == "G")  return "Gr";
    if (residueName == "T")  return "Tr";
    if (residueName == "U")  return "Ur";
    if (residueName == "C")  return "Cr";
    if (residueName == "DA") return "Ad";
    if (residueName == "DG") return "Gd";
    if (residueName == "DT") return "Td";
    if (residueName == "DU") return "Ud";
    if (residueName == "DC") return "Cd";

    std::abort();
}

static
std::string trim(const char *str) {
    size_t len = strlen(str);
    const char *first = &str[0];
    const char *last = &str[len - 1];

    size_t idx = 0;
    while (idx < len) {
        first = &str[idx++];
        if (!std::isspace(first[0]))
            break;
    }

    idx = len - 1;
    while (last > first) {
        last = &str[idx--];
        if (!std::isspace(last[0]))
            break;
    }

    return std::string(first, last - first + 1);
}

std::vector<std::string> all_altconfs(mmdb::Residue *res) {
    mmdb::Atom **atoms;
    int numAtoms;
    res->GetAtomTable(atoms, numAtoms);

    std::vector<std::string> altconfs;
    for (int idx = 0; idx < numAtoms; idx++) {
        std::string altLoc = atoms[idx]->altLoc;
        if (!altLoc.empty() && std::find(altconfs.cbegin(), altconfs.cend(), altLoc) == altconfs.cend()) {
            altconfs.push_back(altLoc);
        }
    }

    return altconfs;
}

mmdb::Residue * clone_mmdb_residue(mmdb::Residue *original, const std::string &onlyAltConf) {
    mmdb::Residue *clone = coot::util::deep_copy_this_residue(original, { !onlyAltConf.empty(), onlyAltConf });
    CopyCharArray(clone->name, original->name);
    CopyCharArray(clone->label_asym_id, original->label_asym_id);
    CopyCharArray(clone->label_comp_id, original->label_comp_id);
    CopyCharArray(clone->insCode, original->insCode);
    clone->label_seq_id = original->label_seq_id;
    clone->label_entity_id = original->label_entity_id;
    clone->SSE = original->SSE;

    return clone;
}

mmdb::Manager * clone_mmdb_structure(mmdb::Manager *original) {
    mmdb::Manager *clone = new mmdb::Manager{};

    for (int modelIdx = 1; modelIdx <= original->GetNumberOfModels(); modelIdx++) {
        mmdb::Model *originalModel = original->GetModel(modelIdx);
        mmdb::Model *clonedModel = new mmdb::Model{};

        for (int chainIdx = 0; chainIdx < originalModel->GetNumberOfChains(); chainIdx++) {
            mmdb::Chain *originalChain = originalModel->GetChain(chainIdx);
            mmdb::Chain *clonedChain = new mmdb::Chain{clonedModel, originalChain->GetChainID()};

            for (int residueIdx = 0; residueIdx < originalChain->GetNumberOfResidues(); residueIdx++) {
                mmdb::Residue *originalResidue = originalChain->GetResidue(residueIdx);
                mmdb::Residue *clonedResidue = clone_mmdb_residue(originalResidue);

                clonedChain->AddResidue(clonedResidue);
            }
        }

        clone->AddModel(clonedModel);
    }

    clone->FinishStructEdit();
    return clone;
}

mmdb::Manager * LLKA_structure_to_mmdb_structure(const LLKA_Structure &llkaStru) {
    mmdb::Model *model = new mmdb::Model{};

    for (size_t idx = 0; idx < llkaStru.nAtoms; idx++) {
        const LLKA_Atom &llkaAtom = llkaStru.atoms[idx];

        mmdb::Chain *chain = model->GetChain(llkaAtom.label_asym_id);
        if (!chain) {
            chain = new mmdb::Chain{};
            chain->SetChainID(llkaAtom.label_asym_id);
            model->AddChain(chain);
        }

        mmdb::Residue *residue = chain->GetResidue(llkaAtom.label_seq_id, llkaAtom.pdbx_PDB_ins_code);
        if (!residue) {
            residue = new mmdb::Residue{chain, llkaAtom.label_comp_id, llkaAtom.label_seq_id, llkaAtom.pdbx_PDB_ins_code};
        }

        mmdb::Atom *mmdbAtom = new mmdb::Atom{residue};
        const char altLoc[2] = { llkaAtom.label_alt_id, '\0' };
        mmdbAtom->SetAtomName(residue->GetNumberOfAtoms(), llkaAtom.id, llkaAtom.label_atom_id, altLoc, "", llkaAtom.type_symbol);
        mmdbAtom->SetCoordinates(llkaAtom.coords.x, llkaAtom.coords.y, llkaAtom.coords.z , 1.0, 1.0);
    }

    mmdb::Manager *mmdbStru = new mmdb::Manager{};
    mmdbStru->AddModel(model);
    mmdbStru->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);
    mmdbStru->FinishStructEdit();

    return mmdbStru;
}

std::vector<LLKA_NtC> make_ntc_range(LLKA_NtC first, LLKA_NtC last) {
    using UT = std::underlying_type<LLKA_NtC>::type;

    std::vector<LLKA_NtC> list{};
    UT _first = first;
    UT _last = last;
    assert(first <= last);

    for (UT ntc = _first; ntc <= _last; ntc++) {
        list.push_back(LLKA_NtC(ntc));
    }

    return list;
}

LLKA_Structure mmdb_structure_to_LLKA_structure(mmdb::Manager *mmdbStru) {
    LLKA_Structure llkaStru{};

    if (mmdbStru->GetNofBiomolecules() > 1) {
        return {};
    }

    mmdb::Model *model = mmdbStru->GetModel(1);
    if (!model) {
        return {};
    }

    mmdb::Chain **chains;
    int NChains;
    model->GetChainTable(chains, NChains);
    for (int chainIdx = 0; chainIdx < NChains; chainIdx++) {
        mmdb::Chain *chain = chains[chainIdx];

        mmdb::Residue **residues;
        int NResidues;
        chain->GetResidueTable(residues, NResidues);
        for (int residueIdx = 0; residueIdx < NResidues; residueIdx++) {
            mmdb::Residue *residue = residues[residueIdx];

            mmdb::Atom **atoms;
            int NAtoms;
            residue->GetAtomTable(atoms, NAtoms);
            for (int atomIdx = 0; atomIdx < NAtoms; atomIdx++) {
                mmdb::Atom *atom = atoms[atomIdx];

                LLKA_Atom llkaAtom = mmdb_atom_to_LLKA_atom(atom, residue, model->GetEntryID());
                LLKA_appendAtom(&llkaAtom, &llkaStru);

                std::cout << llkaAtom.label_atom_id << ", " << llkaAtom.type_symbol << ", " << llkaAtom.label_comp_id << ", " << llkaAtom.label_seq_id << ", " << (llkaAtom.label_alt_id == LLKA_NO_ALTID ? '-' : llkaAtom.label_alt_id) << "\n";
            }
        }
    }

    return llkaStru;
}

void relabel_mmdb_step(mmdb::Manager *relabelee, mmdb::Manager *relabeler, bool relabelAtomNames) {
    // This function assumes that both arguments constitute a valid NtC step
    static atom_selection_container_t standard_residues = read_standard_residues();

    // Relabel chain
    mmdb::Chain *chainR = relabeler->GetChain(1, 0);
    mmdb::Chain *chainE = relabelee->GetChain(1, 0);
    chainE->SetChainID(chainR->GetChainID());

    // Relabel first residue
    mmdb::Residue *residueE = chainE->GetResidue(0);
    mmdb::Residue *residueR = chainR->GetResidue(0);
    residueE->SetResID(residueE->GetResName(), residueR->GetSeqNum(), residueR->GetInsCode()); // Getting name from the relabelee is not a bug, we must not change it here
    if (relabelAtomNames) {
        std::string nameR = refmac_residue_name(residueE->GetResName());
        mmdb::Residue *reference = get_standard_residue_instance(nameR, standard_residues.mol);
        if (reference) {
            fix_up_atom_names(residueE, reference);
        }
    }

    // Relabel second residue
    residueE = chainE->GetResidue(1);
    residueR = chainR->GetResidue(1);
    residueE->SetResID(residueE->GetResName(), residueR->GetSeqNum(), residueR->GetInsCode()); // Getting name from the relabelee is not a bug, we must not change it here
    if (relabelAtomNames) {
        std::string nameR = refmac_residue_name(residueE->GetResName());
        mmdb::Residue *reference = get_standard_residue_instance(nameR, standard_residues.mol);
        if (reference) {
            fix_up_atom_names(residueE, reference);
        }
    }
}

void replace_bases(mmdb::Manager *replacee, mmdb::Manager *replacer) {
    static atom_selection_container_t standard_residues = read_standard_residues();

    for (int modelNo = 1; modelNo <= replacer->GetNumberOfModels(); modelNo++) {
        mmdb::Model *modelE = replacee->GetModel(modelNo);
        mmdb::Model *modelR = replacer->GetModel(modelNo);

        assert(modelE); assert(modelR);

        for (int chainIdx = 0; chainIdx < modelR->GetNumberOfChains(); chainIdx++) {
            mmdb::Chain *chainE = modelE->GetChain(chainIdx);
            mmdb::Chain *chainR = modelR->GetChain(chainIdx);

            assert(chainE); assert(chainR);

            for (int residueIdx = 0; residueIdx < chainR->GetNumberOfResidues(); residueIdx++) {
                mmdb::Residue *residueE = chainE->GetResidue(residueIdx);
                mmdb::Residue *residueR = chainR->GetResidue(residueIdx);

                assert(residueE); assert(residueR);

                std::string nameR = refmac_residue_name(residueR->GetResName());
                mmdb::Residue *standardBase = get_standard_residue_instance(nameR, standard_residues.mol);
                if (!standardBase) {
                    continue;
                }

                coot::util::mutate_base(residueE, standardBase, false);
            }
        }
    }
    replacee->FinishStructEdit();
}

#ifdef COOT_BUILD_WINDOWS
std::wstring string_to_wstring(const std::string &path) {
    int len;

    len = MultiByteToWideChar(CP_ACP, MB_ERR_INVALID_CHARS, path.c_str(), -1, nullptr, 0);
    if (len < 1) {
        return L"";
    }

    wchar_t *buf = new wchar_t[len + 1];
    len = MultiByteToWideChar(CP_ACP, MB_ERR_INVALID_CHARS, path.c_str(), -1, buf, len);
    if (len < 1) {
        return L"";
    }

    std::wstring wstr{buf};
    delete [] buf;

    return wstr;
}
#else
std::wstring string_to_wstring(const std::string &path) {
    // We do not need this right now on UNIXes
    return L"";
}
#endif // COOT_BUILD_WINDOWS
