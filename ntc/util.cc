// vim: set sw=4 ts=4 sts=4 expandtab :

#include "util.hh"

#include "coot-utils/coot-coord-utils.hh"

#include <mmdb2/mmdb_manager.h>

#include <cassert>
#include <cctype>
#include <cstring>

template <size_t N>
static
void _CopyCharArray(char dst[N], const char src[N]) {
    std::strncpy(dst, src, N);
}

template <typename CharArray>
static
void CopyCharArray(CharArray dst, const CharArray src) {
    _CopyCharArray<sizeof(CharArray)>(dst, src);
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
        mmdbAtom->GetLabelSeqID(),
        mmdb_altloc_to_altloc(mmdbAtom->altLoc),
        mmdbResidue->GetLabelSeqID(),
        mmdbResidue->GetInsCode(),
        0,
        &coords
    );
}

mmdb::Residue * clone_mmdb_residue(mmdb::Residue *original, const std::string &onlyAltConf) {
    mmdb::Residue *clone = coot::util::deep_copy_this_residue(original, { !onlyAltConf.empty(), onlyAltConf });
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

    const int NChains = model->GetNumberOfChains();
    for (int chainIdx = 0; chainIdx < NChains; chainIdx++) {
        mmdb::Chain *chain = model->GetChain(chainIdx);
        assert(chain);

        const int NResidues = chain->GetNumberOfResidues();
        for (int residueIdx = 0; residueIdx < NResidues; residueIdx++) {
            mmdb::Residue *residue = chain->GetResidue(residueIdx);
            assert(residue);

            const int NAtoms = residue->GetNumberOfAtoms();
            for (int atomIdx = 0; atomIdx < NAtoms; atomIdx++) {
                mmdb::Atom *atom = residue->GetAtom(atomIdx);

                LLKA_Atom llkaAtom = mmdb_atom_to_LLKA_atom(atom, residue, model->GetEntryID());
                LLKA_appendAtom(&llkaAtom, &llkaStru);

                std::cout << llkaAtom.label_atom_id << ", " << llkaAtom.type_symbol << ", " << llkaAtom.label_comp_id << ", " << llkaAtom.label_seq_id << "\n";
            }
        }
    }

    return llkaStru;
}

void relabel_mmdb_step(mmdb::Manager *relabelee, mmdb::Manager *relabeler) {
    // This function assumes that both arguments constitute a valid NtC step

    // Relabel chain
    mmdb::Chain *chainR = relabeler->GetChain(1, 0);
    mmdb::Chain *chainE = relabelee->GetChain(1, 0);
    chainE->SetChainID(chainR->GetChainID());

    // Relabel first residue
    mmdb::Residue *residueE = chainE->GetResidue(0);
    mmdb::Residue *residueR = chainR->GetResidue(0);
    residueE->SetResID(residueR->GetResName(), residueR->GetSeqNum(), residueR->GetInsCode());

    // Relabel second residue
    residueE = chainE->GetResidue(1);
    residueR = chainR->GetResidue(1);
    residueE->SetResID(residueR->GetResName(), residueR->GetSeqNum(), residueR->GetInsCode());
}
