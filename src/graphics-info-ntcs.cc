// vim: set sw=4 ts=4 sts=4 expandtab :

#include "ntcs.h"
#include "graphics-info.h"
#include "utils/coot-utils.hh"

int graphics_info_t::ntc_conformations_imol = -1;
int graphics_info_t::ntc_conformations_atom_index = -1;

static
mmdb::Manager * makeStepSlice(std::vector<molecule_class_info_t> &molecules, int atom_index, int imol) {
    auto &molecule = molecules[imol];

    mmdb::Atom *atom = molecule.atom_sel.atom_selection[atom_index];
    mmdb::Residue *residue = atom->residue;
    const std::string &altconf = atom->altLoc;
    mmdb::Residue *residue2 = molecule.get_following_residue(coot::residue_spec_t(residue));
    if (!residue2)
        return nullptr;

    if (!IBT::is_nucleotide(residue) || !IBT::is_nucleotide(residue2))
        return nullptr;

    // WARN: This returns a pointer that should probably be freed somewhere!!!
    mmdb::Residue *filteredResidue = coot::util::deep_copy_this_residue(residue, { true, altconf });
    mmdb::Residue *filteredResidue2 = coot::util::deep_copy_this_residue(residue2, { true, altconf });

    // WARN: Another bunch of raw pointers
    mmdb::Manager *mol = new mmdb::Manager;
	mmdb::Model *model = new mmdb::Model;
    mmdb::Chain *chain = new mmdb::Chain;

    chain->AddResidue(filteredResidue);
    chain->AddResidue(filteredResidue2);
    chain->SetChainID(residue->GetChainID());
    model->AddChain(chain);
    mol->AddModel(model);
    mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX);

    return mol;
}

static
void setupDialog(const int imol, GtkBuilder *gtkbuilder) {
    GtkWidget *dialog = widget_from_builder("ntc_conformations_dialog");
    g_object_set_data(G_OBJECT(dialog), "imol", GINT_TO_POINTER(imol));

    GtkComboBox *listOfNtCClasses = GTK_COMBO_BOX(gtk_builder_get_object(gtkbuilder, "ntc_conformations_list_of_ntc_classes"));
    if (gtk_combo_box_get_model(listOfNtCClasses) == nullptr) {
        GtkListStore *listStore = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
        GtkTreeIter iter;

        for (size_t idx = 0; idx < IBT::NTC_CLASSES.size(); idx++) {
            const auto &cls = IBT::NTC_CLASSES[idx];
            gtk_list_store_append(listStore, &iter);
            gtk_list_store_set(listStore, &iter, 0, cls.c_str(), 1, idx, -1);
        }

        gtk_combo_box_set_model(listOfNtCClasses, GTK_TREE_MODEL(listStore));
        g_object_unref(listStore);

        GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
        gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(listOfNtCClasses), renderer, TRUE);
        gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(listOfNtCClasses), renderer, "text", 0, NULL);
    }
    gtk_combo_box_set_active(listOfNtCClasses, 0);

    graphics_info_t g;
    g.ntc_conformations_setup_ntc_combobox(0);

    gtk_widget_show(dialog);
}

bool
graphics_info_t::do_ntc_conformations(int atom_index, int imol) {
    mmdb::Manager *mol = makeStepSlice(molecules, atom_index, imol);
    if (!mol) {
        std::cout << "Cannot make step\n";
        return false;
    }

    ntc_conformations_imol = imol;
    ntc_conformations_atom_index = atom_index;

    setupDialog(imol, gtkbuilder);

    ntc_conformations_generate_moving_atoms(atom_index, imol, -1);

    return true;
}

void
graphics_info_t::ntc_conformations_generate_moving_atoms(int atom_index, int imol, int ntcIdx) {
    /* We could sanitize this earlier but that would sprinkle the "ntcs.h" header
     * all over the code. I assume you have highier aims for your live rather
     * than watch Coot rebuild it self over and over whenever you touch anything */
    if (ntcIdx >= IBT::NTCS.size())
        return;

    mmdb::Manager *mol = makeStepSlice(molecules, atom_index, imol);
    if (!mol) {
        std::cout << "Cannot make step\n";
        return;
    }

    // BEWARE, BEWARE: some global data!
    imol_moving_atoms = imol;
    if (!moving_atoms_asc) {
        moving_atoms_asc = new atom_selection_container_t;
    }

    if (ntcIdx >= 0) {
        try {
            IBT::apply_NtC(mol, Geom_p(), IBT::NTCS[ntcIdx]);
        } catch (const std::runtime_error &rte) {
            std::cout << "Cannot apply NtC: " << rte.what() << "\n";
        }
    }

    *moving_atoms_asc = make_asc(mol);
    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF;

    // Draw stuff
    make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
    graphics_draw();
}

void
graphics_info_t::ntc_conformations_setup_ntc_combobox(int clsIdx) {
    GtkComboBox *listOfNtCs = GTK_COMBO_BOX(gtk_builder_get_object(gtkbuilder, "ntc_conformations_list_of_ntcs"));

    GtkTreeModel *model = gtk_combo_box_get_model(listOfNtCs);
    GtkListStore *listStore;
    if (model == nullptr) {
        listStore = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_INT);
        gtk_combo_box_set_model(listOfNtCs, GTK_TREE_MODEL(listStore));

        GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
        gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(listOfNtCs), renderer, TRUE);
        gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(listOfNtCs), renderer, "text", 0, NULL);
    } else {
        listStore = GTK_LIST_STORE(model);
        gtk_list_store_clear(listStore);
    }

    GtkTreeIter iter;

    const auto &cls = IBT::NTC_CLASSES[clsIdx];
    for (size_t idx = 0; idx < IBT::NTCS.size(); idx++) {
        const auto &ntc = IBT::NTCS[idx];
        if (ntc.NtCClass == cls) {
            gtk_list_store_append(listStore, &iter);
            gtk_list_store_set(listStore, &iter, 0, ntc.name.c_str(), 1, idx, -1);
        }
    }

    if (model == nullptr)
        g_object_unref(listStore);

    gtk_combo_box_set_active(listOfNtCs, 0);
}

#define NTC_ANG_DIST_LABEL_ACTUAL(l) "ntc_conformations_"#l"_actual"
#define NTC_SHOW_ANG_DIST_ACTUAL(mol, tag, param) \
{ \
    GtkLabel *label = GTK_LABEL(gtk_builder_get_object(gtkbuilder, NTC_ANG_DIST_LABEL_ACTUAL(tag))); \
    auto value = IBT::measure_NtC(mol, IBT::NtC::param); \
    char buf[8]; \
    snprintf(buf, 8, "%6.2f", value); \
    gtk_label_set_text(label, buf); \
}

#define NTC_ANG_DIST_LABEL_PRESCRIBED(l) "ntc_conformations_"#l"_prescribed"
#define NTC_CLEAR_ANG_DIST_PRESCRIBED(n) \
{ \
    GtkLabel *label = GTK_LABEL(gtk_builder_get_object(gtkbuilder, NTC_ANG_DIST_LABEL_PRESCRIBED(n))); \
    gtk_label_set_text(label, "-"); \
}
#define NTC_SHOW_ANG_DIST_PRESCRIBED(n, ntc) \
{ \
    GtkLabel *label = GTK_LABEL(gtk_builder_get_object(gtkbuilder, NTC_ANG_DIST_LABEL_PRESCRIBED(n))); \
    char buf[8]; \
    snprintf(buf, 8, "%6.2f", ntc.n); \
    gtk_label_set_text(label, buf); \
}

void
graphics_info_t::ntc_conformations_show_actual(int atom_index, int imol) {
    mmdb::Manager *mol = makeStepSlice(molecules, atom_index, imol);
    if (!mol)
        return;

    NTC_SHOW_ANG_DIST_ACTUAL(mol, delta_1, DELTA_1)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, epsilon_1, EPSILON_1)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, zeta_1, ZETA_1)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, alpha_2, ALPHA_2)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, beta_2, BETA_2)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, gamma_2, GAMMA_2)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, delta_2, DELTA_2)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, chi_1, CHI_1)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, chi_2, CHI_2)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, cc, CC)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, nn, NN)
    NTC_SHOW_ANG_DIST_ACTUAL(mol, mu, MU)
}

void
graphics_info_t::ntc_conformations_clear_prescribed() {
    NTC_CLEAR_ANG_DIST_PRESCRIBED(delta_1)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(epsilon_1)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(zeta_1)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(alpha_2)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(beta_2)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(gamma_2)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(delta_2)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(chi_1)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(chi_2)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(cc)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(nn)
    NTC_CLEAR_ANG_DIST_PRESCRIBED(mu)
}

void
graphics_info_t::ntc_conformations_show_prescribed(int ntcIdx) {
    if (ntcIdx >= IBT::NTCS.size()) {
        ntc_conformations_clear_prescribed();
        return;
    }

    const auto &ntc = IBT::NTCS[ntcIdx];

    NTC_SHOW_ANG_DIST_PRESCRIBED(delta_1, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(epsilon_1, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(zeta_1, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(alpha_2, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(beta_2, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(gamma_2, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(delta_2, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(chi_1, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(chi_2, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(cc, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(nn, ntc)
    NTC_SHOW_ANG_DIST_PRESCRIBED(mu, ntc)
}
