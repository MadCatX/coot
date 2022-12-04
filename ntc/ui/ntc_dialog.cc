// vim: set sw=4 ts=4 sts=4 expandtab :

#include "ntc_dialog.hh"

#include "common.hh"
#include "conn_simil_plots_dialog.hh"
#include "util.hh"

#include <LLKA/llka_classification.h>
#include <LLKA/llka_util.h>

#include <algorithm>
#include <cassert>
#include <memory>
#include <mutex>
#include <string>

#define CONN_SIMIL_DEFAULT_WIDTH 400
#define CONN_SIMIL_DEFAULT_HEIGHT 600

static_assert(
    sizeof(std::underlying_type<LLKA_NtC>::type) == sizeof(gint),
    "LLKA_NtC is represented by a type whose size is different than size of gint. This could break GtkTreeModel "
    "because it can deal directly only with values representable by GValue and we use G_TYPE_INT to store NtC."
);

static const char * NotAvail = "- N/A -";

struct NtCDialog {
    GtkWidget *root;

    GtkLabel *delta_1_actual;
    GtkLabel *epsilon_1_actual;
    GtkLabel *zeta_1_actual;
    GtkLabel *alpha_2_actual;
    GtkLabel *beta_2_actual;
    GtkLabel *gamma_2_actual;
    GtkLabel *delta_2_actual;
    GtkLabel *chi_1_actual;
    GtkLabel *chi_2_actual;
    GtkLabel *cc_actual;
    GtkLabel *nn_actual;
    GtkLabel *mu_actual;

    GtkLabel *delta_1_diff;
    GtkLabel *epsilon_1_diff;
    GtkLabel *zeta_1_diff;
    GtkLabel *alpha_2_diff;
    GtkLabel *beta_2_diff;
    GtkLabel *gamma_2_diff;
    GtkLabel *delta_2_diff;
    GtkLabel *chi_1_diff;
    GtkLabel *chi_2_diff;
    GtkLabel *cc_diff;
    GtkLabel *nn_diff;
    GtkLabel *mu_diff;

    GtkLabel *assigned_ntc;
    GtkLabel *closest_ntc;
    GtkLabel *rmsd;
    GtkComboBox *list_of_altconfs;
    GtkComboBox *list_of_ntcs;
    LLKA_NtC closest_ntc_id;

    GtkButton *toggle_conn_simil_plots;
    NtCConnSimilPlotsDialog *conn_simil_plots_dialog;
    NtCStepAltConfs altconfs;
    GtkSignalConnection altconf_changed_sgc;
    std::shared_ptr<NtCConnectivities> connectivities;
    NtCSimilarities similarities;

    GtkSignalConnection list_of_ntcs_changed_sgc;

    NtCDialogOptions options;

    bool destroyed;
};

NtCDialogOptions::NtCDialogOptions() :
    onAltconfChanged{nullptr},
    onDisplayedNtCChanged{nullptr},
    onAccepted{nullptr},
    onRejected{nullptr},
    connSimilDlgWidth{CONN_SIMIL_DEFAULT_WIDTH},
    connSimilDlgHeight{CONN_SIMIL_DEFAULT_HEIGHT}
{
}

static
void destroy_ui(NtCDialog *dlg);

static
LLKA_NtC get_selected_ntc(NtCDialog *dlg);

static
void switch_list_to_ntc(GtkComboBox *list_of_ntcs, LLKA_NtC targetNtC);

static
void display_connectivities(NtCDialog *dlg, LLKA_NtC ntc) {
    ntc_csp_dialog_update_connectivities(
        dlg->conn_simil_plots_dialog,
        dlg->connectivities,
        ntc
    );
}

static
void on_ntc_dialog_closed(GtkWidget *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog*>(data);

    // Do not call destroy_ui() because our root widget is already gone
    ntc_csp_dialog_destroy(dlg->conn_simil_plots_dialog);
    dlg->destroyed = true;

    if (dlg->options.onRejected)
        dlg->options.onRejected(dlg, dlg->closest_ntc_id);
}

static
void on_cancel_button_clicked(GtkButton *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog*>(data);

    destroy_ui(dlg);

    if (dlg->options.onRejected)
        dlg->options.onRejected(dlg, dlg->closest_ntc_id);
}

static
void on_displayed_ntc_changed(GtkComboBox *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog*>(data);
}

static
void on_ok_button_clicked(GtkButton *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog*>(data);

    if (dlg->options.onAccepted)
        dlg->options.onAccepted(dlg, get_selected_ntc(dlg));

    destroy_ui(dlg);
}

static
void on_list_of_altconfs_changed(GtkComboBox *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog *>(data);
    GtkListStore *store = GTK_LIST_STORE(gtk_combo_box_get_model(self));
    assert(store);

    GtkTreeIter iter;
    int idx;
    if (gtk_combo_box_get_active_iter(self, &iter) == TRUE) {
        gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &idx, -1);
        assert(idx >= 0 && idx < dlg->altconfs.size());

        if (dlg->options.onAltconfChanged) {
            dlg->options.onAltconfChanged(dlg, dlg->altconfs[idx]);
        }
    }
}

static
void on_list_of_ntcs_changed(GtkComboBox *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog*>(data);
    GtkListStore *store = GTK_LIST_STORE(gtk_combo_box_get_model(self));
    assert(store);

    GtkTreeIter iter;
    int ntc;
    if (gtk_combo_box_get_active_iter(self, &iter) == TRUE) {
        gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &ntc, -1);

        if (dlg->options.onDisplayedNtCChanged)
            dlg->options.onDisplayedNtCChanged(dlg, LLKA_NtC(ntc));
    }
}

static
void on_ntc_csp_dialog_destroyed(NtCConnSimilPlotsDialog *self, gpointer data)
{
    NtCDialog *dlg = static_cast<NtCDialog *>(data);
    dlg->conn_simil_plots_dialog = nullptr;
}

static
void on_reset_ntc_clicked(GtkButton *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog *>(data);
    switch_list_to_ntc(dlg->list_of_ntcs, dlg->closest_ntc_id);
}

static
void on_similarity_selected(NtCDialog *dlg, NtCSimilarity similarity) {
    LLKA_NtC ntc = LLKA_nameToNtC(similarity.NtC.c_str());
    if (ntc != LLKA_INVALID_NTC) {
        switch_list_to_ntc(dlg->list_of_ntcs, ntc);
    }
}

static
void on_toggle_conn_simil_plots_clicked(GtkButton *self, gpointer data) {
    NtCDialog *dlg = static_cast<NtCDialog *>(data);
    NtCConnSimilPlotsDialog *cspDlg = dlg->conn_simil_plots_dialog;

    if (cspDlg) {
        if (ntc_csp_dialog_is_valid(cspDlg)) {
            ntc_csp_dialog_destroy(cspDlg);
        }
        dlg->conn_simil_plots_dialog = nullptr;
    } else {
        LLKA_NtC ntc = get_selected_ntc(dlg);
        assert(ntc != LLKA_INVALID_NTC);

        dlg->conn_simil_plots_dialog = ntc_csp_dialog_make(
            [dlg](NtCSimilarity similarity) { on_similarity_selected(dlg, similarity); },
            [dlg](int width, int height) { dlg->options.connSimilDlgWidth = width; dlg->options.connSimilDlgHeight = height; }
        );
        assert(dlg->conn_simil_plots_dialog);

        g_signal_connect(
            GTK_WIDGET(ntc_csp_dialog_widget(dlg->conn_simil_plots_dialog)),
            "destroy",
            G_CALLBACK(on_ntc_csp_dialog_destroyed),
            dlg
        );

        ntc_csp_dialog_update_connectivities(
            dlg->conn_simil_plots_dialog,
            dlg->connectivities,
            ntc
        );
        ntc_csp_dialog_update_similarities(dlg->conn_simil_plots_dialog, dlg->similarities);
        ntc_csp_dialog_show(dlg->conn_simil_plots_dialog, dlg->options.connSimilDlgWidth, dlg->options.connSimilDlgHeight);
    }
}

static
double convert_angle(double ang) {
    return LLKA_fullAngleFromDeg(LLKA_rad2deg(ang));
}

static
void destroy_ui(NtCDialog *dlg) {
    ntc_csp_dialog_destroy(dlg->conn_simil_plots_dialog);
    dlg->conn_simil_plots_dialog = nullptr;
    gtk_widget_destroy(dlg->root);

    dlg->destroyed = true;
}

static
std::string format_decimal_number(double value, unsigned int numDigits, unsigned int numDecimals) {
    size_t len = snprintf(nullptr, 0, "%*.*f", numDigits, numDecimals, value);
    char *buf = new char[len + 1];

    snprintf(buf, len + 1, "%*.*f", numDigits, numDecimals, value);

    std::string str{buf};
    delete [] buf;
    return str;
}

static
LLKA_NtC get_selected_ntc(NtCDialog *dlg) {
    LLKA_NtC ntc = LLKA_INVALID_NTC;
    GtkTreeModel *store = gtk_combo_box_get_model(dlg->list_of_ntcs);
    GtkTreeIter iter;

    gtk_combo_box_get_active_iter(dlg->list_of_ntcs, &iter);
    gtk_tree_model_get(store, &iter, 0, &ntc, -1);

    return ntc;
}

static
void fill_list_of_altconfs(GtkListStore *store, const NtCStepAltConfs &altconfs, GtkSignalConnection sgc) {
    assert(store);

    GtkTreeIter iter;

    sgc.block();

    gtk_list_store_clear(store);

    for (size_t idx = 0; idx < altconfs.size(); idx++) {
        const NtCStepAltConf &ac = altconfs[idx];
        std::string text = (ac.first.empty() ? "(no altconf)" : ac.first) + " | " + (ac.second.empty() ? "(no altconf)" : ac.second);

        gtk_list_store_append(store, &iter);
        gtk_list_store_set(
            store, &iter,
            0, int(idx),
            1, text.c_str(),
            -1
        );
    }

    sgc.unblock();
}

static
void prepare_list_of_altconfs(GtkComboBox *combo, const NtCStepAltConfs &altconfs) {
    GtkListStore *store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);
    gtk_combo_box_set_model(combo, GTK_TREE_MODEL(store));

    GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, TRUE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", 1, NULL);

    if (altconfs.size() > 0) {
        fill_list_of_altconfs(store, altconfs, GtkSignalConnection{});
        gtk_combo_box_set_active(combo, 0);
    }
}

static
void prepare_list_of_ntcs(GtkComboBox *combo) {
    GtkListStore *store = gtk_list_store_new(2, G_TYPE_INT, G_TYPE_STRING);
    GtkTreeIter iter;

    using UT = std::underlying_type<LLKA_NtC>::type;
    for (UT ntc = 0; ntc <= LLKA_LAST_NTC; ntc++) {
        gtk_list_store_append(store, &iter);
        gtk_list_store_set(
            store, &iter,
            0, ntc,
            1, LLKA_NtCToName(LLKA_NtC(ntc)),
            -1
        );
    }

    gtk_combo_box_set_model(combo, GTK_TREE_MODEL(store));

    GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, TRUE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", 1, NULL);

    gtk_combo_box_set_active(combo, 0);
}

static
void set_metrics_label_text(GtkLabel *label, double value) {
    gtk_label_set_text(label, format_decimal_number(value, 5, 2).c_str());
}

static
void set_ntc_text(GtkLabel *label, LLKA_NtC ntc) {
    gtk_label_set_text(label, LLKA_NtCToName(ntc));
}

static
void switch_list_to_ntc(GtkComboBox *list_of_ntcs, LLKA_NtC targetNtC) {
    GtkListStore *store = GTK_LIST_STORE(gtk_combo_box_get_model(list_of_ntcs));
    assert(store);

    LLKA_NtC ntc;
    GtkTreeIter iter;
    gboolean hasNext = gtk_tree_model_get_iter_first(GTK_TREE_MODEL(store), &iter);
    gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &ntc, -1);
    while (ntc != targetNtC && hasNext == TRUE) {
        hasNext = gtk_tree_model_iter_next(GTK_TREE_MODEL(store), &iter);
        gtk_tree_model_get(GTK_TREE_MODEL(store), &iter, 0, &ntc, -1);
    }

    gtk_combo_box_set_active_iter(list_of_ntcs, &iter);
}

void ntc_dialog_change_ntc(NtCDialog *dlg, LLKA_NtC ntc) {
    assert(!dlg->destroyed);

    dlg->list_of_ntcs_changed_sgc.block();
    switch_list_to_ntc(dlg->list_of_ntcs, ntc);
    dlg->list_of_ntcs_changed_sgc.unblock();
}

void ntc_dialog_destroy(NtCDialog *dlg) {
    if (!dlg)
        return;

    if (!dlg->destroyed) {
        destroy_ui(dlg);
    }

    delete dlg;
}

void ntc_dialog_display_classification(NtCDialog *dlg, const NtCMaybe<LLKA_ClassifiedStep> &_classified) {
    assert(!dlg->destroyed);

    if (_classified) {
        const auto &classified = _classified.value();

        dlg->closest_ntc_id = classified.closestNtC;

        set_ntc_text(dlg->assigned_ntc, classified.assignedNtC);
        set_ntc_text(dlg->closest_ntc, classified.closestNtC);

        gtk_label_set_text(dlg->rmsd, format_decimal_number(classified.rmsdToClosestNtC, 5, 2).c_str());

        set_metrics_label_text(dlg->delta_1_actual, convert_angle(classified.metrics.delta_1));
        set_metrics_label_text(dlg->epsilon_1_actual, convert_angle(classified.metrics.epsilon_1));
        set_metrics_label_text(dlg->zeta_1_actual, convert_angle(classified.metrics.zeta_1));
        set_metrics_label_text(dlg->alpha_2_actual, convert_angle(classified.metrics.alpha_2));
        set_metrics_label_text(dlg->beta_2_actual, convert_angle(classified.metrics.beta_2));
        set_metrics_label_text(dlg->gamma_2_actual, convert_angle(classified.metrics.gamma_2));
        set_metrics_label_text(dlg->delta_2_actual, convert_angle(classified.metrics.delta_2));
        set_metrics_label_text(dlg->chi_1_actual, convert_angle(classified.metrics.chi_1));
        set_metrics_label_text(dlg->chi_2_actual, convert_angle(classified.metrics.chi_2));
        set_metrics_label_text(dlg->cc_actual, classified.metrics.CC);
        set_metrics_label_text(dlg->nn_actual, classified.metrics.NN);
        set_metrics_label_text(dlg->mu_actual, convert_angle(classified.metrics.mu));

        set_metrics_label_text(dlg->delta_1_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.delta_1));
        set_metrics_label_text(dlg->epsilon_1_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.epsilon_1));
        set_metrics_label_text(dlg->zeta_1_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.zeta_1));
        set_metrics_label_text(dlg->alpha_2_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.alpha_2));
        set_metrics_label_text(dlg->beta_2_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.beta_2));
        set_metrics_label_text(dlg->gamma_2_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.gamma_2));
        set_metrics_label_text(dlg->delta_2_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.delta_2));
        set_metrics_label_text(dlg->chi_1_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.chi_1));
        set_metrics_label_text(dlg->chi_2_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.chi_2));
        set_metrics_label_text(dlg->cc_diff, classified.differencesFromNtCAverages.CC);
        set_metrics_label_text(dlg->nn_diff, classified.differencesFromNtCAverages.NN);
        set_metrics_label_text(dlg->mu_diff, LLKA_rad2deg(classified.differencesFromNtCAverages.mu));

        switch_list_to_ntc(dlg->list_of_ntcs, dlg->closest_ntc_id);
    } else {
        gtk_label_set_text(dlg->assigned_ntc, NotAvail);
        gtk_label_set_text(dlg->closest_ntc, NotAvail);

        gtk_label_set_text(dlg->rmsd, NotAvail);

        gtk_label_set_text(dlg->delta_1_actual, NotAvail);
        gtk_label_set_text(dlg->epsilon_1_actual, NotAvail);
        gtk_label_set_text(dlg->zeta_1_actual, NotAvail);
        gtk_label_set_text(dlg->alpha_2_actual, NotAvail);
        gtk_label_set_text(dlg->beta_2_actual, NotAvail);
        gtk_label_set_text(dlg->gamma_2_actual, NotAvail);
        gtk_label_set_text(dlg->delta_2_actual, NotAvail);
        gtk_label_set_text(dlg->chi_1_actual, NotAvail);
        gtk_label_set_text(dlg->chi_2_actual, NotAvail);
        gtk_label_set_text(dlg->cc_actual, NotAvail);
        gtk_label_set_text(dlg->nn_actual, NotAvail);
        gtk_label_set_text(dlg->mu_actual, NotAvail);

        gtk_label_set_text(dlg->delta_1_diff, NotAvail);
        gtk_label_set_text(dlg->epsilon_1_diff, NotAvail);
        gtk_label_set_text(dlg->zeta_1_diff, NotAvail);
        gtk_label_set_text(dlg->alpha_2_diff, NotAvail);
        gtk_label_set_text(dlg->beta_2_diff, NotAvail);
        gtk_label_set_text(dlg->gamma_2_diff, NotAvail);
        gtk_label_set_text(dlg->delta_2_diff, NotAvail);
        gtk_label_set_text(dlg->chi_1_diff, NotAvail);
        gtk_label_set_text(dlg->chi_2_diff, NotAvail);
        gtk_label_set_text(dlg->cc_diff, NotAvail);
        gtk_label_set_text(dlg->nn_diff, NotAvail);
        gtk_label_set_text(dlg->mu_diff, NotAvail);
    }
}

void ntc_dialog_display_differences(NtCDialog *dlg, const NtCMaybe<LLKA_StepMetrics> &_differences) {
    if (_differences) {
        const auto &differences = _differences.value();

        set_metrics_label_text(dlg->delta_1_diff, LLKA_rad2deg(differences.delta_1));
        set_metrics_label_text(dlg->epsilon_1_diff, LLKA_rad2deg(differences.epsilon_1));
        set_metrics_label_text(dlg->zeta_1_diff, LLKA_rad2deg(differences.zeta_1));
        set_metrics_label_text(dlg->alpha_2_diff, LLKA_rad2deg(differences.alpha_2));
        set_metrics_label_text(dlg->beta_2_diff, LLKA_rad2deg(differences.beta_2));
        set_metrics_label_text(dlg->gamma_2_diff, LLKA_rad2deg(differences.gamma_2));
        set_metrics_label_text(dlg->delta_2_diff, LLKA_rad2deg(differences.delta_2));
        set_metrics_label_text(dlg->chi_1_diff, LLKA_rad2deg(differences.chi_1));
        set_metrics_label_text(dlg->chi_2_diff, LLKA_rad2deg(differences.chi_2));
        set_metrics_label_text(dlg->cc_diff, differences.CC);
        set_metrics_label_text(dlg->nn_diff, differences.NN);
        set_metrics_label_text(dlg->mu_diff, LLKA_rad2deg(differences.mu));
    } else {
        gtk_label_set_text(dlg->delta_1_diff, NotAvail);
        gtk_label_set_text(dlg->epsilon_1_diff, NotAvail);
        gtk_label_set_text(dlg->zeta_1_diff, NotAvail);
        gtk_label_set_text(dlg->alpha_2_diff, NotAvail);
        gtk_label_set_text(dlg->beta_2_diff, NotAvail);
        gtk_label_set_text(dlg->gamma_2_diff, NotAvail);
        gtk_label_set_text(dlg->delta_2_diff, NotAvail);
        gtk_label_set_text(dlg->chi_1_diff, NotAvail);
        gtk_label_set_text(dlg->chi_2_diff, NotAvail);
        gtk_label_set_text(dlg->cc_diff, NotAvail);
        gtk_label_set_text(dlg->nn_diff, NotAvail);
    }
}

void ntc_dialog_display_rmsd(NtCDialog *dlg, const NtCMaybe<double> &_rmsd) {
    if (_rmsd) {
        double rmsd = _rmsd.value();
        gtk_label_set_text(dlg->rmsd, format_decimal_number(rmsd, 5, 2).c_str());
    } else {
        gtk_label_set_text(dlg->rmsd, NotAvail);
    }
}

NtCStepAltConf ntc_dialog_get_current_step_altconf(NtCDialog *dlg) {
    GtkTreeModel *store = GTK_TREE_MODEL(gtk_combo_box_get_model(dlg->list_of_altconfs));
    assert(store);

    GtkTreeIter iter;
    if (gtk_combo_box_get_active_iter(dlg->list_of_altconfs, &iter) == FALSE) {
        assert(false);
    }

    int idx;
    gtk_tree_model_get(store, &iter, 0, &idx, -1);
    assert(idx >= 0 && idx < dlg->altconfs.size());

    return dlg->altconfs[idx];
}

LLKA_NtC ntc_dialog_get_current_ntc(NtCDialog *dlg) {
    GtkTreeModel *store = GTK_TREE_MODEL(gtk_combo_box_get_model(dlg->list_of_ntcs));
    assert(store);

    GtkTreeIter iter;
    if (gtk_combo_box_get_active_iter(dlg->list_of_ntcs, &iter) == FALSE) {
        assert(false);
    }

    LLKA_NtC ntc;
    gtk_tree_model_get(store, &iter, 0, &ntc, -1);
    assert(ntc != LLKA_INVALID_NTC);

    return ntc;
}

NtCDialogOptions ntc_dialog_get_options(NtCDialog *dlg) {
    assert(dlg);

    return dlg->options;
}

bool ntc_dialog_is_valid(NtCDialog *dlg) {
    return !dlg->destroyed;
}

NtCDialog * ntc_dialog_make(const NtCDialogOptions &options) {
    GtkBuilder *b = gtk_builder_new_from_file(get_glade_file("ntc_dialog.glade").c_str());
    assert(b);

    NtCDialog *dialog = new NtCDialog{};

    dialog->root = GTK_WIDGET(get_widget(b, "dialog"));
    assert(dialog->root);
    dialog->delta_1_actual = GTK_LABEL(get_widget(b, "delta_1_actual"));
    assert(dialog->delta_1_actual);
    dialog->epsilon_1_actual = GTK_LABEL(get_widget(b, "epsilon_1_actual"));
    assert(dialog->epsilon_1_actual);
    dialog->zeta_1_actual = GTK_LABEL(get_widget(b, "zeta_1_actual"));
    assert(dialog->zeta_1_actual);
    dialog->alpha_2_actual = GTK_LABEL(get_widget(b, "alpha_2_actual"));
    assert(dialog->alpha_2_actual);
    dialog->beta_2_actual = GTK_LABEL(get_widget(b, "beta_2_actual"));
    assert(dialog->beta_2_actual);
    dialog->gamma_2_actual = GTK_LABEL(get_widget(b, "gamma_2_actual"));
    assert(dialog->gamma_2_actual);
    dialog->delta_2_actual = GTK_LABEL(get_widget(b, "delta_2_actual"));
    assert(dialog->delta_2_actual);
    dialog->chi_1_actual = GTK_LABEL(get_widget(b, "chi_1_actual"));
    assert(dialog->chi_1_actual);
    dialog->chi_2_actual = GTK_LABEL(get_widget(b, "chi_2_actual"));
    assert(dialog->chi_2_actual);
    dialog->cc_actual = GTK_LABEL(get_widget(b, "cc_actual"));
    assert(dialog->cc_actual);
    dialog->nn_actual = GTK_LABEL(get_widget(b, "nn_actual"));
    assert(dialog->nn_actual);
    dialog->mu_actual = GTK_LABEL(get_widget(b, "mu_actual"));
    assert(dialog->mu_actual);

    dialog->delta_1_diff = GTK_LABEL(get_widget(b, "delta_1_diff"));
    assert(dialog->delta_1_diff);
    dialog->epsilon_1_diff = GTK_LABEL(get_widget(b, "epsilon_1_diff"));
    assert(dialog->epsilon_1_diff);
    dialog->zeta_1_diff = GTK_LABEL(get_widget(b, "zeta_1_diff"));
    assert(dialog->zeta_1_diff);
    dialog->alpha_2_diff = GTK_LABEL(get_widget(b, "alpha_2_diff"));
    assert(dialog->alpha_2_diff);
    dialog->beta_2_diff = GTK_LABEL(get_widget(b, "beta_2_diff"));
    assert(dialog->beta_2_diff);
    dialog->gamma_2_diff = GTK_LABEL(get_widget(b, "gamma_2_diff"));
    assert(dialog->gamma_2_diff);
    dialog->delta_2_diff = GTK_LABEL(get_widget(b, "delta_2_diff"));
    assert(dialog->delta_2_diff);
    dialog->chi_1_diff = GTK_LABEL(get_widget(b, "chi_1_diff"));
    assert(dialog->chi_1_diff);
    dialog->chi_2_diff = GTK_LABEL(get_widget(b, "chi_2_diff"));
    assert(dialog->chi_2_diff);
    dialog->cc_diff = GTK_LABEL(get_widget(b, "cc_diff"));
    assert(dialog->cc_diff);
    dialog->nn_diff = GTK_LABEL(get_widget(b, "nn_diff"));
    assert(dialog->nn_diff);
    dialog->mu_diff = GTK_LABEL(get_widget(b, "mu_diff"));
    assert(dialog->mu_diff);

    dialog->assigned_ntc = GTK_LABEL(get_widget(b, "assigned_ntc"));
    assert(dialog->assigned_ntc);
    dialog->closest_ntc = GTK_LABEL(get_widget(b, "closest_ntc"));
    assert(dialog->closest_ntc);
    dialog->rmsd = GTK_LABEL(get_widget(b, "rmsd"));
    assert(dialog->rmsd);

    dialog->list_of_altconfs = GTK_COMBO_BOX(get_widget(b, "list_of_altconfs"));
    assert(dialog->list_of_altconfs);
    dialog->altconfs = {};
    prepare_list_of_altconfs(dialog->list_of_altconfs, dialog->altconfs);
    dialog->altconf_changed_sgc = GtkSignalConnection{dialog->list_of_altconfs, "changed", G_CALLBACK(on_list_of_altconfs_changed), dialog};

    dialog->list_of_ntcs = GTK_COMBO_BOX(get_widget(b, "list_of_ntcs"));
    assert(dialog->list_of_ntcs);
    prepare_list_of_ntcs(dialog->list_of_ntcs);
    dialog->list_of_ntcs_changed_sgc = GtkSignalConnection{dialog->list_of_ntcs, "changed", G_CALLBACK(on_list_of_ntcs_changed), dialog};

    dialog->toggle_conn_simil_plots = GTK_BUTTON(get_widget(b, "toggle_conn_simil_plots"));
    assert(dialog->toggle_conn_simil_plots);
    dialog->conn_simil_plots_dialog = nullptr;
    g_signal_connect(dialog->toggle_conn_simil_plots, "clicked", G_CALLBACK(on_toggle_conn_simil_plots_clicked), dialog);
    dialog->closest_ntc_id = LLKA_INVALID_NTC;

    dialog->connectivities = std::shared_ptr<NtCConnectivities>(new NtCConnectivities{});

    dialog->options = options;

    dialog->destroyed = false;

    GtkButton *resetNtC = GTK_BUTTON(get_widget(b, "reset_displayed_ntc"));
    assert(resetNtC);
    g_signal_connect(resetNtC, "clicked", G_CALLBACK(on_reset_ntc_clicked), dialog);

    GtkButton *cancel= GTK_BUTTON(get_widget(b, "cancel_button"));
    assert(cancel);
    g_signal_connect(cancel, "clicked", G_CALLBACK(on_cancel_button_clicked), dialog);

    GtkButton *ok = GTK_BUTTON(get_widget(b, "ok_button"));
    assert(ok);
    g_signal_connect(ok, "clicked", G_CALLBACK(on_ok_button_clicked), dialog);

    g_signal_connect(dialog->root, "destroy", G_CALLBACK(on_ntc_dialog_closed), dialog);

    g_object_unref(b);

    return dialog;
}

void ntc_dialog_show(NtCDialog *dlg) {
    assert(!dlg->destroyed);

    gtk_widget_show(dlg->root);
}

void ntc_dialog_update_connectivities(NtCDialog *dlg, LLKA_NtC ntc, NtCConnectivities connectivities) {
    assert(!dlg->destroyed);

    *dlg->connectivities = std::move(connectivities);
    if (dlg->conn_simil_plots_dialog) {
        display_connectivities(dlg, ntc);
    }
}

void ntc_dialog_update_similarities(NtCDialog *dlg, std::vector<NtCSimilarity> similarities) {
    assert(!dlg->destroyed);

    dlg->similarities = std::move(similarities);

    if (dlg->conn_simil_plots_dialog) {
        ntc_csp_dialog_update_similarities(dlg->conn_simil_plots_dialog, dlg->similarities);
    }
}

void ntc_dialog_update_step_altconfs(NtCDialog *dlg, const NtCStepAltConfs &altconfs) {
    assert(!dlg->destroyed);
    assert(altconfs.size() > 0);

    dlg->altconfs = altconfs;
    fill_list_of_altconfs(GTK_LIST_STORE(gtk_combo_box_get_model(dlg->list_of_altconfs)), dlg->altconfs, dlg->altconf_changed_sgc);
    gtk_combo_box_set_active(dlg->list_of_altconfs, 0);
}
