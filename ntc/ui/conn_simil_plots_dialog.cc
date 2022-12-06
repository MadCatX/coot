// vim: set sw=4 ts=4 sts=4 expandtab :

#include "conn_simil_plots_dialog.hh"

#include "common.hh"
#include "plot_widget.hh"
#include "util.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

struct NtCConnSimilPlotsDialog {
    GtkDialog *root;

    GtkLabel *connectivity_previous_caption;
    GtkLabel *connectivity_next_caption;

    PlotWidget *connectivity_previous_plot;
    PlotWidget *connectivity_next_plot;
    PlotWidget *similarity_plot;

    GtkComboBox *previous_altconfs;
    GtkComboBox *next_altconfs;
    GtkSignalConnection previous_altconf_changed_sgc;
    GtkSignalConnection next_altconf_changed_sgc;

    OnSimilaritySelected on_similarity_point_selected;
    OnWidgetResized on_widget_resized;

    std::shared_ptr<NtCConnectivities> connectivities;

    bool destroyed;
};

static const gchar *PrevAltconfTooltipText = "Alternate configuration of the previous step if the previous step has any atoms in alternate positions. "
                                             "Alternate configuration of the current step is set in the main NtC window.";
static const gchar *NextAltconfTooltipText = "Alternate configuration of the next step if the next step has any atoms in alternate positions. "
                                             "Alternate configuration of the current step is set in the main NtC window.";

static
double calculate_axis_maximum(const std::vector<DotPlotPoint> &points, double DotPlotPoint::* coord, double scaleFactor);

static
std::string get_current_altconf(GtkComboBox *combo);

static
void display_connectivities(NtCConnSimilPlotsDialog *dlg) {
    auto previousAltconf = get_current_altconf(dlg->previous_altconfs);
    auto prevIt = std::find_if(
        dlg->connectivities->previous.cbegin(),
        dlg->connectivities->previous.cend(),
        [&previousAltconf](const AltConfNtCConnectivities &conn) { return previousAltconf == conn.altconf; }
    );

    auto nextAltconf = get_current_altconf(dlg->next_altconfs);
    auto nextIt = std::find_if(
        dlg->connectivities->next.cbegin(),
        dlg->connectivities->next.cend(),
        [&nextAltconf](const AltConfNtCConnectivities &conn) { return nextAltconf == conn.altconf; }
    );

    std::vector<DotPlotPoint> prevPoints;
    std::vector<DotPlotPoint> nextPoints;

    if (prevIt != dlg->connectivities->previous.cend()) {
        for (const auto &prev : prevIt->conns) {
            prevPoints.emplace_back(prev.connectivity.C5PrimeDistance, prev.connectivity.O3PrimeDistance, 1.0, 1.0, 0.0, prev.NtC);
        }
    }

    if (nextIt != dlg->connectivities->next.cend()) {
        for (const auto &next : nextIt->conns) {
            nextPoints.emplace_back(next.connectivity.C5PrimeDistance, next.connectivity.O3PrimeDistance, 0.0, 1.0, 1.0, next.NtC);
        }
    }

    ntc_plot_widget_set_points(dlg->connectivity_previous_plot, std::move(prevPoints));
    ntc_plot_widget_set_points(dlg->connectivity_next_plot, std::move(nextPoints));

    ntc_plot_widget_reset_zoom(dlg->connectivity_previous_plot);
    ntc_plot_widget_reset_zoom(dlg->connectivity_next_plot);
}

static
void fill_altconfs_list(GtkListStore *store, const std::vector<std::string> &altconfs, GtkSignalConnection sgc) {
    assert(store);
    assert(!altconfs.empty());

    // Block the "changed" signal so that we never see an empty box
    sgc.block();

    gtk_list_store_clear(store);

    GtkTreeIter iter;
    for (size_t idx = 0; idx < altconfs.size(); idx++) {
        const auto &ac = altconfs[idx];
        std::string text = ac.empty() ? "(no altconfs) " : ac;

        gtk_list_store_append(store, &iter);
        gtk_list_store_set(
            store, &iter,
            0, int(idx),
            1, text.c_str(),
            2, ac.c_str(),
            -1
        );
    }

    sgc.unblock();
}

static
std::vector<std::string> gather_altconfs(const std::vector<AltConfNtCConnectivities> &conns) {
    if (conns.empty()) {
        return { "" };
    } else {
        std::vector<std::string> altconfs;
        std::transform(conns.cbegin(), conns.cend(), std::back_inserter(altconfs), [](const AltConfNtCConnectivities &ac) { return ac.altconf; });

        return altconfs;
    }
}

static
std::string get_current_altconf(GtkComboBox *combo) {
    GtkTreeModel *store = GTK_TREE_MODEL(gtk_combo_box_get_model(combo));
    assert(store);

    GtkTreeIter iter;
    if (gtk_combo_box_get_active_iter(combo, &iter) == FALSE) {
        assert(false);
    }

    gchar *str;
    gtk_tree_model_get(store, &iter, 2, &str, -1);
    std::string altconf{str};
    g_free(str);

    return altconf;
}

static
void on_altconf_changed(GtkComboBox *, gpointer data) {
    NtCConnSimilPlotsDialog *dlg = static_cast<NtCConnSimilPlotsDialog *>(data);

    display_connectivities(dlg);
}

static
void on_dialog_closed(GtkWidget *self, gpointer data) {
    NtCConnSimilPlotsDialog *dlg = static_cast<NtCConnSimilPlotsDialog *>(data);
    dlg->destroyed = true;
}

static
void on_close_button_clicked(GtkButton *self, gpointer data) {
    NtCConnSimilPlotsDialog *dlg = static_cast<NtCConnSimilPlotsDialog *>(data);
    dlg->destroyed = true;

    gtk_widget_destroy(GTK_WIDGET(dlg->root));
}

static
gboolean on_configure_event(GtkWidget *self, GdkEvent *event, gpointer data) {
    NtCConnSimilPlotsDialog *dlg = static_cast<NtCConnSimilPlotsDialog *>(data);

    if (dlg->on_widget_resized) {
        dlg->on_widget_resized(event->configure.width, event->configure.height);
    }

    return FALSE;
}

static
void on_connectivity_point_selected(NtCConnSimilPlotsDialog *dlg, const DotPlotPoint &pt) {
}

static
void on_similarity_point_selected(NtCConnSimilPlotsDialog *dlg, const DotPlotPoint &pt) {
    if (dlg->on_similarity_point_selected) {
        dlg->on_similarity_point_selected(NtCSimilarity{{ pt.x, pt.y }, pt.caption});
    }
}

static
void prepare_altconfs_list(GtkComboBox *combo, const std::vector<std::string> &altconfs) {
    assert(combo);

    GtkListStore *store = gtk_list_store_new(3, G_TYPE_INT, G_TYPE_STRING, G_TYPE_STRING);
    gtk_combo_box_set_model(combo, GTK_TREE_MODEL(store));

    GtkCellRenderer *renderer = gtk_cell_renderer_text_new();
    gtk_cell_layout_pack_start(GTK_CELL_LAYOUT(combo), renderer, TRUE);
    gtk_cell_layout_set_attributes(GTK_CELL_LAYOUT(combo), renderer, "text", 1, NULL);

    fill_altconfs_list(store, altconfs, GtkSignalConnection{});

    gtk_combo_box_set_active(combo, 0);
}

static
PlotWidget * make_plot_widget(GtkDrawingArea *area, const std::string &xAxisTitle, const std::string &yAxisTitle, OnPointSelected onPointSelected) {
    return ntc_plot_widget_make(area, xAxisTitle, yAxisTitle, 1, onPointSelected);
}

static
void update_connectivity_caption(GtkLabel *label, std::string text, LLKA_NtC ntc) {
    static const std::string REPL("{}");
    assert(label);

    std::string ntcName = LLKA_NtCToName(ntc);
    size_t pos;
    while ((pos = text.find_first_of(REPL)) != std::string::npos) {
        text.replace(pos, REPL.length(), ntcName);
    }

    gtk_label_set_text(label, text.c_str());
}

static
void update_connectivity_next_caption(NtCConnSimilPlotsDialog *dlg, LLKA_NtC ntc) {
    update_connectivity_caption(dlg->connectivity_previous_caption, "Connectivity with previous step and {} used as reference NtC", ntc);
}

static
void update_connectivity_previous_caption(NtCConnSimilPlotsDialog *dlg, LLKA_NtC ntc) {
    update_connectivity_caption(dlg->connectivity_next_caption, "Connectivity with next step and {} used as reference NtC", ntc);
}

static
std::pair<double, double> rmsd_to_semaphore(double rmsd) {
    const double MinRmsd = 0.0;
    const double MaxRmsd = 1.0;
    const double Half = MaxRmsd / 2.0;

    double normalized = (rmsd + MinRmsd) / MaxRmsd;
    if (normalized > MaxRmsd)
        normalized = MaxRmsd;
    double red = (1.0 * (2 * normalized < 1 ? 2 * normalized : 1));
    double green = (1.0 * (1 - 2 * (normalized - Half > 0 ? normalized - Half : 0)));

    return { red, green };
}

void ntc_csp_dialog_destroy(NtCConnSimilPlotsDialog *dlg) {
    if (!dlg) {
        return;
    }

    if (!dlg->destroyed) {
        gtk_widget_destroy(GTK_WIDGET(dlg->root));
    }

    delete dlg;
}

bool ntc_csp_dialog_is_valid(NtCConnSimilPlotsDialog *dlg) {
    return !dlg->destroyed;
}

NtCConnSimilPlotsDialog * ntc_csp_dialog_make(OnSimilaritySelected onSimilaritySelected, OnWidgetResized onWidgetResized) {
    GtkBuilder *b = gtk_builder_new_from_file(get_glade_file("conn_simil_plots_dialog.glade").c_str());
    assert(b);

    NtCConnSimilPlotsDialog *dlg = new NtCConnSimilPlotsDialog{};

    dlg->on_similarity_point_selected = onSimilaritySelected;
    dlg->on_widget_resized = onWidgetResized;

    dlg->root = GTK_DIALOG(get_widget(b, "dialog"));
    assert(dlg->root);

    dlg->connectivity_previous_caption = GTK_LABEL(get_widget(b, "connectivity_previous_caption"));
    assert(dlg->connectivity_previous_caption);
    dlg->connectivity_next_caption = GTK_LABEL(get_widget(b, "connectivity_next_caption"));
    assert(dlg->connectivity_next_caption);

    GtkDrawingArea *similarity_area = GTK_DRAWING_AREA(get_widget(b, "similarity_plot"));
    assert(similarity_area);
    dlg->similarity_plot = make_plot_widget(similarity_area, "Cartesian RMSD [\xE2\x84\xAB]", "Euclidean distance", [dlg](const DotPlotPoint &pt) { on_similarity_point_selected(dlg, pt); } );
    assert(dlg->similarity_plot);

    GtkDrawingArea *connectivity_previous_area = GTK_DRAWING_AREA(get_widget(b, "connectivity_previous_plot"));
    assert(connectivity_previous_area);
    dlg->connectivity_previous_plot = make_plot_widget(connectivity_previous_area, "C5 [\xE2\x84\xAB]", "O3 [\xE2\x84\xAB]", [dlg](const DotPlotPoint &p) {});
    assert(dlg->connectivity_previous_plot);

    GtkDrawingArea *connectivity_next_area = GTK_DRAWING_AREA(get_widget(b, "connectivity_next_plot"));
    assert(connectivity_next_area);
    dlg->connectivity_next_plot = make_plot_widget(connectivity_next_area, "C5 [\xE2\x84\xAB]", "O3 [\xE2\x84\xAB]", [dlg](const DotPlotPoint &p) {});
    assert(dlg->connectivity_next_plot);

    dlg->previous_altconfs = GTK_COMBO_BOX(get_widget(b, "previous_step_altconfs"));
    assert(dlg->previous_altconfs);
    prepare_altconfs_list(dlg->previous_altconfs, { "" });
    dlg->previous_altconf_changed_sgc = GtkSignalConnection{dlg->previous_altconfs, "changed", G_CALLBACK(on_altconf_changed), dlg};

    dlg->next_altconfs = GTK_COMBO_BOX(get_widget(b, "next_step_altconfs"));
    assert(dlg->next_altconfs);
    prepare_altconfs_list(dlg->next_altconfs, { "" });
    dlg->next_altconf_changed_sgc = GtkSignalConnection{dlg->next_altconfs, "changed", G_CALLBACK(on_altconf_changed), dlg};

    GtkButton *close = GTK_BUTTON(get_widget(b, "close_button"));
    assert(close);
    g_signal_connect(close, "clicked", G_CALLBACK(on_close_button_clicked), dlg);

    GtkWidget *prev_altconf_cap = GTK_WIDGET(get_widget(b, "connectivity_previous_altconf_caption"));
    assert(prev_altconf_cap);
    gtk_widget_set_tooltip_text(prev_altconf_cap, PrevAltconfTooltipText);

    GtkWidget *next_altconf_cap = GTK_WIDGET(get_widget(b, "connectivity_next_altconf_caption"));
    assert(next_altconf_cap);
    gtk_widget_set_tooltip_text(next_altconf_cap, NextAltconfTooltipText);

    g_signal_connect(dlg->root, "destroy", G_CALLBACK(on_dialog_closed), dlg);

    dlg->destroyed = false;

    g_object_unref(b);

    return dlg;
}

void ntc_csp_dialog_show(NtCConnSimilPlotsDialog *dlg, int width, int height) {
    if (width > 0 && height > 0) {
        gtk_window_set_default_size(GTK_WINDOW(dlg->root), width, height);
    }
    gtk_widget_show(GTK_WIDGET(dlg->root));
}

void ntc_csp_dialog_update_connectivities(NtCConnSimilPlotsDialog *dlg, const std::shared_ptr<NtCConnectivities> &connectivities, LLKA_NtC ntc) {
    dlg->connectivities = connectivities;

    update_connectivity_previous_caption(dlg, ntc);
    update_connectivity_next_caption(dlg, ntc);

    fill_altconfs_list(GTK_LIST_STORE(gtk_combo_box_get_model(dlg->previous_altconfs)), gather_altconfs(dlg->connectivities->previous), dlg->previous_altconf_changed_sgc);
    gtk_combo_box_set_active(dlg->previous_altconfs, 0);

    fill_altconfs_list(GTK_LIST_STORE(gtk_combo_box_get_model(dlg->next_altconfs)), gather_altconfs(dlg->connectivities->next), dlg->next_altconf_changed_sgc);
    gtk_combo_box_set_active(dlg->next_altconfs, 0);

    display_connectivities(dlg);
}

void ntc_csp_dialog_update_similarities(NtCConnSimilPlotsDialog *dlg, const std::vector<NtCSimilarity> &similarities) {
    std::vector<DotPlotPoint> points;

    for (const NtCSimilarity &s : similarities) {
        std::pair<double, double> redGreen = rmsd_to_semaphore(s.similarity.rmsd);
        points.emplace_back(s.similarity.rmsd, s.similarity.euclideanDistance, redGreen.first, redGreen.second, 0.0, s.NtC);
    }

    ntc_plot_widget_set_points(dlg->similarity_plot, std::move(points));
    ntc_plot_widget_reset_zoom(dlg->similarity_plot);
}

void * ntc_csp_dialog_widget(NtCConnSimilPlotsDialog *dlg) {
    if (!dlg || dlg->destroyed) {
        return nullptr;
    }

    return dlg->root;
}
