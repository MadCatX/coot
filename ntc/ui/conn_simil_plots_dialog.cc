// vim: set sw=4 ts=4 sts=4 expandtab :

#include "conn_simil_plots_dialog.hh"

#include "common.hh"
#include "plot_widget.hh"

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

    OnSimilaritySelected on_similarity_point_selected;
    OnWidgetResized on_widget_resized;

    bool destroyed;
};

static
double calculate_axis_maximum(const std::vector<DotPlotPoint> &points, double DotPlotPoint::* coord, double scaleFactor);

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

    g_signal_connect(dlg->root, "configure-event", G_CALLBACK(on_configure_event), dlg);

    GtkButton *close = GTK_BUTTON(get_widget(b, "close_button"));
    assert(close);
    g_signal_connect(close, "clicked", G_CALLBACK(on_close_button_clicked), dlg);

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

void ntc_csp_dialog_update_connectivities(NtCConnSimilPlotsDialog *dlg, const NtCConnectivities &connectivities, LLKA_NtC ntc) {
    std::vector<DotPlotPoint> prevPoints;
    std::vector<DotPlotPoint> nextPoints;

    for (const auto &prev : connectivities.previous) {
        prevPoints.emplace_back(prev.connectivity.C5PrimeDistance, prev.connectivity.O3PrimeDistance, 1.0, 1.0, 0.0, prev.NtC);
    }

    for (const auto &next: connectivities.next) {
        nextPoints.emplace_back(next.connectivity.C5PrimeDistance, next.connectivity.O3PrimeDistance, 0.0, 1.0, 1.0, next.NtC);
    }

    update_connectivity_previous_caption(dlg, ntc);
    update_connectivity_next_caption(dlg, ntc);

    ntc_plot_widget_set_points(dlg->connectivity_previous_plot, std::move(prevPoints));
    ntc_plot_widget_set_points(dlg->connectivity_next_plot, std::move(nextPoints));

    ntc_plot_widget_reset_zoom(dlg->connectivity_previous_plot);
    ntc_plot_widget_reset_zoom(dlg->connectivity_next_plot);
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
