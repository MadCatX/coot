// vim: set sw=4 ts=4 sts=4 expandtab :

#include "plot_widget.hh"

#include "dot_plot.hh"

#include <cassert>
#include <limits>
#include <vector>

struct PlotWidget {
    std::vector<DotPlotPoint> data_points;
    double x_minimum;
    double x_maximum;
    double y_minimum;
    double y_maximum;

    DotPlotAxis x_axis;
    DotPlotAxis y_axis;

    DotPlotOptions options;

    // Pointer tracking
    int pointer_x;
    int pointer_y;

    OnPointSelected on_point_selected;

    GtkEventController *scroller;
    GtkDrawingArea *area;
};

static
void plot_zoom(PlotWidget *pw, gdouble delta);

static
double zoom(const DotPlotAxis &axis, gdouble delta) {
    assert(delta != 0.0);

    constexpr double SCALE = 0.25;

    double factor = SCALE * delta;
    double p = axis.maximum + axis.maximum * factor;

    return p > axis.minimum ? p : axis.maximum;
}

static
void zoom_reset(PlotWidget *pw) {
    pw->x_axis.minimum = pw->x_minimum;
    pw->x_axis.maximum = pw->x_maximum;
    pw->y_axis.minimum = pw->y_minimum;
    pw->y_axis.maximum = pw->y_maximum;
}

// Event handlers
static
gboolean on_plot_click_event(GtkDrawingArea *area, GdkEventButton *event, gpointer data) {
    PlotWidget *pw = static_cast<PlotWidget *>(data);

    if (event->type == GDK_BUTTON_PRESS && event->button == 1) {
        std::pair<bool, DotPlotPoint> ret = dot_plot_point_at_pixel(event->x, event->y, area, pw->x_axis, pw->y_axis, pw->options, pw->data_points);
        if (ret.first == true && pw->on_point_selected) {
            const DotPlotPoint &pt = ret.second;
            pw->on_point_selected(pt);
        }
    } else if (event->type == GDK_2BUTTON_PRESS) {
        zoom_reset(pw);
        gtk_widget_queue_draw(GTK_WIDGET(pw->area));
    }

    return TRUE;
}

static
gboolean on_plot_draw(GtkDrawingArea *area, cairo_t *cr, gpointer data) {
    PlotWidget *pw = static_cast<PlotWidget *>(data);

    dot_plot_draw(area, cr, pw->x_axis, pw->y_axis, pw->options, pw->data_points);

    return TRUE;
}

static
gboolean on_plot_scroll_event(GtkEventControllerScroll *scroller, gdouble dx, gdouble dy, gpointer data) {
    PlotWidget *pw = static_cast<PlotWidget *>(data);

    plot_zoom(pw, dy);

    return TRUE;
}

static
gboolean on_plot_pointer_motion_event(GtkDrawingArea *area, GdkEventMotion *motion, gpointer data) {
    PlotWidget *pw = static_cast<PlotWidget *>(data);

    gboolean handled = FALSE;

    if (motion->state & GDK_BUTTON1_MASK) {
        int dx = motion->x - pw->pointer_x;
        int dy = motion->y - pw->pointer_y;

        DotPlotAxisShift shift = dot_plot_axis_shift(dx, dy, pw->area, pw->x_axis, pw->y_axis, pw->options);

        pw->x_axis.minimum -= shift.sx;
        pw->x_axis.maximum -= shift.sx;

        pw->y_axis.minimum += shift.sy;
        pw->y_axis.maximum += shift.sy;

        gtk_widget_queue_draw(GTK_WIDGET(pw->area));

        handled = TRUE;
    }

    pw->pointer_x = motion->x;
    pw->pointer_y = motion->y;

    return handled;
}

static
void plot_zoom(PlotWidget *pw, gdouble delta) {
    if (delta == 0.0)
        return;

    pw->x_axis.maximum = zoom(pw->x_axis, delta);
    pw->y_axis.maximum = zoom(pw->y_axis, delta);

    gtk_widget_queue_draw(GTK_WIDGET(pw->area));
}

PlotWidget * ntc_plot_widget_make(GtkDrawingArea *area, const std::string &xAxisTitle, const std::string &yAxisTitle, int tickStep, OnPointSelected onPointSelected) {
    PlotWidget *pw = new PlotWidget;

    pw->area = area;

    pw->x_minimum = 0.0;
    pw->x_maximum = 1.0;
    pw->y_minimum = 0.0;
    pw->y_maximum = 1.0;

    pw->x_axis.minimum = pw->x_minimum;
    pw->x_axis.maximum = pw->x_maximum;
    pw->x_axis.tickStep = tickStep;
    pw->x_axis.title = xAxisTitle.c_str();

    pw->y_axis.minimum = pw->y_minimum;
    pw->y_axis.maximum = pw->y_maximum;
    pw->y_axis.tickStep = tickStep;
    pw->y_axis.title = yAxisTitle.c_str();

    pw->options.leftOffsetPx = 50;
    pw->options.rightOffsetPx = 10;
    pw->options.topOffsetPx = 10;
    pw->options.bottomOffsetPx = 50;
    pw->options.dotRadius = 6;
    pw->options.tickLength = 5;

    pw->on_point_selected = onPointSelected;

    pw->scroller = gtk_event_controller_scroll_new(GTK_WIDGET(area), GTK_EVENT_CONTROLLER_SCROLL_VERTICAL);

    // Wire up event handlers
    gtk_widget_add_events(GTK_WIDGET(area), GDK_BUTTON_PRESS_MASK | GDK_SCROLL_MASK | GDK_POINTER_MOTION_MASK);
    g_signal_connect(pw->scroller, "scroll", G_CALLBACK(on_plot_scroll_event), pw);
    g_signal_connect(area, "button-press-event", G_CALLBACK(on_plot_click_event), pw);
    g_signal_connect(area, "motion-notify-event", G_CALLBACK(on_plot_pointer_motion_event), pw);

    g_signal_connect(area, "draw", G_CALLBACK(on_plot_draw), pw);

    return pw;
}

void ntc_plot_widget_reset_zoom(PlotWidget *pw) {
    zoom_reset(pw);

    gtk_widget_queue_draw(GTK_WIDGET(pw->area));
}

void ntc_plot_widget_set_points(PlotWidget *pw, std::vector<DotPlotPoint> points) {
    pw->data_points = std::move(points);

    if (pw->data_points.empty()) {
        pw->x_minimum = 0.0;
        pw->x_maximum = 1.0;
        pw->y_minimum = 0.0;
        pw->y_maximum = 1.0;

        zoom_reset(pw);
    } else {
        pw->x_minimum = std::numeric_limits<double>::max();
        pw->x_maximum = std::numeric_limits<double>::min();
        pw->y_minimum = std::numeric_limits<double>::max();
        pw->y_maximum = std::numeric_limits<double>::min();

        for (const DotPlotPoint &pt : pw->data_points) {
            pw->x_minimum = (pw->x_minimum > pt.x) * pt.x + !(pw->x_minimum > pt.x) * pw->x_minimum;
            pw->x_maximum = (pw->x_maximum < pt.x) * pt.x + !(pw->x_maximum < pt.x) * pw->x_maximum;

            pw->y_minimum = (pw->y_minimum > pt.y) * pt.y + !(pw->y_minimum > pt.y) * pw->y_minimum;
            pw->y_maximum = (pw->y_maximum < pt.y) * pt.y + !(pw->y_maximum < pt.y) * pw->y_maximum;
        }
    }

    gtk_widget_queue_draw(GTK_WIDGET(pw->area));
}
