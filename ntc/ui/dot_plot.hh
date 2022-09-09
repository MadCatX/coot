// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_DOT_PLOT_HH
#define _NTC_UI_DOT_PLOT_HH

#include <gtk/gtk.h>

#include <string>
#include <vector>

struct DotPlotAxis {
    double minimum;
    double maximum;

    int tickStep;

    std::string title;
};

struct DotPlotAxisShift {
    double sx;
    double sy;
};

struct DotPlotOptions {
    int leftOffsetPx;
    int rightOffsetPx;

    int topOffsetPx;
    int bottomOffsetPx;

    int tickLength;
    int dotRadius;
};

struct DotPlotPoint {
    DotPlotPoint() : x{0}, y{0}, color{0, 0, 0} {}
    DotPlotPoint(
        double x, double y,
        double red, double green, double blue,
        std::string caption
    ) : x{x}, y{y}, color{red, green, blue}, caption{std::move(caption)}
    {}

    double x;
    double y;
    struct {
        double red;
        double green;
        double blue;
    } color;
    std::string caption;
};

DotPlotAxisShift dot_plot_axis_shift(int dxPx, int dxPy, GtkDrawingArea *area, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options);
void dot_plot_draw(GtkDrawingArea *area, cairo_t *cr, DotPlotAxis xAxis, DotPlotAxis yAxis, const DotPlotOptions &options, const std::vector<DotPlotPoint> &points);
std::pair<bool, DotPlotPoint> dot_plot_point_at_pixel(int x, int y, GtkDrawingArea *area, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options, const std::vector<DotPlotPoint> &points);

#endif // _NTC_UI_DOT_PLOT_HH
