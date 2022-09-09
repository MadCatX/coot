// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_PLOT_WIDGET_HH
#define _NTC_UI_PLOT_WIDGET_HH

#include "dot_plot.hh"

#include <gtk/gtk.h>

#include <functional>

struct PlotWidget;
using OnPointSelected = std::function<void (DotPlotPoint)>;

struct PlotWidget * ntc_plot_widget_make(GtkDrawingArea *area, const std::string &xAxisTitle, const std::string &yAxisTitle, int tickStep, OnPointSelected onPointSelected = nullptr);
void ntc_plot_widget_reset_zoom(PlotWidget *pw);
void ntc_plot_widget_set_points(PlotWidget *pw, std::vector<DotPlotPoint> points);

#endif // _NTC_UI_PLOT_WIDGET_HH

