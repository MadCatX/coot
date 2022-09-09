// vim: set sw=4 ts=4 sts=4 expandtab :

#include "dot_plot.hh"

#include <cassert>
#include <cmath>
#include <string>

struct Coordinate {
    double x;
    double y;
};

struct Pixel {
    int x;
    int y;
};

struct Color {
    double red;
    double green;
    double blue;
};

static Color TEXT_COLOR = { 0.0, 0.0, 0.0 };

static
std::string format_number(double value);

static
int axis_length(int dimension, int leadOffset, int tailOffset) {
    return dimension - leadOffset - tailOffset;
}

static
double tick_step(const DotPlotAxis &axis) {
    assert(axis.maximum > axis.minimum);

    double fullRange = axis.maximum - axis.minimum;

    // Cut off all digits except the most significant digit
    double scale = std::pow(10.0, std::floor(std::log10(fullRange)));
    double step = int(fullRange / scale) * scale;

    return step;
}

static
double first_tick_value(int oneTickPx, double oneTickRange, const DotPlotAxis &axis) {
    assert(oneTickRange > 0.0);

    double a = axis.minimum / oneTickRange;
    double value = (std::floor(a) + 1.0) * oneTickRange;

    return value;
}

static
int first_tick_offset_in_pixels(int oneTickPx, double oneTickRange, const DotPlotAxis &axis) {
    double ftv = first_tick_value(oneTickPx, oneTickRange, axis);
    double diff = ftv - axis.minimum;

    assert(diff >= 0.0);

    int result = int(oneTickPx * diff / oneTickRange);

    return result;
}

static
double one_tick_in_pixels(int dimension, double oneTickRange, const DotPlotAxis &axis) {
    assert(axis.maximum - axis.minimum);

    int rng = int(dimension * oneTickRange / (axis.maximum - axis.minimum));

    return rng > 1 ? rng : 1;
}

static
double one_tick_range(const DotPlotAxis &axis) {
    double stepping = std::pow(10.0, axis.tickStep);
    double step = tick_step(axis);
    double result = step / stepping;

    return result;
}

static
double pixel_scaling(int canvasDimension, int leadAxisOffset, int tailAxisOffset, const DotPlotAxis &axis) {
    double axisRange = axis.maximum - axis.minimum;
    int axisLength = axis_length(canvasDimension, leadAxisOffset, tailAxisOffset);
    double oneTickRange = one_tick_range(axis);
    int oneTickPx = one_tick_in_pixels(axisLength, oneTickRange, axis);
    double nTicks = axisRange / oneTickRange;
    int pxMax = oneTickPx * nTicks;

    return pxMax / axisRange;
}

static
Pixel to_pixel_coords(const Coordinate &pt, int canvasHeightPx, double pixelScaleX, double pixelScaleY, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options);

static
void draw_background(cairo_t *cr, guint width, guint height, const Color &color) {
    cairo_save(cr);

    cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.0);
    cairo_rectangle(cr, 0, 0, width, height);
    cairo_fill(cr);

    cairo_restore(cr);
}

static
void draw_points(cairo_t *cr, const std::vector<DotPlotPoint> &points, int canvasWidthPx, int canvasHeightPx, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options) {
    double pixelScaleX = pixel_scaling(canvasWidthPx, options.leftOffsetPx, options.rightOffsetPx, xAxis);
    double pixelScaleY = pixel_scaling(canvasHeightPx, options.topOffsetPx, options.bottomOffsetPx, yAxis);

    cairo_save(cr);

    // Set up clipping so that we do not draw outside the axes
    cairo_new_sub_path(cr);
    cairo_move_to(cr, options.leftOffsetPx, options.topOffsetPx);
    cairo_line_to(cr, canvasWidthPx - options.rightOffsetPx, options.topOffsetPx);
    cairo_line_to(cr, canvasWidthPx - options.rightOffsetPx, canvasHeightPx - options.bottomOffsetPx);
    cairo_line_to(cr, options.leftOffsetPx, canvasHeightPx - options.bottomOffsetPx);
    cairo_close_path(cr);
    cairo_clip(cr);

    // Draw points first
    for (const DotPlotPoint &pt : points) {
        Pixel pixel = to_pixel_coords({ pt.x, pt.y }, canvasHeightPx, pixelScaleX, pixelScaleY, xAxis, yAxis, options);
        cairo_set_source_rgba(cr, pt.color.red, pt.color.green, pt.color.blue, 1.0);
        cairo_arc(cr, pixel.x, pixel.y, options.dotRadius, 0, 2.0 * M_PI);
        cairo_fill(cr);
    }

    // Draw text labels on top
    cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_source_rgba(cr, TEXT_COLOR.red, TEXT_COLOR.green, TEXT_COLOR.blue, 1.0);
    for (const DotPlotPoint &pt : points) {
        Pixel pixel = to_pixel_coords({ pt.x, pt.y }, canvasHeightPx, pixelScaleX, pixelScaleY, xAxis, yAxis, options);

        cairo_text_extents_t extents;
        cairo_text_extents(cr, pt.caption.c_str(), &extents);
        cairo_move_to(cr, pixel.x - extents.width / 2, pixel.y - extents.height - 2);
        cairo_show_text(cr, pt.caption.c_str());
    }

    cairo_restore(cr);
}

static
void draw_x_axis(cairo_t *cr, int canvasWidthPx, int canvasHeightPx, const DotPlotAxis &axis, const DotPlotOptions &options, const Color &color) {
    int y = canvasHeightPx - options.bottomOffsetPx;

    int axisLength = axis_length(canvasWidthPx, options.leftOffsetPx, options.rightOffsetPx);
    int xStop = axisLength + options.leftOffsetPx;

    double oneTickRange = one_tick_range(axis);
    int oneTickPx = one_tick_in_pixels(axisLength, oneTickRange, axis);
    int firstTickOffsetPx = first_tick_offset_in_pixels(oneTickPx, oneTickRange, axis);
    double firstTickValue = first_tick_value(oneTickPx, oneTickRange, axis);

    cairo_save(cr);

    cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.0);
    cairo_set_line_width(cr, 1.0);

    // Main line
    cairo_move_to(cr, options.leftOffsetPx, y);
    cairo_line_to(cr, options.leftOffsetPx + axisLength, y);

    // Ticks
    int x = options.leftOffsetPx + firstTickOffsetPx;
    while (x < xStop) {
        if (x >= options.leftOffsetPx) {
            cairo_move_to(cr, x, y);
            cairo_line_to(cr, x, y + options.tickLength);
        }

        x += oneTickPx;
    }

    // Text captions
    cairo_text_extents_t extents;
    cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    int n = 0;
    x = options.leftOffsetPx + firstTickOffsetPx;
    while (x < xStop) {
        if (x >= options.leftOffsetPx) {
            double value = n * oneTickRange + firstTickValue;

            std::string text = format_number(value);
            cairo_text_extents(cr, text.c_str(), &extents);
            int textX = x - extents.width / 2;

            cairo_move_to(cr, textX, y + extents.height + options.tickLength + 2);
            cairo_show_text(cr, text.c_str());
        }

        x += oneTickPx;
        n++;
    }

    // Axis title
    cairo_text_extents(cr, axis.title.c_str(), &extents);
    cairo_move_to(cr, (canvasWidthPx - extents.width) / 2, canvasHeightPx - extents.height - 1);
    cairo_show_text(cr, axis.title.c_str());

    cairo_stroke(cr);

    cairo_restore(cr);
}

static
void draw_y_axis(cairo_t *cr, int canvasWidthPx, int canvasHeightPx, const DotPlotAxis &axis, const DotPlotOptions &options, const Color &color) {
    int x = options.leftOffsetPx;
    int yOrigin = canvasHeightPx - options.bottomOffsetPx;

    int axisLength = axis_length(canvasHeightPx, options.topOffsetPx, options.bottomOffsetPx);
    int yStop = options.topOffsetPx;

    double oneTickRange = one_tick_range(axis);
    int oneTickPx = one_tick_in_pixels(axisLength, oneTickRange, axis);
    int firstTickOffsetPx = first_tick_offset_in_pixels(oneTickPx, oneTickRange, axis);
    double firstTickValue = first_tick_value(oneTickPx, oneTickRange, axis);

    cairo_save(cr);

    cairo_set_source_rgba(cr, color.red, color.green, color.blue, 1.0);
    cairo_set_line_width(cr, 1.0);

    // Main line
    cairo_move_to(cr, x, yOrigin);
    cairo_line_to(cr, x, yOrigin - axisLength);

    // Ticks
    int y = yOrigin - firstTickOffsetPx;
    while (y > yStop) {
        cairo_move_to(cr, x, y);
        cairo_line_to(cr, x - options.tickLength, y);

        y -= oneTickPx;
    }

    // Text captions
    cairo_text_extents_t extents;
    cairo_select_font_face(cr, "monospace", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    int n = 0;
    y = yOrigin - firstTickOffsetPx;
    while (y > yStop) {
        double value = n * oneTickRange + firstTickValue;

        std::string text = format_number(value);
        cairo_text_extents(cr, text.c_str(), &extents);
        int textY = y + extents.height / 2;

        cairo_move_to(cr, x - options.tickLength - 2 - extents.width, textY);
        cairo_show_text(cr, text.c_str());

        y -= oneTickPx;
        n++;
    }

    // Axis title
    cairo_rotate(cr, -M_PI / 2);
    cairo_text_extents(cr, axis.title.c_str(), &extents);
    cairo_move_to(cr, -(canvasHeightPx + extents.width) / 2, extents.height + 1);
    cairo_show_text(cr, axis.title.c_str());

    cairo_stroke(cr);

    cairo_restore(cr);
}

static
std::string format_number(double value) {
    size_t len = snprintf(nullptr, 0, "%g", value);
    char *buf = new char[len + 1];

    snprintf(buf, len + 1, "%g", value);

    std::string str{buf};
    delete [] buf;
    return str;
}

static
Pixel to_pixel_coords(const Coordinate &pt, int canvasHeightPx, double pixelScaleX, double pixelScaleY, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options) {
    int xPx = pixelScaleX * (pt.x - xAxis.minimum) + options.leftOffsetPx;
    int yPx = canvasHeightPx - pixelScaleY * (pt.y - yAxis.minimum) - options.bottomOffsetPx;

    return { xPx, yPx };
}

DotPlotAxisShift dot_plot_axis_shift(int dxPx, int dyPx, GtkDrawingArea *area, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options) {
    assert(xAxis.maximum > xAxis.minimum);
    assert(yAxis.maximum > yAxis.minimum);

    guint canvasWidthPx = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    guint canvasHeightPx = gtk_widget_get_allocated_height(GTK_WIDGET(area));

    double xAxisRange = xAxis.maximum - xAxis.minimum;
    double sx = xAxisRange * dxPx / (canvasWidthPx - options.leftOffsetPx - options.rightOffsetPx);

    double yAxisRange = yAxis.maximum - yAxis.minimum;
    double sy = yAxisRange * dyPx / (canvasHeightPx - options.topOffsetPx - options.bottomOffsetPx);

    return { sx, sy };
}

void dot_plot_draw(GtkDrawingArea *area, cairo_t *cr, DotPlotAxis xAxis, DotPlotAxis yAxis, const DotPlotOptions &options, const std::vector<DotPlotPoint> &points) {
    assert(area);

    guint width = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    guint height = gtk_widget_get_allocated_height(GTK_WIDGET(area));
    GtkStyleContext *ctx = gtk_widget_get_style_context(GTK_WIDGET(area));
    assert(ctx);

    Color bgColor = { 1.0, 1.0, 1.0 };
    Color axisColor = { 0.0, 0.0, 0.0 };

    draw_background(cr, width, height, bgColor);
    draw_x_axis(cr, width, height, xAxis, options, axisColor);
    draw_y_axis(cr, width, height, yAxis, options, axisColor);

    draw_points(cr, points, width, height, xAxis, yAxis, options);
}

std::pair<bool, DotPlotPoint> dot_plot_point_at_pixel(int x, int y, GtkDrawingArea *area, const DotPlotAxis &xAxis, const DotPlotAxis &yAxis, const DotPlotOptions &options, const std::vector<DotPlotPoint> &points) {
    guint canvasWidthPx = gtk_widget_get_allocated_width(GTK_WIDGET(area));
    guint canvasHeightPx = gtk_widget_get_allocated_height(GTK_WIDGET(area));

    double pixelScaleX = pixel_scaling(canvasWidthPx, options.leftOffsetPx, options.rightOffsetPx, xAxis);
    double pixelScaleY = pixel_scaling(canvasHeightPx, options.topOffsetPx, options.bottomOffsetPx, yAxis);

    for (const DotPlotPoint &pt : points) {
        auto pixel = to_pixel_coords({ pt.x, pt.y }, canvasHeightPx, pixelScaleX, pixelScaleY, xAxis, yAxis, options);

        int dxx = (x - pixel.x) * (x - pixel.x);
        int dyy = (y - pixel.y) * (y - pixel.y);

        if (dxx + dyy < options.dotRadius * options.dotRadius) {
            return { true, pt };
        }
    }

    return { false, {} };
}
