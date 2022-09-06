// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_COMMON_HH
#define _NTC_UI_COMMON_HH

#include "utils/coot-utils.hh"

#include <gtk/gtk.h>

#include <cassert>

static
GtkWidget * get_widget(GtkBuilder *b, const char *name) {
    GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(b, name));
    assert(b);

    return w;
}

static
std::string get_glade_file(const std::string &fileName) {
    std::string dir = coot::package_data_dir();
    std::string dir_glade = coot::util::append_dir_dir(dir, "glade");
    return coot::util::append_dir_file(dir_glade, fileName);
}

#endif // _NTC_UI_COMMON_HH
