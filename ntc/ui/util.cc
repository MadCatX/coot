// vim: set sw=4 ts=4 sts=4 expandtab :

#include "common.hh"

#include <gtk/gtk.h>

std::string pick_ntc_parameters_directory() {
    std::string path{};

    GtkWidget *fileChooser = gtk_file_chooser_dialog_new(
        "Specify path to the directory with NtC parameters file",
        nullptr,
        GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER,
        "Cancel",
        GTK_RESPONSE_CANCEL,
        "Open",
        GTK_RESPONSE_ACCEPT,
        NULL
    );

    gint ret = gtk_dialog_run(GTK_DIALOG(fileChooser));
    if (ret == GTK_RESPONSE_ACCEPT) {
        char *cpath = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(fileChooser));
        path = std::string{cpath};
        g_free(cpath);
    }

    gtk_widget_destroy(fileChooser);

    return path;
}
