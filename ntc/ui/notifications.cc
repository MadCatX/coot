// vim: set sw=4 ts=4 sts=4 expandtab :

#include "notifications.hh"

#include <gtk/gtk.h>

static const std::string NTC_MODULE_INFO = "NtC module information";

static
void run_modal_dialog(const std::string &title, const std::string &text, GtkMessageType type, GtkWindow *parent = nullptr) {
    GtkWidget *dlg = gtk_message_dialog_new(
        parent,
        GTK_DIALOG_MODAL,
        type,
        GTK_BUTTONS_OK,
        "%s",
        text.c_str()
    );
    gtk_window_set_title(GTK_WINDOW(dlg), title.c_str());

    gtk_dialog_run(GTK_DIALOG(dlg));

    gtk_widget_destroy(dlg);
}

void ntc_notify_error(const std::string &text) {
    run_modal_dialog(NTC_MODULE_INFO, text, GTK_MESSAGE_ERROR);
}

void ntc_notify_info(const std::string &text) {
    run_modal_dialog(NTC_MODULE_INFO, text, GTK_MESSAGE_INFO);
}

void ntc_notify_warning(const std::string &text) {
    run_modal_dialog(NTC_MODULE_INFO, text, GTK_MESSAGE_WARNING);
}
