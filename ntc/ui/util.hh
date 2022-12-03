// vim: set sw=4 ts=4 sts=4 expandtab :

#ifndef _NTC_UI_UTIL_HH
#define _NTC_UI_UTIL_HH

#include <gtk/gtk.h>

#include <cassert>
#include <string>

struct GtkSignalConnection {
public:
    GtkSignalConnection() :
        id{0UL},
        obj{nullptr}
    {}

    GtkSignalConnection(gulong id, GObject *obj) :
        id{id},
        obj{obj}
    {}

    template <typename TGObject>
    GtkSignalConnection(
        TGObject *obj,
        const char *detailed_signal,
        GCallback handler,
        gpointer data
    ) {
        id = g_signal_connect(obj, detailed_signal, handler, data);
        this->obj = G_OBJECT(obj);
        assert(this->obj);
    }

    operator bool() const { return obj != nullptr; }

    void block() {
        if (*this) {
            g_signal_handler_block(obj, id);
        }
    }

    void unblock() {
        if (*this) {
            g_signal_handler_unblock(obj, id);
        }
    }

private:
    gulong id;
    GObject *obj;
};

std::string pick_ntc_parameters_directory();

#endif // _NTC_UI_UTIL_HH
