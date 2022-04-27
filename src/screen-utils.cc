
#include <string>
#include <iostream>
#include <sys/stat.h>
#include <gtk/gtk.h>
#include "utils/coot-utils.hh"

#include "compat/coot-sysdep.h"

// this returns the effective screen height if possible otherwise an estimate
int
get_max_effective_screen_height() {

    // using properties
    gboolean ok;
    guchar* raw_data = NULL;
    gint data_len = 0;
    gint width = 200;
    gint height = 200;
    int max_height;
    max_height = -1;

// no gdk_property get on windows (at the moment)
#if !defined WINDOWS_MINGW && !defined _MSC_VER
    ok = gdk_property_get(gdk_get_default_root_window(),  // a gdk window
                          gdk_atom_intern("_NET_WORKAREA", FALSE),  // property
                          gdk_atom_intern("CARDINAL", FALSE),  // property type
                          0,  // byte offset into property
                          0xff,  // property length to retrieve
                          false,  // delete property after retrieval?
                          NULL,  // returned property type
                          NULL,  // returned data format
                          &data_len,  // returned data len
                          &raw_data);  // returned data

    if (ok) {

        // We expect to get four longs back: x, y, width, height.
        if (data_len >= static_cast<gint>(4 * sizeof(glong))) {
            glong* data = reinterpret_cast<glong*>(raw_data);
            gint x = data[0];
            gint y = data[1];
            width = data[2];
            height = data[3];
            max_height = height;
        }
        g_free(raw_data);
    }
#endif // MINGW
    if (max_height < 0) {
        GdkScreen *screen;
        screen = gdk_screen_get_default();
        if (screen) {
           // width = gdk_screen_get_width(screen);
           // height = gdk_screen_get_height(screen);

#ifdef WINDOWS_MINGW
            max_height = int(height * 0.95);
#else
            max_height = int(height * 0.9);
#endif // MINGW
        } else {
            g_print ("BL ERROR:: couldnt get gdk screen; should never happen\n");
        }
    }
    return max_height;
}

int
setup_screen_size_settings() {

   int ret = 0;
   int max_height;
   max_height = get_max_effective_screen_height();

   // adjust the icons size of the refinement toolbar icons
   if (max_height <= 620) {
       max_height = 620;
       std::cout << "setup_screen_size_settings() Fix these screen settings " << std::endl;
       // gtk_rc_parse_string("gtk-icon-sizes=\"gtk-large-toolbar=10,10:gtk-button=10,10\"");
       // gtk_rc_parse_string("class \"GtkLabel\" style \"small-font\"");
       ret = 1;
   } else if (max_height <= 720) {
       int icon_size = 12 + (max_height - 620) / 25;
       // std::cout << "BL INFO:: screen has " << max_height << " height, will make "
       //           << "icons to " <<  icon_size <<std::endl;
       std::string toolbar_txt = "gtk-icon-sizes = \"gtk-large-toolbar=";
       toolbar_txt += coot::util::int_to_string(icon_size);
       toolbar_txt += ",";
       toolbar_txt += coot::util::int_to_string(icon_size);
       toolbar_txt += ":gtk-button=";
       toolbar_txt += coot::util::int_to_string(icon_size);
       toolbar_txt += ",";
       toolbar_txt += coot::util::int_to_string(icon_size);
       toolbar_txt += "\"";

       // deprecated - use GtkStyleContext
       // gtk_rc_parse_string (toolbar_txt.c_str());
   }
   return ret;
}


void setup_application_icon(GtkWindow *window) {

   std::string splash_screen_pixmap_dir = coot::package_data_dir();
   splash_screen_pixmap_dir += "/";
   splash_screen_pixmap_dir += "pixmaps";

   // over-ridden by user?
   char *s = getenv("COOT_PIXMAPS_DIR");
   if (s) {
      splash_screen_pixmap_dir = s;
   }

   // now add the application icon
   std::string app_icon_path = coot::util::append_dir_file(splash_screen_pixmap_dir, "coot-icon.png");

   struct stat buf;
   int status = stat(app_icon_path.c_str(), &buf);
   if (status == 0) { // icon file was found

      GdkPixbuf *icon_pixbuf =
	 gdk_pixbuf_new_from_file (app_icon_path.c_str(), NULL);
      if (icon_pixbuf) {
	 gtk_window_set_icon (GTK_WINDOW (window), icon_pixbuf);
	 g_object_unref(icon_pixbuf); // - what does this do?  I mean,
	 // why do I want to unref this pixbuf now? What does
	 // gtk_window_set_icon() do to the refcount? How do I know
	 // the refcount on a widget?
      }
   }

   // load svg/png files to antialias icons
   // maybe should go somewhere else.
   GtkIconSet* iconset;

   GtkIconFactory *iconfactory = gtk_icon_factory_new();
   GtkIconTheme *icon_theme = gtk_icon_theme_new();
   gtk_icon_theme_set_custom_theme(icon_theme, "coot");

   // GSList* stock_ids = gtk_stock_list_ids(); deprecated 
   GError *error = NULL;
   GdkPixbuf* pixbuf;

   std::vector<std::string> iconFiles = coot::gather_files_by_patterns(splash_screen_pixmap_dir, { "*.svg", "*.png" });
   for (const std::string &filepath : iconFiles) {
      pixbuf = gdk_pixbuf_new_from_file(filepath.c_str(), &error);
      if (error) {
	 g_print("Error loading icon: %s\n", error->message);
	 g_error_free(error);
	 error = NULL;
      } else if (pixbuf) {
         iconset = gtk_icon_set_new_from_pixbuf(pixbuf);
	 g_object_unref(pixbuf);

         std::string filename = coot::util::file_name_non_directory(filepath);
         if (!filename.empty()) {
            gtk_icon_factory_add(iconfactory, filename.c_str(), iconset);
            gtk_icon_factory_add_default(iconfactory);
         }
      } else {
         g_print("WARNING: No error was set but pixbuf was NULL.\n");
      }
   }
}
