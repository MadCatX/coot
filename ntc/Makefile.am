lib_LTLIBRARIES = libcoot-ntc.la

AM_CPPFLAGS = \
	@GTK_CFLAGS@ \
	-I$(top_srcdir) \
	-D_USE_MATH_DEFINES

AM_CXXFLAGS = \
	@GTK_CFLAGS@ \
	-I$(top_srcdir) \
	-D_USE_MATH_DEFINES

libcoot_ntc_la_LDFLAGS = $(SHARED_LDFLAGS)

libcoot_ntc_la_LIBADD = -lmmdb2 -lclipper-core \
                         $(top_builddir)/coords/libcoot-coords.la  \
                         $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
                         $(top_builddir)/utils/libcoot-utils.la \
			 $(GTK_LIBS)

if COND_LLKA_CUSTOM_PREFIX
libcoot_ntc_la_LIBADD += -L$(LLKA_PREFIX)/bin/ -lLLKA
endif

pkginclude_HEADERS = ntc.hh ui/ntc_dialog.hh ui/conn_simil_plots_dialog.hh ui/notifications.hh
pkgincludedir = $(includedir)/coot/ntc

libcoot_ntc_la_SOURCES = ntc.cc util.cc \
	ui/ntc_dialog.cc ui/util.cc ui/conn_simil_plots_dialog.cc \
	ui/dot_plot.cc ui/plot_widget.cc ui/notifications.cc
