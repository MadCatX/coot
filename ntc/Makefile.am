lib_LTLIBRARIES = libcoot-ntc.la

AM_CPPFLAGS = \
	@GTK_CFLAGS@ \
	-I$(top_srcdir)

AM_CXXFLAGS = \
	@GTK_CFLAGS@ \
	-I$(top_srcdir)


libcoot_ntc_la_LDFLAGS = $(SHARED_LDFLAGS)

if COND_LLKA_CUSTOM_PREFIX
libcoot_ntc_la_LIBADD = -L$(LLKA_PREFIX)/lib/ -lLLKA
endif

pkginclude_HEADERS = ntc.hh
pkgincludedir = $(includedir)/coot/ntc

libcoot_ntc_la_SOURCES = ntc.cc util.cc \
	ui/ntc_dialog.cc ui/util.cc
