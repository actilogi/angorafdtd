
AM_CPPFLAGS = \
	-static \
	-I$(blitzdir)/include \
	-I$(configdir)/include \
	-I$(top_srcdir)/src

AM_LDFLAGS = -static

noinst_LTLIBRARIES = libnetwork.la
libnetwork_la_SOURCES = Cnetwork.cpp Cnetwork.h \
	Cport.cpp Cport.h \
	read_network.cpp read_network.h

