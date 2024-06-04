##
## TODO: Link only to netcdf if needed!
##

if ENABLE_DEBUG
AM_CFLAGS = -DDEBUG -O0 -g  -ansi -pedantic
DEFS += -DDEBUG -O0 -g  -ansi -pedantic
endif

if ENABLE_STATIC
AM_LDFLAGS = $(STATIC_LIBS)
endif

bin_PROGRAMS = grd2stream cdist2grd
noinst_HEADERS = debug_printf.h grdio.h log.h gmtopts.h

## make sure we keep them
EXTRA_DIST = own_getopt.c own_getopt.h

##
## grd2stream
##
grd2stream_CFLAGS = 
grd2stream_LDFLAGS = $(STATIC_LIBS)
grd2stream_LDADD =
grd2stream_SOURCES = grdio.c debug_printf.c log.c  grd2stream.c
if !SYSTEM_GETOPT
grd2stream_SOURCES += own_getopt.c
noinst_HEADERS += own_getopt.h
endif

##
## cdist2grd
##
cdist2grd_CFLAGS =
cdist2grd_LDFLAGS = $(STATIC_LIBS)
cdist2grd_LDADD =
cdist2grd_SOURCES  = cdist2grd.c gmtopts.c grdio.c debug_printf.c
if !SYSTEM_GETOPT
cdist2grd_SOURCES += own_getopt.c
noinst_HEADERS += own_getopt.h
endif


CLEANFILES = *.c~ *.h~

