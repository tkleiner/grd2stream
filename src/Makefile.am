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
EXTRA_DIST = own_getopt.c own_getopt.h own_getopt.h

##
## grd2stream
##
grd2stream_CFLAGS = 
grd2stream_LDFLAGS = $(STATIC_LIBS)
grd2stream_LDADD =
grd2stream_SOURCES = grdio.c debug_printf.c log.c  grd2stream.c own_getopt.c

##
## cdist2grd
##
cdist2grd_CFLAGS =
cdist2grd_LDFLAGS = $(STATIC_LIBS)
cdist2grd_LDADD =
cdist2grd_SOURCES  = cdist2grd.c gmtopts.c grdio.c debug_printf.c own_getopt.c



CLEANFILES = *.c~ *.h~

