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

bin_PROGRAMS = grd2stream
noinst_HEADERS = debug_printf.h grdio.h log.h #options.h

##
## grd2stream
##
grd2stream_CFLAGS = 
grd2stream_LDFLAGS = $(STATIC_LIBS)
grd2stream_LDADD =   #-lnetcdf
grd2stream_SOURCES  = grdio.c debug_printf.c log.c options.c grd2stream.c

if GGO
BUILT_SOURCES = options.c
options.c options.h: options.ggo
#	$(GENGETOPT) -u -N -i options.ggo -f options -F options
	$(GENGETOPT) --unamed-opts --no-version --input=options.ggo --func-name=options --file-name=options
endif



CLEANFILES = *.c~ *.h~

