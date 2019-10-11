########################################################################
##
## Init
##
########################################################################
AC_INIT([grd2stream], [0.2.6], [tkleiner@awi.de])
AC_CONFIG_SRCDIR([configure.ac])
# Put helper scripts here to reduce top level clutter
AC_CONFIG_AUX_DIR([admin])
# Tell system about local m4 macros to keep libtool in sync
AC_CONFIG_MACRO_DIR([m4])
# Put preprocessor symbols in config.h
AC_CONFIG_HEADERS([config.h])
# Use automake to create Makefiles
AM_INIT_AUTOMAKE([tar-ustar])


dnl canonical host-system: host, host_cpu, host_vendor, and host_os 
AC_CANONICAL_HOST

dnl Add support for the Dmalloc package. 
dnl If the user runs configure with --with-dmalloc, then define WITH_DMALLOC and add -ldmalloc to LIBS. 
AM_WITH_DMALLOC

dnl
dnl Set the default install prefix to the installers home directory
dnl
AC_PREFIX_DEFAULT(${HOME})
AC_MSG_CHECKING([--prefix])
if test "x$prefix" = "xNONE" ; then
   prefix=$ac_default_prefix
fi
AC_MSG_RESULT([$prefix])


dnl 
dnl What should I build?
dnl 
AC_ARG_ENABLE(debug, AC_HELP_STRING([--enable-debug],
    [Compile in debug mode @<:@default=no@:>@]),
    enable_debug=$enableval,
    enable_debug=no)
AM_CONDITIONAL(ENABLE_DEBUG, test "x${enable_debug}" == "xyes" )

dnl
dnl Static libs
dnl 
AC_ARG_ENABLE(static,
        AC_HELP_STRING([--enable-static],
        [force static linking @<:@default=no@:>@]),
        [enable_static=$enableval],
        [enable_static=no])
if test "x${enable_static}" = "xyes"; then
        STATIC_LIBS=-static
fi
AC_SUBST(STATIC_LIBS)
AM_CONDITIONAL(ENABLE_STATIC, test "x${enable_static}" = "xyes")



# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB

# Checks for libraries.
#netcdf first!
AC_CHECK_LIB([m], [sqrt])
dnl AC_CHECK_LIB([netcdf], [nc_open])

# Checks for header files.
AC_HEADER_STDC
dnl AC_CHECK_HEADERS([stdlib.h string.h getopt.h netcdf.h])
AC_CHECK_HEADERS([stdlib.h string.h getopt.h])


########################################################################
## netcdf 
## Begin mandatory link test with libnetcdf.a
## see /cvsroot/nco/nco/configure.in,v 1.364 2012/08/03 00:29:44 zender
########################################################################

# Make precious variables for ./configure --help
AC_ARG_VAR(NETCDF_INC,Location of netCDF headers (compile-time))
AC_ARG_VAR(NETCDF_LIB,Location of netCDF library (compile-time))

# User-specified netCDF locations, if any, take precedence
if test "${NETCDF_DIR}" != ''; then
    if test -d "${NETCDF_DIR}"; then
	CPPFLAGS="-I${NETCDF_DIR}/include ${CPPFLAGS}";
	LDFLAGS="-L${NETCDF_DIR}/lib ${LDFLAGS}";
    else
	echo "WARNING: NETCDF_DIR location \"${NETCDF_DIR}\" does not exist!"
    fi
fi # !NETCDF_DIR

if test "${NETCDF_LIB}"; then
    if test -d "${NETCDF_LIB}"; then
	LDFLAGS="-L${NETCDF_LIB} ${LDFLAGS} "
    else
	echo "WARNING: NETCDF_LIB location \"${NETCDF_LIB}\" does not exist!"
    fi
fi # !NETCDF_LIB

if test "${NETCDF_INC}"; then
    if test -d "${NETCDF_INC}"; then
	CPPFLAGS="-I${NETCDF_INC} ${CPPFLAGS} "
    else
	echo "WARNING: NETCDF_INC location \"${NETCDF_INC}\" does not exist!"
    fi
fi # !NETCDF_INC

# no longer support for netCDF2
AC_CHECK_HEADERS([netcdf.h],[],[AC_MSG_FAILURE([cannot find netCDF header])])
AC_DEFINE(NO_NETCDF_2,1,No netCDF version 2.x API)
AC_CHECK_LIB([netcdf],[nc_open],[LIBS="-lnetcdf ${LIBS}"],[AC_MSG_FAILURE([cannot find netCDF library])])

#
#  Test for netCDF4 support
#

#echo "DEBUG: Before netCDF4 section NETCDF_DIR = ${NETCDF_DIR}"
#echo "DEBUG: Before netCDF4 section NETCDF4_DIR = ${NETCDF4_DIR}"
#echo "DEBUG: Before netCDF4 section LDFLAGS = ${LDFLAGS}"
AC_ARG_ENABLE(netcdf4,AS_HELP_STRING([--enable-netcdf4],[Enable netCDF Version 4 features (same as enable-netcdf-4) [[default=yes]]]),enable_netcdf4=${enableval},enable_netcdf4=yes)
dnl # 20100105: make --enable-netcdf-4 synonym for --enable-netcdf4
dnl if test "x${enable_netcdf4}" = 'xno' ; then
dnl     AC_ARG_ENABLE(netcdf-4,AS_HELP_STRING([--enable-netcdf-4],[Enable netCDF Version 4 features (same as enable-netcdf4) [[default=yes]]]),enable_netcdf4=${enableval},enable_netcdf4=yes)  
dnl fi # !netCDF4

AC_ARG_VAR([NETCDF4_DIR],[Root of netCDF4 installation] (compile-time))
# Test if user set NETCDF4_DIR and it exists
if test "x${enable_netcdf4}" = 'xyes'; then
#echo "DEBUG: IN netCDF4 enable_netcdf4 WAS set"
    if test -z "${NETCDF_DIR}"; then
	NETCDF_CONFIG='nc-config'
	NETCDF_DIR="`nc-config --prefix`"
 #       echo "DEBUG: no NETCDF_DIR setting to = ${NETCDF_DIR}"
    else
	NETCDF_CONFIG="${NETCDF_DIR}/bin/nc-config"
 #       echo "DEBUG: yes NETCDF_DIR = ${NETCDF_DIR}"
    fi
    LDFLAGS="`${NETCDF_CONFIG} --libs` ${LDFLAGS}";
    CPPFLAGS="`${NETCDF_CONFIG} --cflags` ${CPPFLAGS}";
#else
#echo "DEBUG: IN netCDF4 enable_netcdf4 was NOT set"
fi
#echo "DEBUG: After netCDF4 section NETCDF_DIR = ${NETCDF_DIR}"
#echo "DEBUG: After netCDF4 section NETCDF4_DIR = ${NETCDF4_DIR}"
#echo "DEBUG: After netCDF4 section LDFLAGS = ${LDFLAGS}"

# ENABLE_NETCDF4 refers to presence of netCDF4-enabled library
# ENABLE_NETCDF4 implies HAVE_NETCDF4_H
# The converse is not true! Some distributions (e.g., Debian stable)
# build _only_ netCDF3 library from netCDF4 source code
AC_CHECK_DECL([NC_CHUNKED],[header_defines_nc_chunked='yes'],[header_defines_nc_chunked='no'],[[#include <netcdf.h>]])
if test "x${header_defines_nc_chunked}" = 'xyes' ; then
    AC_DEFINE([HAVE_NETCDF4_H],[1],[Define to 1 if netcdf.h contains netCDF4 definitions])
    AC_MSG_NOTICE([netcdf.h contains netCDF4 definitions])
fi
# echo "DEBUG: After netCDF4 section header_defines_nc_chunked = ${header_defines_nc_chunked}"
if (test "x${enable_netcdf4}" = 'xyes' && test "x${header_defines_nc_chunked}" = 'xyes') ; then
    AC_CHECK_LIB([netcdf],[nc_def_var_chunking],[LIBS="-lnetcdf ${LIBS}"
	    AC_DEFINE([ENABLE_NETCDF4],[1],[Define to 1 if libnetcdf.a contains netCDF4 functions])
	    AC_MSG_NOTICE([Good news: Simple program compiles and links to netCDF4 library])
	    ],[
	    AC_MSG_NOTICE([Bad news: Simple program does not compile and link with netCDF4 library])
	    ])
fi



dnl
dnl Checks for typedefs, structures, and compiler characteristics.
dnl


AC_CONFIG_FILES([Makefile \
        m4/Makefile \
			  src/Makefile \
			  example/Makefile \
			   ])
dnl data/Makefile
AC_OUTPUT



echo "*******************************************************************"
echo "Configuration complete. Ensure \"make\" points to GNU Make."
echo ""
echo "You are now ready to compile grd2stream for $host"
echo "Type \"make\" to compile grd2stream"
echo "Type \"make install\" to install in: $prefix"
echo "*******************************************************************"
