########################################################################
##
## Init
##
########################################################################
AC_INIT([grd2stream],[0.2.14],[tkleiner@awi.de])
AC_CONFIG_SRCDIR([configure.ac])
# Put helper scripts here to reduce top level clutter
AC_CONFIG_AUX_DIR([admin])
# Tell system about local m4 macros to keep libtool in sync
AC_CONFIG_MACRO_DIR([m4])
# Put preprocessor symbols in config.h
AC_CONFIG_HEADERS([config.h])
# Use automake to create Makefiles
AM_INIT_AUTOMAKE([tar-pax])


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
AC_ARG_ENABLE(debug, AS_HELP_STRING([--enable-debug],[Compile in debug mode @<:@default=no@:>@]),
    enable_debug=$enableval,
    enable_debug=no)
AM_CONDITIONAL(ENABLE_DEBUG, test "x${enable_debug}" == "xyes" )

dnl this could dbe removed later, if api calls are tested
AC_ARG_ENABLE(gmt-api,
        AS_HELP_STRING([--enable-gmt-api],[try linking to gmt api@<:@default=no@:>@]),
        [enable_gmt_api=$enableval],
        [enable_gmt_api=no])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O
LT_INIT
AC_PROG_INSTALL

# netcdf.h and gmt_resources.h require C99 extensions
AX_CHECK_COMPILE_FLAG([-std=c11],
  [AX_APPEND_FLAG([-std=c11])],
  [AC_MSG_WARN([-std=c11 not supported, you may have to set CFLAGS to enable C11 support.])
])

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
# getopt() is in unistd.h
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])


#-------------------------------------------------------------------------------#
#  Find NetCDF library and header files                                         #
#-------------------------------------------------------------------------------#
AC_ARG_WITH([nc-config],
    AS_HELP_STRING([--without-nc-config],
        [do not use nc-config to get netcdf configuration]),
    [],
    [with_nc_config=yes])

AS_IF([test "x$with_nc_config" != xno],
      [AC_CHECK_PROG(have_nc_config, nc-config, yes, no, [], [])]
)

AS_IF([test "x$have_nc_config" == xyes],
  [
    # Prepend linker flags to LDFLAGS:
    AC_MSG_CHECKING(netcdf linker flags)
    NETCDF_LIBS=`nc-config --libs`
    AC_MSG_RESULT($NETCDF_LIBS)
    LDFLAGS="$NETCDF_LIBS $LDFLAGS"
    # Prepend preprocessor and compiler flags to CPPFLAGS:
    AC_MSG_CHECKING(netcdf preprocessor and compiler flags)
    NETCDF_CFLAGS=`nc-config --cflags`
    AC_MSG_RESULT($NETCDF_CFLAGS)
    CPPFLAGS="$NETCDF_CFLAGS $CPPFLAGS"
  ]
)

# Check that netcdf header files can be compiled:
AC_CHECK_HEADERS(netcdf.h, [],
    AC_MSG_ERROR([netcdf.h was not compiled - defining CPPFLAGS may help]))

# Check that netcdf library can be found.
# Linker flags are prepended to LIBS if needed.
AC_SEARCH_LIBS(nc_open, netcdf, [],
    AC_MSG_ERROR([netcdf library was not linked - defining LDFLAGS may help]))




#-------------------------------------------------------------------------------#
#  Find GMT5/6 library and header files                                         #
#-------------------------------------------------------------------------------#
AC_ARG_WITH([gmt-config],
    AS_HELP_STRING([--without-gmt-config],
        [do not use gmt-config to get GMT API configuration]),
    [],
    [with_gmt_config=yes])

AS_IF([test "x$with_gmt_config" != xno],
      [AC_CHECK_PROG(have_gmt_config, gmt-config, yes, no, [], [])]
)

AS_IF([test "x$have_gmt_config" == xyes],
  [
    # Prepend linker flags to LDFLAGS:
    AC_MSG_CHECKING(GMT linker flags)
    GMT_LIBS=`gmt-config --libs`
    AC_MSG_RESULT($GMT_LIBS)
    LDFLAGS="$GMT_LIBS $LDFLAGS"
    # Prepend preprocessor and compiler flags to CPPFLAGS:
    AC_MSG_CHECKING(GMT preprocessor and compiler flags)
    GMT_CFLAGS=`gmt-config --cflags`
    AC_MSG_RESULT($GMT_CFLAGS)
    CPPFLAGS="$GMT_CFLAGS $CPPFLAGS"

    
    gmt_major_version=`gmt-config --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gmt_major_version}" = "x" ; then
       gmt_major_version=0
    fi

    gmt_minor_version=`gmt-config --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gmt_minor_version}" = "x" ; then
       gmt_minor_version=0
    fi

    gmt_micro_version=`gmt-config --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gmt_micro_version}" = "x" ; then
       gmt_micro_version=0
    fi


  ]
)

# Check that gmt header files can be compiled:
AC_CHECK_HEADER([gmt.h],AC_DEFINE([HAVE_GMT_H],1,[Define to 1 if you have the <gmt.h> header file.]),enable_gmt=no)

# Check that gmt library can be found.
AC_SEARCH_LIBS(GMT_Create_Session, gmt, [],enable_gmt=no)

if test "x${enable_gmt_api}" = "xyes"; then
   AC_MSG_NOTICE([GMT support requested by user)])
   if test "${enable_gmt}" = 'no'; then
      AC_MSG_WARN([GMT support disabled (missing/unusable library or header file)])
   else
      if test "x${gmt_major_version}" = "x" ; then
         AC_MSG_WARN([GMT support disabled (missing/unusable library or header file)])
      else 
         CPPFLAGS="$CPPFLAGS -DENABLE_GMT_API=${gmt_major_version}"
      fi
   fi                                    
fi

dnl
dnl check for netcdf only if we don't use the gmt api
dnl





dnl
dnl Checks for typedefs, structures, and compiler characteristics.
dnl

AC_CONFIG_FILES([Makefile \
        m4/Makefile \
	    src/Makefile \
		example/Makefile \
		])
AC_OUTPUT



echo "*******************************************************************"
echo "Configuration complete. Ensure \"make\" points to GNU Make."
echo ""
echo "You are now ready to compile grd2stream for $host"
echo "Type \"make\" to compile grd2stream"
echo "Type \"make install\" to install in: $prefix"
echo "*******************************************************************"
