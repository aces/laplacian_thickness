dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([solveLaplace.cc])
AM_INIT_AUTOMAKE(laplacian_thickness, 0.1)
AM_CONFIG_HEADER(config.h)

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

dnl Build only static libs by default

smr_REQUIRED_LIB(m,sqrt,math.h)
smr_REQUIRED_LIB(netcdf,ncopen,netcdf.h)
smr_REQUIRED_LIB(minc,miicv_create,minc.h)
smr_REQUIRED_LIB(volume_io,input_volume,volume_io.h)
smr_REQUIRED_LIB(bicpl,input_graphics_file,bicpl.h)

AC_LANG_PUSH(C++)

AC_MSG_CHECKING([for -loobicpl])
LIBS="-loobicpl $LIBS"
AC_TRY_LINK([#include <mniVolume.h>],
            [mniVolume::mniVolume();]
            AC_MSG_ERROR(cannot find oobicpl library))
AC_MSG_RESULT([yes])

AC_LANG_POP

AC_CONFIG_FILES([Makefile])
AC_OUTPUT