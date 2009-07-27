dnl Process this file with autoconf to produce a configure script.

# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([laplacian_thickness],[1.2.0],[Claude Lepage <claude@bic.mni.mcgill.ca>])
AC_CONFIG_SRCDIR([src/laplacianThickness.cc])

AM_INIT_AUTOMAKE

AC_CONFIG_HEADER([config.h])

smr_WITH_BUILD_PATH

AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_CC
AC_PROG_CXX

AC_PATH_PROGS(PERL, [perl], no)

dnl Build only static libs by default

mni_REQUIRE_BICPL

AC_LANG_PUSH(C++)

AC_MSG_CHECKING([for -loobicpl])
LIBS="-loobicpl $LIBS"
AC_TRY_LINK([#include <mniVolume.h>],
            [mniVolume::mniVolume();]
            AC_MSG_ERROR(cannot find oobicpl library))
AC_MSG_RESULT([yes])

AC_LANG_POP

AC_CONFIG_FILES([Makefile src/Makefile scripts/Makefile linear_thickness/Makefile])
AC_OUTPUT(
scripts/laplacian_initialisation
epm-header)
