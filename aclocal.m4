dnl aclocal.m4 generated automatically by aclocal 1.4

dnl Copyright (C) 1994, 1995-8, 1999 Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY, to the extent permitted by law; without
dnl even the implied warranty of MERCHANTABILITY or FITNESS FOR A
dnl PARTICULAR PURPOSE.

# Do all the work for Automake.  This macro actually does too much --
# some checks are only needed if your package does certain things.
# But this isn't really a big deal.

# serial 1

dnl Usage:
dnl AM_INIT_AUTOMAKE(package,version, [no-define])

AC_DEFUN(AM_INIT_AUTOMAKE,
[AC_REQUIRE([AC_PROG_INSTALL])
PACKAGE=[$1]
AC_SUBST(PACKAGE)
VERSION=[$2]
AC_SUBST(VERSION)
dnl test to see if srcdir already configured
if test "`cd $srcdir && pwd`" != "`pwd`" && test -f $srcdir/config.status; then
  AC_MSG_ERROR([source directory already configured; run "make distclean" there first])
fi
ifelse([$3],,
AC_DEFINE_UNQUOTED(PACKAGE, "$PACKAGE", [Name of package])
AC_DEFINE_UNQUOTED(VERSION, "$VERSION", [Version number of package]))
AC_REQUIRE([AM_SANITY_CHECK])
AC_REQUIRE([AC_ARG_PROGRAM])
dnl FIXME This is truly gross.
missing_dir=`cd $ac_aux_dir && pwd`
AM_MISSING_PROG(ACLOCAL, aclocal, $missing_dir)
AM_MISSING_PROG(AUTOCONF, autoconf, $missing_dir)
AM_MISSING_PROG(AUTOMAKE, automake, $missing_dir)
AM_MISSING_PROG(AUTOHEADER, autoheader, $missing_dir)
AM_MISSING_PROG(MAKEINFO, makeinfo, $missing_dir)
AC_REQUIRE([AC_PROG_MAKE_SET])])

#
# Check to make sure that the build environment is sane.
#

AC_DEFUN(AM_SANITY_CHECK,
[AC_MSG_CHECKING([whether build environment is sane])
# Just in case
sleep 1
echo timestamp > conftestfile
# Do `set' in a subshell so we don't clobber the current shell's
# arguments.  Must try -L first in case configure is actually a
# symlink; some systems play weird games with the mod time of symlinks
# (eg FreeBSD returns the mod time of the symlink's containing
# directory).
if (
   set X `ls -Lt $srcdir/configure conftestfile 2> /dev/null`
   if test "[$]*" = "X"; then
      # -L didn't work.
      set X `ls -t $srcdir/configure conftestfile`
   fi
   if test "[$]*" != "X $srcdir/configure conftestfile" \
      && test "[$]*" != "X conftestfile $srcdir/configure"; then

      # If neither matched, then we have a broken ls.  This can happen
      # if, for instance, CONFIG_SHELL is bash and it inherits a
      # broken ls alias from the environment.  This has actually
      # happened.  Such a system could not be considered "sane".
      AC_MSG_ERROR([ls -t appears to fail.  Make sure there is not a broken
alias in your environment])
   fi

   test "[$]2" = conftestfile
   )
then
   # Ok.
   :
else
   AC_MSG_ERROR([newly created file is older than distributed files!
Check your system clock])
fi
rm -f conftest*
AC_MSG_RESULT(yes)])

dnl AM_MISSING_PROG(NAME, PROGRAM, DIRECTORY)
dnl The program must properly implement --version.
AC_DEFUN(AM_MISSING_PROG,
[AC_MSG_CHECKING(for working $2)
# Run test in a subshell; some versions of sh will print an error if
# an executable is not found, even if stderr is redirected.
# Redirect stdin to placate older versions of autoconf.  Sigh.
if ($2 --version) < /dev/null > /dev/null 2>&1; then
   $1=$2
   AC_MSG_RESULT(found)
else
   $1="$3/missing $2"
   AC_MSG_RESULT(missing)
fi
AC_SUBST($1)])

# Like AC_CONFIG_HEADER, but automatically create stamp file.

AC_DEFUN(AM_CONFIG_HEADER,
[AC_PREREQ([2.12])
AC_CONFIG_HEADER([$1])
dnl When config.status generates a header, we must update the stamp-h file.
dnl This file resides in the same directory as the config header
dnl that is generated.  We must strip everything past the first ":",
dnl and everything past the last "/".
AC_OUTPUT_COMMANDS(changequote(<<,>>)dnl
ifelse(patsubst(<<$1>>, <<[^ ]>>, <<>>), <<>>,
<<test -z "<<$>>CONFIG_HEADERS" || echo timestamp > patsubst(<<$1>>, <<^\([^:]*/\)?.*>>, <<\1>>)stamp-h<<>>dnl>>,
<<am_indx=1
for am_file in <<$1>>; do
  case " <<$>>CONFIG_HEADERS " in
  *" <<$>>am_file "*<<)>>
    echo timestamp > `echo <<$>>am_file | sed -e 's%:.*%%' -e 's%[^/]*$%%'`stamp-h$am_indx
    ;;
  esac
  am_indx=`expr "<<$>>am_indx" + 1`
done<<>>dnl>>)
changequote([,]))])

dnl ####################### -*- Mode: M4 -*- ###########################
dnl smr.m4 -- 
dnl 
dnl Copyright (C) 1999 Matthew D. Langston <langston@SLAC.Stanford.EDU>
dnl Copyright (C) 1998 Steve Robbins <stever@cs.mcgill.ca>
dnl
dnl This file is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This file is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this file; if not, write to:
dnl
dnl   Free Software Foundation, Inc.
dnl   Suite 330
dnl   59 Temple Place
dnl   Boston, MA 02111-1307, USA.
dnl ####################################################################


dnl NOTE: The macros in this file are extensively documented in the
dnl       accompanying `smr_macros.texi' Texinfo file.  Please see the
dnl       Texinfo documentation for the definitive specification of how
dnl       these macros are supposed to work.  If the macros work
dnl       differently than the Texinfo documentation says they should,
dnl       then the macros (and not the Texinfo documentation) has the
dnl       bug(s).

dnl This is a convenient macro which translates illegal characters for
dnl bourne shell variables into legal characters.  It has the same
dnl functionality as sed 'y%./+-:%__p__%'.
AC_DEFUN([smr_safe_translation], [patsubst(patsubst([$1], [+], [p]), [./-:], [_])])

AC_DEFUN(smr_SWITCH,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_help_string], $2)
  pushdef([smr_default],     $3)
  pushdef([smr_yes_define],  $4)
  pushdef([smr_no_define],   $5)

  dnl Do some sanity checking of the arguments.
  ifelse([regexp(smr_default, [^\(yes\|no\)$])], -1, [AC_MSG_ERROR($0: third arg must be either yes or no)])

  dnl Create the help string
  pushdef([smr_lhs], [--ifelse(smr_default, yes, disable, enable)-smr_name])dnl
  pushdef([smr_rhs], [ifelse(smr_default, yes, disable, enable) smr_help_string (default is smr_default)])dnl

  AC_HELP_STRING([smr_lhs], [smr_rhs], smr_name[]_switch_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_ENABLE (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_ENABLE([smr_name],
                [$]smr_name[]_switch_help_string,
                smr_cv_enable_[]smr_name=$enableval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([whether to enable smr_help_string],
                 smr_cv_enable_[]smr_name,
                 smr_cv_enable_[]smr_name=smr_default)

    ifelse(smr_yes_define, , , test x"[$]smr_cv_enable_[]smr_name" = xyes && AC_DEFINE(smr_yes_define))
    ifelse(smr_no_define, , ,  test x"[$]smr_cv_enable_[]smr_name" = xno  && AC_DEFINE(smr_no_define))

  dnl Sanity check the value assigned to smr_cv_enable_$1 to force it to
  dnl be either `yes' or `no'.
  if test ! x"[$]smr_cv_enable_[]smr_name" = xyes; then
    if test ! x"[$]smr_cv_enable_[]smr_name" = xno; then
      AC_MSG_ERROR([smr_lhs must be either yes or no])
    fi
  fi

  popdef([smr_name])
  popdef([smr_help_string])
  popdef([smr_default])
  popdef([smr_yes_define])
  popdef([smr_no_define])
  popdef([smr_lhs])
  popdef([smr_rhs])
])


AC_DEFUN(smr_ARG_WITHLIB,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_libname],     ifelse($2, , smr_name, $2))
  pushdef([smr_help_string], $3)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Create the help string
  AC_HELP_STRING([--with-smr_safe_name-library[[=PATH]]],
                 [use smr_name library ifelse(smr_help_string, , , (smr_help_string))],
                 smr_safe_name[]_library_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_WITH (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_WITH(smr_safe_name-library,
              [$]smr_safe_name[]_library_help_string,
              smr_cv_with_[]smr_safe_name[]_library=$withval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([whether to use smr_name library],
                 smr_cv_with_[]smr_safe_name[]_library,
                 smr_cv_with_[]smr_safe_name[]_library=maybe)


  case x"[$]smr_cv_with_[]smr_safe_name[]_library" in
      xyes | xmaybe)
          smr_safe_name[]_LIBS="-l[]smr_libname"
          with_[]smr_safe_name=[$]smr_cv_with_[]smr_safe_name[]_library
          ;;
      xno)
          smr_safe_name[]_LIBS=
          with_[]smr_safe_name=no
          ;;
      *)
          if test -f "[$]smr_cv_with_[]smr_safe_name[]_library"; then
            smr_safe_name[]_LIBS=[$]smr_cv_with_[]smr_safe_name[]_library
          elif test -d "[$]smr_cv_with_[]smr_safe_name[]_library"; then
            smr_safe_name[]_LIBS="-L[$]smr_cv_with_[]smr_safe_name[]_library -l[]smr_libname"
          else
            AC_MSG_ERROR([argument must be boolean, file, or directory])
          fi
          with_[]smr_safe_name=yes
          ;;
  esac

  popdef([smr_name])
  popdef([smr_libname])
  popdef([smr_help_string])
  popdef([smr_safe_name])
])


AC_DEFUN(smr_ARG_WITHINCLUDES,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_header],      $2)
  pushdef([smr_extra_flags], $3)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Create the help string
  AC_HELP_STRING([--with-smr_safe_name-includes[[=DIR]]],
                 [set directory for smr_name headers],
                 smr_safe_name[]_includes_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_WITH (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_WITH(smr_safe_name-includes,
              [$]smr_safe_name[]_includes_help_string,
              smr_cv_with_[]smr_safe_name[]_includes=$withval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([where to find the smr_name header files],
                 smr_cv_with_[]smr_safe_name[]_includes,
                 smr_cv_with_[]smr_safe_name[]_includes=)

  if test ! x"[$]smr_cv_with_[]smr_safe_name[]_includes" = x; then
    if test -d "[$]smr_cv_with_[]smr_safe_name[]_includes"; then
      smr_safe_name[]_CFLAGS="-I[$]smr_cv_with_[]smr_safe_name[]_includes"
    else
      AC_MSG_ERROR([argument must be a directory])
    fi
  else
    smr_safe_name[]_CFLAGS=
  fi

  dnl This bit of logic comes from the autoconf AC_PROG_CC macro.  We
  dnl need to put the given include directory into CPPFLAGS temporarily,
  dnl but then restore CPPFLAGS to its old value.
  smr_test_CPPFLAGS="${CPPFLAGS+set}"
  smr_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS [$]smr_safe_name[]_CFLAGS smr_extra_flags"

  dnl If the header file smr_header exists, then define
  dnl HAVE_[]smr_header (in all capitals).
  AC_CHECK_HEADERS([smr_header],
                   smr_have_[]smr_safe_name[]_header=yes,
                   smr_have_[]smr_safe_name[]_header=no)

  if test x"$smr_test_CPPFLAGS" = xset; then
    CPPFLAGS=$smr_save_CPPFLAGS
  else
    unset CPPFLAGS
  fi

  popdef([smr_name])
  popdef([smr_header])
  popdef([smr_extra_flags])
  popdef([smr_safe_name])
])


AC_DEFUN(smr_CHECK_LIB,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_libname],     ifelse($2, , smr_name, $2))
  pushdef([smr_help_string], $3)
  pushdef([smr_function],    $4)
  pushdef([smr_header],      $5)
  pushdef([smr_extra_libs],  $6)
  pushdef([smr_extra_flags], $7)
  pushdef([smr_prototype],   $8)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Give the user (via "configure --help") an interface to specify
  dnl whether we should use the library or not, and possibly where we
  dnl should find it.
  smr_ARG_WITHLIB([smr_name], [smr_libname], [smr_help_string])

  if test ! x"$with_[]smr_safe_name" = xno; then

    # If we got this far, then the user didn't explicitly ask not to use
    # the library.

    dnl If the caller of smr_CHECK_LIB specified a header file for this
    dnl library, then give the user (via "configure --help") an
    dnl interface to specify where this header file can be found (if it
    dnl isn't found by the compiler by default).
    ifelse(smr_header, , , [smr_ARG_WITHINCLUDES(smr_name, smr_header, smr_extra_flags)])

    # We need only look for the library if the header has been found
    # (or no header is needed).
    if test [$]smr_have_[]smr_safe_name[]_header != no; then

       mdl_CHECK_LIB(smr_libname,
                    smr_function,
                    smr_have_[]smr_safe_name[]_library=yes,
                    smr_have_[]smr_safe_name[]_library=no,
                    [$]smr_safe_name[]_CFLAGS [smr_extra_flags] [$]smr_safe_name[]_LIBS [smr_extra_libs],
                    [ifelse(smr_prototype, , , [[#]include <smr_header>])],
                    smr_prototype)
    fi

    if test x"[$]smr_have_[]smr_safe_name[]_library" = xyes; then
      AC_MSG_RESULT([using smr_name library])
    else
      smr_safe_name[]_LIBS=
      smr_safe_name[]_CFLAGS=

      if test x"$with_[]smr_safe_name" = xmaybe; then
        AC_MSG_RESULT([not using smr_name library])
      else
        AC_MSG_WARN([requested smr_name library not found!])
      fi
    fi
  fi

  popdef([smr_name])
  popdef([smr_libname])
  popdef([smr_help_string])
  popdef([smr_function])
  popdef([smr_header])
  popdef([smr_extra_libs])
  popdef([smr_extra_flags])
  popdef([smr_prototype])
  popdef([smr_safe_name])
])


dnl This is a small fix to the standard `AC_CHECK_LIB' macro which
dnl allows checking for C++ libraries.
dnl
dnl AC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, INCLUDES [, FUNC-PTR-PROTOTYPE]]]]])
pushdef([AC_CHECK_LIB],
[AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-:%__p__%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="-l$1 $5 $LIBS"
AC_TRY_LINK(dnl
ifelse([$6], , dnl If there is no INCLUDES then concoct our own declaration.
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])), [$6]),
dnl If the user passed in a prototype for a pointer to FUNCTION
dnl (i.e. $7) then instead of just calling FUNCTION, construct a
dnl statement that takes the address of FUNCTION.  This is particularly
dnl helpful in checking for C++ class libraries.
ifelse([$7], ,[$2()], [$7 = &$2]),
	    eval "ac_cv_lib_$ac_lib_var=yes",
	    eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])




dnl Copy of the above with one small modification.
dnl We suppress adding "-l$1" to the LIBS, since that is already
dnl done by smr_CHECK_LIB.  Yes, this is a gross hack.  I'm counting
dnl on autoconf 2.50 allowing me to clean this up, he says vainly.
dnl
dnl AC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, INCLUDES [, FUNC-PTR-PROTOTYPE]]]]])
AC_DEFUN(mdl_CHECK_LIB,
[AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-:%__p__%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="$5 $LIBS"
AC_TRY_LINK(dnl
ifelse([$6], , dnl If there is no INCLUDES then concoct our own declaration.
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])), [$6]),
dnl If the user passed in a prototype for a pointer to FUNCTION
dnl (i.e. $7) then instead of just calling FUNCTION, construct a
dnl statement that takes the address of FUNCTION.  This is particularly
dnl helpful in checking for C++ class libraries.
ifelse([$7], ,[$2()], [$7 = &$2]),
	    eval "ac_cv_lib_$ac_lib_var=yes",
	    eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])

dnl ####################### -*- Mode: M4 -*- ###########################
dnl Copyright (C) 98, 99, 2000 Matthew D. Langston <langston@SLAC.Stanford.EDU>
dnl
dnl This file is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This file is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this file; if not, write to:
dnl
dnl   Free Software Foundation, Inc.
dnl   Suite 330
dnl   59 Temple Place
dnl   Boston, MA 02111-1307, USA.
dnl ####################################################################


dnl NOTE: The AC_HELP_STRING macro has been accepted for inclusion in
dnl       Autoconf 2.15, but this version of Autoconf hasn't been
dnl       released yet as of the time that I write this.  Therefore, I
dnl       am including it here for the time being


dnl AC_HELP_STRING
dnl --------------
dnl
dnl usage: AC_HELP_STRING(LHS, RHS, HELP-STRING)
dnl
dnl Format an Autoconf macro's help string so that it looks pretty when
dnl the user executes "configure --help".  This macro take three
dnl arguments, a "left hand side" (LHS), a "right hand side" (RHS), and
dnl a variable (HELP-STRING) to set to the pretty-printed concatenation
dnl of LHS and RHS (the new, pretty-printed "help string").
dnl
dnl The resulting string in HELP-STRING is suitable for use in other
dnl macros that require a help string (e.g. AC_ARG_WITH).
dnl 
dnl AC_DEFUN(AC_HELP_STRING,
pushdef([AC_HELP_STRING],
[
dnl 
dnl Here is the sample string from the Autoconf manual (Node: External
dnl Software) which shows the proper spacing for help strings.
dnl 
dnl    --with-readline         support fancy command line editing
dnl  ^ ^                       ^ 
dnl  | |                       |
dnl  | column 2                column 26     
dnl  |
dnl  column 0
dnl 
dnl A help string is made up of a "left hand side" (LHS) and a "right
dnl hand side" (RHS).  In the example above, the LHS is
dnl "--with-readline", while the RHS is "support fancy command line
dnl editing".
dnl 
dnl If the LHS extends past column 24, then the LHS is terminated with a
dnl newline so that the RHS is on a line of its own beginning in column
dnl 26.
dnl 
dnl Therefore, if the LHS were instead "--with-readline-blah-blah-blah",
dnl then the AC_HELP_STRING macro would expand into:
dnl
dnl
dnl    --with-readline-blah-blah-blah
dnl  ^ ^                       support fancy command line editing
dnl  | |                       ^ 
dnl  | column 2                |
dnl  column 0                  column 26     

dnl We divert everything to AC_DIVERSION_NOTICE (which gets output very
dnl early in the configure script) because we want the user's help
dnl string to be set before it is used.

AC_DIVERT_PUSH(AC_DIVERSION_NOTICE)dnl
# This is from AC_HELP_STRING
lhs="$1"
rhs="$2"

lhs_column=25
rhs_column=`expr $lhs_column + 1`

# Insure that the LHS begins with exactly two spaces.
changequote(, )dnl
lhs=`echo "$lhs" | sed -n -e "s/[ ]*\(.*\)/  \1/p"`
changequote([, ])dnl

# Is the length of the LHS less than $lhs_column?
if `echo "$lhs" | grep -v ".\{$lhs_column\}" > /dev/null 2>&1`; then

  # Pad the LHS with spaces.  Note that padding the LHS is an
  # "expensive" operation (i.e. expensive in the sense of there being
  # multiple calls to `grep') only the first time AC_HELP_STRING is
  # called.  Once this macro is called once, subsequent calls will be
  # nice and zippy.
  : ${lhs_pad=""}
changequote(, )dnl
  while `echo "$lhs_pad" | grep -v "[ ]\{$lhs_column\}" > /dev/null 2>&1`; do
changequote([, ])dnl
    lhs_pad=" $lhs_pad"
  done

  lhs="${lhs}${lhs_pad}"
changequote(, )dnl
$3=`echo "$lhs" | sed -n -e "/.\{$lhs_column\}[ ][ ]*$/ s/\(.\{$rhs_column\}\).*/\1$rhs/p"`
changequote([, ])dnl

else

  # Build up a string of spaces to pad the left-hand-side of the RHS
  # with.  Note that padding the RHS is an "expensive" operation
  # (i.e. expensive in the sense of there being multiple calls to
  # `grep') only the first time AC_HELP_STRING is called.  Once this
  # macro is called once, subsequent calls will be nice and zippy.
  : ${rhs_pad=""}
changequote(, )dnl
  while `echo "$rhs_pad" | grep -v "[ ]\{$rhs_column\}" > /dev/null 2>&1`; do
changequote([, ])dnl
    rhs_pad=" $rhs_pad"
  done

  # Strip all leading spaces from the RHS.
changequote(, )dnl
  rhs=`echo "$rhs" | sed -n -e "s/[ ]*\(.*\)/\1/p"`
changequote([, ])dnl

$3="$lhs
${rhs_pad}${rhs}"
fi 
AC_DIVERT_POP()dnl
])

