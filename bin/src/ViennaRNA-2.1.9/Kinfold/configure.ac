#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

AC_INIT([kinfold], [1.3], [rna@tbi.univie.ac.at], [Kinfold])
AC_CONFIG_SRCDIR([cache_util.h])

AM_SILENT_RULES([yes])
dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE([1.11])

dnl configure options
AM_WITH_DMALLOC

dnl Checks for programs.
AC_PROG_CC
dnl AC_PROG_MAKE_SET

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS(config.h)

AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl Checks for library functions.
AC_CHECK_FUNCS([strdup memset strchr])

AC_ARG_WITH(ViennaRNA-include,
  [ --with-ViennaRNA-include	where to search for ViennaRNA header files],
  ac_VRNA_includes="$withval",
  )

AC_ARG_WITH(ViennaRNA-lib,
  [ --with-ViennaRNA-lib	where to search for libRNA.a],
  ac_VRNA_lib="$withval",
  )
AC_OPENMP

AC_PATH_VRNA

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(stdlib.h string.h strings.h unistd.h fold.h)

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_CHECK_LIB(m, exp)
AC_FUNC_MALLOC
AC_FUNC_REALLOC

dnl checking for -lRNA will fail if both are built simultaneously
dnl AC_CHECK_LIB(RNA, fold)
LIBS=["-lRNA ${LIBS}"]

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

AC_CONFIG_FILES([Example/Makefile Makefile])

AC_OUTPUT
