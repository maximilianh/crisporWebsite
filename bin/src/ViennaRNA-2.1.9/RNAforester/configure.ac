dnl Process this file with autoconf to produce a configure script.
AC_REVISION($Revision: 1.13 $)
AC_INIT(RNAforester,1.5,mhoechsm@techfak.uni-bielefeld.de,RNAforester)
AC_PREREQ([2.59])

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.11])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER(config.h)
		
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)

dnl Checks for libraries.
dnl AC_CHECK_LIB(g2,g2_open_vd,,AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))

AC_OPENMP

AC_ARG_WITH(ViennaRNA-lib,
  [--with-ViennaRNA-lib        where to search for libRNA.a],
  [LIBS=["$withval/libRNA.a ${LIBS}"] && AC_DEFINE(HAVE_LIBRNA,1,Define to 1 if you have the `RNA' library (-lRNA).)
  LDFLAGS="-L$withval -R$withval ${LDFLAGS}"],
  LDFLAGS="-L../../lib -L../../ViennaRNA ${LDFLAGS} ${OPENMP_CXXFLAGS}" 
  if test -r "../lib/fold.c"; then
    AC_DEFINE(HAVE_LIBRNA)
    LIBS=["-lRNA ${LIBS}"]
    AC_MSG_RESULT([using (to be built) libRNA.a in ../lib])
  else
    if test -r "../ViennaRNA/fold.c"; then
      AC_DEFINE(HAVE_LIBRNA)
      LIBS=["-lRNA ${LIBS}"]
      AC_MSG_RESULT([using (to be built) libRNA.a in ../ViennaRNA])
    fi
  fi
 ) 
dnl AC_CHECK_LIB(RNA, fold,,AC_MSG_WARN(Vienna RNA library libRNA.a was not found. RNAforester will be compiled without the structure prediction option -p and graphical output))

AC_ARG_WITH([xml2-lib],
            [AS_HELP_STRING([--with-xml2-lib],[where to search for libxml2.so])],
	    [LDFLAGS="-L$withval/lib -R$withval/lib ${LDFLAGS}"
	    AC_CHECK_LIB([xml2],[main],AC_DEFINE([HAVE_LIBXML2],[1],[Define if you have libxml2]),AC_MSG_FAILURE([--with-xml2-lib was given but test for libxml2 failed]))
	    LIBS="-lxml2 ${LIBS}"
	    CPPFLAGS="-I$withval/include/libxml2/ ${CPPFLAGS}"],
	    [])

AC_ARG_WITH([xmlplusplus-lib],
            [AS_HELP_STRING([--with-xmlplusplus-lib],[where to search for libxml++-1.0.so])],
	    [LDFLAGS="-L$withval/lib -R$withval/lib ${LDFLAGS}"
	    AC_CHECK_LIB([xml++-1.0],[main],AC_DEFINE([HAVE_LIBXMLPLUSPLUS],[1],[Define if you have libxml++-1.0]),[AC_MSG_FAILURE([--with-xmlplusplus-lib was given but test for libxml++-1.0 failed])],[])
	    LIBS="-lxml++-1.0 ${LIBS}"
	    CPPFLAGS="-I$withval/include/libxml++-1.0/ ${CPPFLAGS}"],
	    [])

save_CPPFLAGS=$CPPFLAGS	
CPPFLAGS=["-Ig2-0.70/src/ ${CPPFLAGS}"]

dnl AC_CHECK_LIB(gd,main,,)
		  
dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADER(g2.h,AC_DEFINE(HAVE_LIBG2,1,Define to 1 if you have the `g2' library (-lg2).),AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))
dnl AC_CHECK_HEADER(g2_PS.h,,AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

CPPFLAGS="$save_CPPFLAGS -I../g2-0.70/include/ "
LDFLAGS=["-L../g2-0.70 ${LDFLAGS}"]

if ${HAVE_LIBG2}; then
  LIBS=["-lg2 ${LIBS}"]
fi	
											
AC_SUBST(CXX)
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)

dnl Checks for library functions.

AC_CONFIG_SUBDIRS(g2-0.70)
AC_OUTPUT(Makefile src/Makefile man/Makefile)














