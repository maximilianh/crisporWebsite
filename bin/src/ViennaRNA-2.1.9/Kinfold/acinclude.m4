dnl -*-autoconf-*-

AC_DEFUN([AC_PATH_VRNA],
[AC_MSG_CHECKING([for ViennaRNA package])
if test -z "$ac_VRNA_includes"; then
  for ac_dir in               \
    ../H                      \
    /usr/local/include/ViennaRNA \
    /usr/local/include        \
    /usr/include/ViennaRNA    \
    /usr/local/ViennaRNA/H    \
    /usr/local/share/ViennaRNA/include \
    /opt/ViennaRNA/include \
    ;\
  do
    if test -r "$ac_dir/part_func.h"; then
      ac_VRNA_includes=$ac_dir
      break
    fi
  done
fi
if test $ac_VRNA_includes; then
  CPPFLAGS="$CPPFLAGS -I$ac_VRNA_includes"
fi

if [[ -d ../lib ] && [ "$ac_VRNA_includes" = "../H" ]]; then
  ac_VRNA_lib=../lib
fi

if test -z "$ac_VRNA_lib"; then
for ac_dir in `echo "$ac_VRNA_includes" | sed -e s/include/lib/ -e s/H$/lib/` \
   /usr/local/lib \
   ; \
do
  for ac_extension in a so sl; do
    if test -r $ac_dir/libRNA.$ac_extension; then
      ac_VRNA_lib=$ac_dir
      break 2
    fi
  done
done
fi # $ac_VRNA_lib = NO
if test $ac_VRNA_lib; then
  LDFLAGS="-L$ac_VRNA_lib $LDFLAGS"
fi
AC_MSG_RESULT([ headers in "$ac_VRNA_includes" and library... "$ac_VRNA_lib"])

dnl So far we've only set up paths, we could also check for
dnl usability of headers and library like so
dnl AC_CHECK_HEADER(part_func.h, [],
dnl		[AC_MSG_ERROR([Cannot find ViennaRNA headers])])
dnl AC_CHECK_LIB(RNA, fold)
])
