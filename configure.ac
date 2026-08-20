dnl -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(UWPKFCVM, 1.1.0, software@scec.org)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign tar-pax no-exeext])

AC_PROG_MKDIR_P
LT_INIT
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MKDIR_P

dnl TAR
AC_PATH_PROG(TAR, tar)
if test -z "$TAR" ; then
  AC_MSG_FAILURE([cannot find 'tar' program.])
  TAR=`echo "Error: tar is not installed." ; false`
fi

dnl
dnl Dependencies
dnl

dnl SQLITE
AC_ARG_ENABLE([sqlite],
    [AC_HELP_STRING([--enable-sqlite],
        [install Sqlite @<:@default=no@:>@])],
	[if test "$enableval" = "yes" ; then install_sqlite=yes; else install_sqlite=no; fi],
	[install_sqlite=no])
AM_CONDITIONAL([INSTALL_SQLITE], [test "$install_sqlite" = yes])
AC_ARG_WITH([sqlite-incdir],
    [AC_HELP_STRING([--with-sqlite-incdir],
        [location of sqlite header files @<:@default=no@:>@])],
	[with_sqlite_incdir=$withval],
	[with_sqlite_incdir=no])
AC_SUBST(with_sqlite_incdir)
AC_ARG_WITH([sqlite-libdir],
    [AC_HELP_STRING([--with-sqlite-libdir],
        [location of sqlite library @<:@default=no@:>@])],
	[with_sqlite_libdir=$withval],
	[with_sqlite_libdir=no])
AC_SUBST(with_sqlite_libdir)

dnl TIFF
AC_ARG_ENABLE([tiff],
    [AC_HELP_STRING([--enable-tiff],
        [install Tiff @<:@default=no@:>@])],
	[if test "$enableval" = "yes" ; then install_tiff=yes; else install_tiff=no; fi],
	[install_tiff=no])
AM_CONDITIONAL([INSTALL_TIFF], [test "$install_tiff" = yes])
AC_ARG_WITH([tiff-incdir],
    [AC_HELP_STRING([--with-tiff-incdir],
        [location of tiff header files @<:@default=no@:>@])],
	[with_tiff_incdir=$withval],
	[with_tiff_incdir=no])
AC_SUBST(with_tiff_incdir)
AC_ARG_WITH([tiff-libdir],
    [AC_HELP_STRING([--with-tiff-libdir],
        [location of tiff library @<:@default=no@:>@])],
	[with_tiff_libdir=$withval],
	[with_tiff_libdir=no])
AC_SUBST(with_tiff_libdir)

dnl PROJ
AC_ARG_ENABLE([proj],
    [AC_HELP_STRING([--enable-proj],
        [install Proj @<:@default=yes@:>@])],
	[if test "$enableval" = "yes" ; then install_proj=yes; else install_proj=no; fi],
	[install_proj=yes])
AM_CONDITIONAL([INSTALL_PROJ], [test "$install_proj" = yes])
AC_ARG_WITH([proj-incdir],
    [AC_HELP_STRING([--with-proj-incdir], [directory containing Proj header files @<:@default=no@:>@])],
	[with_proj_incdir="$withval"],
	[with_proj_incdir=no])
AC_SUBST(with_proj_incdir)
AC_ARG_WITH([proj-libdir],
    [AC_HELP_STRING([--with-proj-libdir], [directory containing Proj library @<:@default=no@:>@])],
	[with_proj_libdir="$withval"],
	[with_proj_libdir=no])
AC_SUBST(with_proj_libdir)



dnl
dnl Setup environment so dependencies are used in build
dnl


CFLAGS="$CFLAGS"
CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib -L$prefix/lib64 $LDFLAGS"

dnl
dnl Verify configuration
dnl

dnl proj (required)
if test "$install_proj" = yes ; then
  dnl sqlite (required by proj)
  if test "$install_sqlite" = yes ; then
    SQLITE3_INCLUDES="-I$prefix/include"
    SQLITE3_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_sqlite_incdir" != no; then
      SQLITE3_INCLUDES="-I$with_sqlite_incdir"
    fi
    if test "$with_sqlite_libdir" != no; then
      SQLITE3_LDFLAGS="-L$with_sqlite_libdir"
    fi
  fi
  AC_SUBST(SQLITE3_INCLUDES)
  AC_SUBST(SQLITE3_LDFLAGS)

  dnl tiff (required by proj)
  if test "$install_tiff" = yes ; then
    TIFF_INCLUDES="-I$prefix/include"
    TIFF_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_tiff_incdir" != no; then
      TIFF_INCLUDES="-I$with_tiff_incdir"
    fi
    if test "$with_tiff_libdir" != no; then
      TIFF_LDFLAGS="-L$with_tiff_libdir"
    fi
  fi
  AC_SUBST(TIFF_INCLUDES)
  AC_SUBST(TIFF_LDFLAGS)

  PROJ_INCLUDES="-I$prefix/include"
  PROJ_LDFLAGS="-L$prefix/lib"
else
  SCEC_PROJ_HEADER
  SCEC_PROJ_LIB
fi
if test "$with_proj_incdir" != no ; then
  PROJ_INCLUDES="-I$with_proj_incdir"
fi
if test "$with_proj_libdir" != no ; then
  PROJ_LDFLAGS="-L$with_proj_libdir -lproj"
fi
AC_SUBST(PROJ_INCLUDES)
AC_SUBST(PROJ_LDFLAGS)



##check optional large data path
##UWPKFCVM_LARGEDATA_DIR=$CVM_LARGEDATA_DIR/model/uwpkfcvm
if test x"$CVM_LARGEDATA_DIR" != x; then
   # test directory existence
   UWPKFCVM_LARGEDATA_DIR=$CVM_LARGEDATA_DIR/model/uwpkfcvm
# In docker container building.. this is nto accessible  yet
   if test x"$CVM_IN_DOCKER" != x; then
     AM_CONDITIONAL([WITH_UWPKFCVM_LARGEDATA_DIR], true)
     AC_SUBST(UWPKFCVM_LARGEDATA_DIR)
   else
     AC_CHECK_FILE([$UWPKFCVM_LARGEDATA_DIR/vp.dat], 
                 [AM_CONDITIONAL([WITH_UWPKFCVM_LARGEDATA_DIR], true) AC_SUBST(UWPKFCVM_LARGEDATA_DIR)],
                 [AM_CONDITIONAL(WITH_UWPKFCVM_LARGEDATA_DIR, false)])
   fi
else
   AM_CONDITIONAL(WITH_UWPKFCVM_LARGEDATA_DIR, false)
fi

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  data/Makefile
  test/Makefile
  ])

AC_OUTPUT

AM_CONDITIONAL([BUILD_SHARED_LIBRARY], [test "$enable_shared" = yes])
if test "$enable_shared" = yes; then
  AC_MSG_NOTICE([Building shared libraries. Shared libraries enabled.])
else
  AC_MSG_WARN([Building static libraries. Shared libraries are not enabled.])
fi


dnl End of file
