#!/bin/sh -x

# TODO:
# 1. clean all CVS files
# 2. clean all Makefiles from makedepend-entry

CleanMakedepend() {
    # usage: $0 makefile
    dir=`dirname $1`
    cp $1 $dir/Makefile.bck
    awk 'BEGIN {nodepend=1;}
/^# DO NOT DELETE/ { nodepend=0; }
/a*/ { if (nodepend) print; }' $dir/Makefile.bck > $1
    rm -f $dir/Makefile.bck
}
    
if test $# -lt 2 ; then
    echo "
  Usage: $0 topdir in-tarball [suffix]
"
    exit 1
fi

TOPDIR=$1
INTAR=$2
SUFFIX=$3
VERSION=`cat $TOPDIR/version`

if test x$SUFFIX != x; then
    DIRNAME=xcrysden-${VERSION}-${SUFFIX}
else
    DIRNAME=xcrysden-${VERSION}
fi

if test ! -d $TOPDIR ; then
    echo "TOPDIR directory \"$TOPDIR\" does not exists !!!"
    exit 1
fi
if test ! -f $INTAR ; then
    echo "tar-ball file \"$INTAR\" does not exists !!!"
    exit 1
fi


cd $TOPDIR
mkdir $DIRNAME
if test ! -d $DIRNAME ; then
    echo "DIRNAME directory \"$DIRNAME\" does not exists !!!"
    exit 1
fi
mv $INTAR $DIRNAME/
cd $DIRNAME
tar xvf $INTAR; rm -f $INTAR

# clean CVS files
find . -name CVS -exec /bin/rm -r {} \;

# clean all Makefiles from makedepend-entry
for make in `find . | grep Makefile`; do
    CleanMakedepend $make
done

cd ..
tar cvf $DIRNAME.tar $DIRNAME/
gzip $DIRNAME.tar
rm -rf $DIRNAME/
