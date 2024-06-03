#!/bin/sh

InstallError () {
    echo "
*************************************************
An ERROR occured during configuration of XCrySDen

Code: $1
Exit 1
"
    exit 1
}

#
# Usage: xcConfigure [topdir]
#

if test "$#" -lt 1 ; then

    if test "$0" = "./xcConfigure.sh"; then
	    XCRYSDEN_TOPDIR=`cd ..; pwd`
	    export XCRYSDEN_TOPDIR
    elif test "$0" = "./scripts/xcConfigure.sh"; then
	    XCRYSDEN_TOPDIR=`pwd`
	    export XCRYSDEN_TOPDIR
    else
	echo "
Please goto XCrySDen toplevel (root) directory and execute the xcConfigure.sh
as \"./scripts/xcConfigure.sh\" or from any other place as \"xcConfigure.sh topdir\"
"
	exit 1
    fi
else
    if test \( -d "$1" \) -a \( -x $1/xcrysden\); then
	XCRYSDEN_TOPDIR=$1
	export XCRYSDEN_TOPDIR
    else
	InstallError "XCRYSDEN_TOPDIR $1 does not exists"
    fi
fi


# check if XCRYSDEN_TOPDIR exists:
xcv=`env | grep XCRYSDEN_TOPDIR`

cd $XCRYSDEN_TOPDIR

#
# assuming user shell is csh or tcsh or sh
#

if test "`echo -n a`" = "-n a" ; then
    ECHO_n=echo
else
    ECHO_n="echo -n"
fi


if [ \( -d /tmp \) -a \( -w /tmp \) ]; then
    tempdir=/tmp
else
    # set tempdir to $HOME
    tempdir=$HOME
fi

CLEAR=$XCRYSDEN_TOPDIR/scripts/dummy.sh
MORE="cat -"
which clear > $tempdir/install.$$ 2>&1
if test "$?" -eq 0 ; then CLEAR=clear; fi
which more > $tempdir/install.$$ 2>&1
if test "$?" -eq 0 ; then MORE=more; fi

rm -f $tempdir/install.$$

date=`date`
export XCRYSDEN_TOPDIR
. $XCRYSDEN_TOPDIR/scripts/xcConfigure_updateProfile.sh
. $XCRYSDEN_TOPDIR/scripts/xcConfigure_definitions.sh

rm -f $tempdir/xcInstall-sh.$$ $tempdir/xcInstall-csh.$$;
