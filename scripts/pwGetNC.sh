#!/bin/sh
#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwGetNC.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# set locales to C
LANG=C 
LC_ALL=C
export LANG LC_ALL

#
# Purpose: creates a "nuclei.charges" file
#
# Usage:   pwGetNC.sh  output|input pw-file
#

#######################################
if test "x`which readlink`" = "x"; then
    # no readlink cmd; make a function-substitute
    readlink_f() {
        file=`ls -l "$1" | awk '{print $NF}'`
        while test -h "$file"; do
            file=`ls -l "$file" | awk '{print $NF}'`
        done
        echo $file
    }
else
    readlink_f() {
        readlink -f "$1"
    }
fi


pathname() {
    file=`which "$1"`    
    if test $? -gt 0; then
        file=`type "$1"`
	if test $? -gt 0; then
	    # give-up
	    file="$1"
        else
            file=`echo $file | awk 'BEGIN {FS="is "} {print $NF}'`
        fi
    fi
    echo $file
}


pathdir() {
    file=`pathname "$1"`
    
    while test -h "$file"; do
	file=`readlink_f "$file"`
    done

    dir=`dirname "$file"`
    ( cd $dir; pwd )
}


if test -z $XCRYSDEN_TOPDIR; then
    # XCRYSDEN_TOPDIR does not exists, guess it from the process
    scriptdir=`pathdir $0`
    export XCRYSDEN_TOPDIR=`(cd $scriptdir/..; pwd)`
fi

if test -f $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh ; then
    . $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh
else
	echo "
ERROR: cannot create nuclei.charges file, because loading of pwLib_old.sh failed
"
	exit 1
fi
#######################################

if test $# -ne 2 ; then
    pwError "Usage:  pwGetNC.sh  output|input pw-file" 1
fi

input=1
if test "$1" = "output" ; then
    input=0
fi

#
# check if it is really PW-INPUT file
#
if test { "x`grep '&input' $2`" = "x" } -a { $input -eq 1 } ; then
    pwError "File \"$2\" is not a PW-INPUT file !!!" 1
fi

#
# check if it is really PW-OUTPUT file
#
if test { "x`grep 'Program PWSCF' $2`" = "x" } -a { $input -eq 0 } ; then
    pwError "File \"$2\" is not a PW-OUTPUT file !!!" 1
fi
  
pwNucleiCharges $2 /dev/null
exit 0
