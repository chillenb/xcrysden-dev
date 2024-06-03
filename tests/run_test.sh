#!/bin/sh

if test -z $XCRYSDEN_TOPDIR; then
    if test -f tests.sh; then
	export XCRYSDEN_TOPDIR=`(cd ../; pwd)`
	PATH=$XCRYSDEN_TOPDIR:$PATH	
    else
	echo "Please define the XCRYSDEN_TOPDIR variable or execute the script as ./make_all_tests.sh"
	exit 1
    fi
fi

. $XCRYSDEN_TOPDIR/tests/tests.sh

if test $# -ne 1 ; then
    echo "
Usage: $0 structures | wien | crystal | pwscf | scripting
"
    exit 1
fi

echo "
 ========================================================================
  
 Running test: $1 ...

 ========================================================================
" > $message_file

xterm -e tail -f $message_file &
#pid=$!

# execute the test
$1

echo "
Tests completed. Please press Ctrl-C.
" >> $message_file

# kill the xterm
#kill $pid
