#!/bin/sh

if test -z $XCRYSDEN_TOPDIR; then
    if test -f tests.sh; then
	XCRYSDEN_TOPDIR=`(cd ../; pwd)`
    else
	echo "Please define the XCRYSDEN_TOPDIR variable or execute the script as ./make_all_tests.sh"
	exit 1
    fi
fi

. $XCRYSDEN_TOPDIR/tests/tests.sh

echo "Making all XCRYSDEN tests ...
" > $message_file

xterm -e tail -f $message_file &

# now run the examples
structures
wien
crystal
pwscf
scripting

echo "
Tests completed. Please press Ctrl-C.
" >> $message_file