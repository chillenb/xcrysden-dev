#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/xcLib.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

xcSignalHandler () {
    # usage: xcSignalHandler program signal ID
    if [ "$1" = "EnScan" ]; then
	rm -rf EnScan1.$$ EnScan2.$$ EnScan3.$$ EnScan_*.$!
	if [ -f EnSCan22.$$ ]; then
	    rm -f EnScan22.$$
	fi
    fi
    case $2 in
	1)  signal='Hangup (1)';;
	2)  signal='Interupt (2)';;
	3)  signal='Quit (3)';;
	15) signal='Software termination (kill) (15)';;
	*)  signal='Unknown signal (x)';;
    esac
    echo "$1 has received the following signal: $signal"
    kill -15 $!
    if [ -d $XCRYSDEN_SCRATCH/xc_$! ]; then
	echo "Deleting directory: $XCRYSDEN_SCRATCH/xc_$! ..."
	rm -rf $XCRYSDEN_SCRATCH/xc_$!
    fi	    
    if [ -f $XCRYSDEN_SCRATCH/STDIN.$$ ]; then
	echo "Deleting file:      $XCRYSDEN_SCRATCH/STDIN.$$ ..."	
	rm -f $XCRYSDEN_SCRATCH/STDIN.$$
    fi
    if [ -f $PWD/core ]; then
	rm -f core
    fi
    echo "Quit !"
    exit $2
}
