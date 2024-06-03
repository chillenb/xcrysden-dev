#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

ForLoop() {
    # Usege: $0 from to incr   
    from=$1
    to=$2
    if test $# -eq 2 ; then
	incr=1
    else
	incr=$3
    fi
    op=-le
    if test $incr -lt 0 ; then
	op=-ge
    fi

    i=$from
    while [ $i $op $to ]
    do
	printf "%d " $i
	i=`expr $i + $incr`
    done
}

pwError() {
    # Usage: $0 message status
    echo "
 ========================================================================
    $1
 ========================================================================
"
    if [ "$2" -ge 0 ]; then
	exit $2
    fi
}

# --------------------------------------------------------------------------
# pwNucleiCharges --
#
# Purpose:    ityp->nat conversion data 
#
# Usage:      pwNucleiCharges pw_input|pw_output outfile
#
# Side efect: creates nuclei.charges file
# --------------------------------------------------------------------------

pwNucleiCharges() {
    #
    # if file nuclei.charges does not exists prompt for ityp->nat conversion !!
    #

    if [ \( "$1" = "" \) -o \( "$2" = "" \) ]; then
	pwError "Usage:  pwNucleiCharges  pw_input|pw_output  outfile" 1
    fi
    
    # do we have PW-INPUT or PW-OUTPUT file ???
    
    if [ "`cat $1 | egrep -i '&input|&system'`" != "" ]; then
	# it is PW-INPUT
	ntyp=`cat "$1" | awk '{gsub(",","\n"); print}' | grep ntyp \
		| awk '{split($0,a,"=|,"); print a[2];}'`
    else
	# PW-OUTPUT
	ntyp=`cat "$1" | grep 'number of atomic types' | \
	    head -1 | awk '{print $NF}'`
	#echo 'NTYP=$ntyp'
	if [ "$ntyp" = "" ]; then
	    # some older PWSCF versions didn't have "number of atomic
	    # types" printout -> user will have to make nuclei.charges
	    # file by himself/herself !!!
	    pwError "This is either non PW-output file or is a PW-output file 
    produced with some old PWSCF version" -1
	    echo -n "How many ityp->nat replacements ? "
	    read ntyp
	fi
    fi		     

    if [ ! -f nuclei.charges ]; then
	echo -n "Please enter $ntyp ityp->nat replacements !!! "
	echo $ntyp > nuclei.charges
	i=0
	while [ $i -lt "$ntyp" ]
	do
	    i=`echo "$i + 1"|bc`
	    echo ""
	    echo "Replacement #${i}: ityp->nat"
	    echo -n "ityp[$i]=$i; nat[$i]="; read nat
	    echo "$i $nat" >> nuclei.charges
	done
    fi
	
    cat nuclei.charges > "$2"
}
