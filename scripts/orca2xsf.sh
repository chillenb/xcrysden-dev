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
# Source: $XCRYSDEN_TOPDIR/scripts/orca2xsf.sh
# ------                                                                    #
# Copyright (c) 2016 by Anton Kokalj                                        #
#############################################################################

# set locales to C
LANG=C 
LC_ALL=C
export LANG LC_ALL

if [ $# -eq 0 ]; then
    input=-
elif [ $# -eq 1 ]; then
    input=$1
else
    echo "
Usage:    orca2xsf.sh orca-output > XSF-file
      or
          orca2xsf.sh < orca-output > XSF-file
"
    exit 1
fi


#------------------------------------------------------------------------
# This is an experimental Orca To XSF converter. Use at your own risk
#------------------------------------------------------------------------
#
# The script greps the Cartesian that follows the following record:
#
# ^---------------------------------
# ^CARTESIAN COORDINATES (ANGSTROEM)
# ^---------------------------------
#
# where "^" marks the beggining of line.
#
# This coordinates are augmented with atomic forces if the following record is found:
#
# ^------------------
# ^CARTESIAN GRADIENT
# ^------------------
#
# ------------------------------------------------------------------------
forc_str='^CARTESIAN GRADIENT'
coor_str='^CARTESIAN COORDINATES \(ANGSTRO'
forces_exists=`egrep -a "$forc_str" $input`
coord_exists=`egrep -a "$coor_str" $input`

if [ "x$forces_exists" != "x"  -a  "x$coord_exists" != "x" ]; then
    
    ncoord=`egrep -a "$coor_str" $input | nl | tail -1 | awk '{print $1}'`

    cat $input | awk -v ncoord=$ncoord '
BEGIN { 
   nimages=0; 
   if (ncoord>1) print "ANIMSTEPS", ncoord;
}
{
  if ( $0 ~ /^CARTESIAN COORDINATES \(ANGSTRO/ ) {    
    i=0;
    getline;
    while ( NF != 4 ) { getline; }
    while ( NF == 4 ) {
      atn[i] = $1;
      x[i]   = $2;
      y[i]   = $3;
      z[i++] = $4;
      getline;
    }
    natoms = i;
    atoms_printed = 0;
  }
}
{
  if ( $0 ~ /^CARTESIAN GRADIENT/ ) {
    i=0;
    getline;
    while ( NF != 6) { getline; }
    while ( NF == 6) {
      fx[i]   = $4;
      fy[i]   = $5;
      fz[i++] = $6;
      getline;
     }
     nforc = i;

     if ( natoms != nforc ) {
        print "### ERROR parsing output: number of atoms and atomic forces does not match"
        exit 1;
     } else { 
        if (ncoord>1) { print "ATOMS", ++nimages; }
        else { print "ATOMS"; }

        for (i=0; i<natoms; i++) 
          {
            printf "%3s   %15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f\n", atn[i], x[i], y[i], z[i], fx[i], fy[i], fz[i];
          }
        atoms_printed = 1;
     }
  }
}
END {
  if (!atoms_printed) {

    # it appears the last coordinates are w/o forces

    if (ncoord>1) { print "ATOMS", ++nimages; }
    else { print "ATOMS"; }

    for (i=0; i<natoms; i++) 
        printf "%3s   %15.10f %15.10f %15.10f\n", atn[i], x[i], y[i], z[i];
  }
}'


elif [ "x$coord_exists" != "x" ]; then

    # is it possible to have more than one coordinates in the output w/o forces?
    # Anyway just in case ...

    ncoord=`egrep -a "$coor_str" $input | nl | tail -1 | awk '{print $1}'`

        cat $input | awk -v ncoord=$ncoord '
BEGIN { 
   nimages=0; 
   if (ncoord>1) print "ANIMSTEPS", ncoord;
}
{
  if ( $0 ~ /^CARTESIAN COORDINATES \(ANGSTRO/ ) {    
    i=0;
    getline;
    while ( NF != 4 ) { getline; }

    if (ncoord>1) { print "ATOMS", ++nimages; }
    else { print "ATOMS"; }

    while ( NF == 4 ) {
      printf "%3s   %15.10f %15.10f %15.10f\n", $1, $2, $3, $4;
      getline;
    }
  }
}'
	
fi
