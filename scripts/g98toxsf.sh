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
# Source: $XCRYSDEN_TOPDIR/scripts/g98toxsf.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
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
Usage:    g98toxsf.sh G98-output > XSF-file
      or
          g98toxsf.sh < G98-output > XSF-file
"
    exit 1
fi


#--------------------------------------------------------------------------
# this is an experimental Gaussian98 To XSF converter. Use at your own risk
#
# Usage:   g98toxsf.sh G98-output > XSF-file
#        or
#          g98toxsf.sh < G98-output > XSF-file
#
# The order of coordinates searching is the following:
#
# ----
# 1.) first check if the following is present
#
# 0 Center     Atomic                   Forces (Hartrees/Bohr)
# 1 Number     Number              X              Y              Z
# 2 -------------------------------------------------------------------
# 3    1          8           -.051649003     .000000000    -.039990931
# 4    2          1            .017695138     .000000000     .030494684
# 5    3          1            .033953865     .000000000     .009496247
# 6 -------------------------------------------------------------------
#
#     if YES, then look for: 
# 
# 0                          Input orientation:                          
# 1 ---------------------------------------------------------------------
# 2 Center     Atomic     Atomic              Coordinates (Angstroms)
# 3 Number     Number      Type              X           Y           Z
# 4 ---------------------------------------------------------------------
# 5    1          8             0         .000000     .000000     .000000
# 6    2          1             0         .000000     .000000    1.000000
# 7    3          1             0         .968148     .000000    -.250380
# 8 ---------------------------------------------------------------------
#
# ----
# 2.) if forces are not present, then look for:
#
#                         Standard orientation:                         
# ---------------------------------------------------------------------
# Center     Atomic     Atomic              Coordinates (Angstroms)
# Number     Number      Type              X           Y           Z
# ---------------------------------------------------------------------
#    1         46             0         .000000     .000000    1.077941
#    2          8             0         .000000     .000000   -1.088059
#    3          7             0         .000000     .000000   -2.335059
#    4          7             0         .000000     .000000   -3.505059
# ---------------------------------------------------------------------
#
# ----
# 3.) if Standard orientation is not present, then look for:
#
#                          Z-Matrix orientation:                         
# ---------------------------------------------------------------------
# Center     Atomic     Atomic              Coordinates (Angstroms)
# Number     Number      Type              X           Y           Z
# ---------------------------------------------------------------------
#    1         46             0         .000000     .000000     .000000
#    2          8             0         .000000     .000000    2.166000
#    3          7             0         .000000     .000000    3.413000
#    4          7             0         .000000     .000000    4.583000
# ---------------------------------------------------------------------
#
#------------------------------------------------------------------------


forces_exists=`egrep 'Forces \(Hartrees\/Bohr\)' $input`
input_exists=`grep ' Input orientation:' $input`
standard_exists=`grep ' Standard orientation:' $input`
Zmatrix_exist=`grep ' Z-Matrix orientation:' $input`

if [ "x$forces_exists" != "x"  -a  "x$input_exists" != "x" ]; then
    #
    # search for:
    # 
    # 0 Center     Atomic                   Forces (Hartrees/Bohr)
    # 1 Number     Number              X              Y              Z
    # 2 -------------------------------------------------------------------
    # 3    1          8           -.051649003     .000000000    -.039990931
    # 4    2          1            .017695138     .000000000     .030494684
    # 5    3          1            .033953865     .000000000     .009496247
    # 6 -------------------------------------------------------------------
    #
    #     if YES, then look for: 
    # 
    # 0                          Input orientation:                          
    # 1 ---------------------------------------------------------------------
    # 2 Center     Atomic     Atomic              Coordinates (Angstroms)
    # 3 Number     Number      Type              X           Y           Z
    # 4 ---------------------------------------------------------------------
    # 5    1          8             0         .000000     .000000     .000000
    # 6    2          1             0         .000000     .000000    1.000000
    # 7    3          1             0         .968148     .000000    -.250380
    # 8 ---------------------------------------------------------------------
    #
    
    cat $input | awk '
{
  if ( $0 ~ / +Input orientation:|Z-Matrix orientation:/ ) {    
    i=0;
    getline;
    getline;
    getline;
    getline;
    getline;
    while ( $1 !~ /\-\-/ ) {
      #print "COOR:", $0;
      atn[i] = $2;
      x[i]   = $4;
      y[i]   = $5;
      z[i++] = $6;
      getline;
    }
    atoms_printed = 0;
  }
}
{
  if ( $0 ~ /Forces \(Hartrees\/Bohr\)/ ) {
    #print "LINE:", $0;
    i=0;
    getline;
    getline;
    getline;
    printf " %s\n", "ATOMS"; 
    while ( $1 !~ /\-\-/ ) {
      printf "%3d   %15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f\n", 
	atn[i], x[i], y[i], z[i], $3, $4, $5;
      i++;
      getline;
    }
    atoms_printed = 1;
  }
}
END {
  if ( !atoms_printed ) {
    # print OPTIMIZED or the LAST coordinates
    printf " %s\n", "ATOMS"; 
    for (ia=0; ia<i; ia++) 
      printf "%3d   %15.10f %15.10f %15.10f   %15.10f %15.10f %15.10f\n", 
	atn[ia], x[ia], y[ia], z[ia], 0.0, 0.0, 0.0;
  }
}' > g98toxsf.$$


elif [ "$Zmatrix_exist" != "" ]; then
    #
    # search for:
    #
    #                          Z-Matrix orientation:                         
    # ---------------------------------------------------------------------
    # Center     Atomic     Atomic              Coordinates (Angstroms)
    # Number     Number      Type              X           Y           Z
    # ---------------------------------------------------------------------
    #    1         46             0         .000000     .000000     .000000
    #    2          8             0         .000000     .000000    2.166000
    #    3          7             0         .000000     .000000    3.413000
    #    4          7             0         .000000     .000000    4.583000
    # ---------------------------------------------------------------------

    cat $input | awk '
{
  if ( $0 ~ / +Z-Matrix orientation:/ ) {
    getline;
    getline;
    getline;
    getline;
    getline;
    printf " %s\n", "ATOMS"; 
    while ( $1 !~ /\-\-/ ) {
      printf "%3d   %15.10f %15.10f %15.10f\n", $2, $4, $5, $6;
      getline;
    }
  }
}' > g98toxsf.$$

elif [ "$standard_exists" != "" ]; then
    #
    # search for:
    #
    #                         Standard orientation:                         
    # ---------------------------------------------------------------------
    # Center     Atomic     Atomic              Coordinates (Angstroms)
    # Number     Number      Type              X           Y           Z
    # ---------------------------------------------------------------------
    #    1         46             0         .000000     .000000    1.077941
    #    2          8             0         .000000     .000000   -1.088059
    #    3          7             0         .000000     .000000   -2.335059
    #    4          7             0         .000000     .000000   -3.505059
    # ---------------------------------------------------------------------

    cat $input | awk '
{
  if ( $0 ~ / +Standard orientation:/ ) {
    getline;
    getline;
    getline;
    getline;
    getline;
    printf " %s\n", "ATOMS"; 
    while ( $1 !~ /\-\-/ ) {
      printf "%3d   %15.10f %15.10f %15.10f\n", $2, $4, $5, $6;
      getline;
    }
  }
}' > g98toxsf.$$

fi

if [ -f g98toxsf.$$ ]; then
    # this is common to all !!!
    nstep=`grep ATOMS g98toxsf.$$ | wc | awk '{print $1}'`
    if [ $nstep -gt 1 ]; then
# make a AXSF file
	cat g98toxsf.$$ | awk -v ns=$nstep '
BEGIN {
  printf " %s %d\n", "ANIMSTEPS", ns;
}
/ATOMS/ { printf " %s %d\n", "ATOMS", ++i; next; }
/a*/    { print; }'
    else
	cat g98toxsf.$$
    fi
    
    rm -f g98toxsf.$$
    exit 0
else
    # coordinates were not extracted ...
    exit 1
fi



    
