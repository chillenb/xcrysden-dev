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
# Source: $XCRYSDEN_TOPDIR/scripts/pwi2xsf.sh
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# set locales to C
LANG=C 
LC_ALL=C
export LANG LC_ALL

#
# pwi2xsf.sh: PW-input to XSF converison
#
# Usage: pwi2xsf [-r] pw-input-file
#
# Written by Tone Kokalj on Tue May  8 20:43:44 2001
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

if [ "$#" -lt 1 ]; then
    echo "
Usage: pwi2xsf.sh [-r] pw-input

Option for PWscf version < 1.2:

-r ... for pw.x < 1.2 one must spefify the ityp->nat conversion, and
       the corresponding data are writen to file nuclei.charges. The
       -r flag deletes this file as to force new ityp->nat specification.
"
    exit 1
fi

r=0
if [ "$1" = "-r" ]; then
    r=1
    shift
fi


#######################################
if test "x`which tee`" = "x"; then
    # no tee cmd; make a function-substitute
    tee() {
	cat - > $1
	cat $1
    }
fi

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
    script_dir=`pathdir $0`
    export XCRYSDEN_TOPDIR=`(cd $script_dir/..; pwd)`
fi

if test -f $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh ; then
    . $XCRYSDEN_TOPDIR/scripts/pwLib_old.sh
    load_old_lib=1
else
    load_old_lib=0
fi

#######################################


# ------------------------------------------------------------------------
# Function: pwi2xsf_via_pwo
#
# Purpose: for "ATOMIC_POSITIONS crystal_sg", get the coordinates from
# pw.x output by making a dry pw.x run
# ------------------------------------------------------------------------
pwi2xsf_via_pwo() {
    PWO_XSF2XSF=pwo_xsf2xsf
    if test -f $XCRYSDEN_TOPDIR/bin/pwo_xsf2xsf ; then
        PWO_XSF2XSF=$XCRYSDEN_TOPDIR/bin/pwo_xsf2xsf
    elif test -f $XCRYSDEN_LIB_BINDIR/pwo_xsf2xsf ; then
        PWO_XSF2XSF=$XCRYSDEN_LIB_BINDIR/pwo_xsf2xsf
    fi

    cat $1 | awk 'BEGIN {RS="[,\n]"; elec_nml=0; elec_print=0; }
toupper($0) ~ /&ELECTRONS/ { elec_nml = 1; }

toupper($0) ~ /&END|^\/|^ +\// { 
  if (elec_nml && !elec_print) { 
     print "   electron_maxstep = 0";
     elec_nml = 0;
  }
}

/=/ { 
  if ( toupper($1) ~ /PREFIX/ ) { 
     # skip prefix, because we will make pwscf.EXIT stop file
     next;
  }
  if ( toupper($1) ~ /ELECTRON_MAXSTEP/ ) { 
     print "   electron_maxstep = 0"; elec_print = 1; next;
  }
}

/a*/ {
   print $0;
}' > pwi.$$

    # do we have pw.x ?
    
    file=`which pw.x`
    if test $? -gt 0; then
	file=`type pw.x`
	if test $? -gt 0; then
	    # give-up
            echo "pw.x does not exist; cannot parse pw.x input with \"ATOMIC_POSITIONS crystal_sg\""
            rm -f pwi.$$
            exit 1
	fi
    fi

    touch pwscf.EXIT
    pw.x < pwi.$$ > pwo.$$ 2> /dev/null
    $XCRYSDEN_TOPDIR/scripts/pwo2xsf.sh -ic pwo.$$ | tee pwi2xsf.xsf_out
    rm -f pwi.$$ pwo.$$
    exit 0
}


#
# check if we have OLD or NEW PW.X input format
#
new_format1=`grep 'ATOMIC_POSITIONS' $1`
new_format2=`grep -i '&system' $1`

if [ "$new_format1" != ""  -a  "$new_format2" != "" ]; then
    #
    # we have NEW PW.X input format
    #
    # check if ATOMIC_POSITIONS are specified in crystal_sg units
    sg=`cat $1 | grep ATOMIC_POSITIONS | grep -i crystal_sg`

    if test "x$sg" != "x"; then
        # coordinates are specified in crystal_sg units
        pwi2xsf_via_pwo $1       
    else
        #
        cat $1 | awk 'BEGIN {RS="[,\n]";} {print $0}' | awk '
BEGIN {
  calculation="";
  num_of_images="";
  nml_end=0;
  nml_end_string="";
  new_neb=0;      
}

toupper($0) ~ /BEGIN_ENGINE_INPUT/ { new_neb=1; }

$1 ~ /^BEGIN$|^END$|^BEGIN_PATH_INPUT$|^END_PATH_INPUT$|^BEGIN_ENGINE_INPUT$|^END_ENGINE_INPUT$|^BEGIN_POSITIONS$|^END_POSITIONS$/ { next; }

toupper($0) ~ /&SYSTEM/          { print; }

/=/ { 
  if ( toupper($1) ~ /^IBRAV($|=)|^CELLDM\([1-6]\)($|=)|^NAT($|=)|^A($|=)|^B($|=)|^C($|=)|^COSAB($|=)|^COSAC($|=)|^COSBC($|=)/ ) { print; } 
  
  if ( toupper($1) ~ /^CALCULATION($|=)/ ) { calculation=toupper($0); }

  if ( toupper($1) ~ /^NUM_OF_IMAGES($|=)/ ) { num_of_images=toupper($0); }
}

$1 ~ /^INTERMEDIATE_IMAGE$|^LAST_IMAGE$/ { print; next; }

$1 ~ /^FIRST_IMAGE$|^ATOMIC_POSITIONS$|^CELL_PARAMETERS$/ {
  if ( !nml_end) {
     # first finish the &SYSTEM namelist
     nml_end=1;
     if (new_neb) { 
        print "   calculation = \"PATH\"";
     } else if (calculation != "") {
        print calculation;
     }       
     if (num_of_images != "") print num_of_images;
     print nml_end_string;
  }
  # now print the current record
  print_line=1;
  print toupper($0);   
  next;
}

# old neb stuff:
NF == 1 && /first_image|intermediate_image|last_image/ { print toupper($1); next; }

toupper($0) ~ /&END|^\/|^ +\// { 
  nml_end_string=$0;
}

/a*/ {
  if ( print_line == 1 ) {
    print $0;
  }
}'> pw.$$

        PWI2XSF=pwi2xsf
    fi
else
    #
    # we have OLD PW.X input format	
    #
    
    if test $load_old_lib = 0; then
	echo "
ERROR: cannot convert to XSF, because loading of pwLib_old.sh failed
"
	exit 1
    fi

    pwNucleiCharges $1 /dev/null

    cat $1 | awk 'BEGIN {RS="[,\n]";} {print}' | awk '
BEGIN {
    end=0;
}
toupper($0) ~ /&INPUT|CELLDM|NAT|LTAUCRY/ { print; }
toupper($0) ~ /IBRAV/ { 
    print;
    split($0,a,"=");
    split(a[1],b,",");
    ibrav = b[1];
}
toupper($0) ~ /&END|^\/|^ \// { end=1; }
/a*/ {
    if ( end == 1 ) print;
}' > pw.$$
    PWI2XSF=pwi2xsf_old
fi

#
# execute $PWI2XSF fortran program and print the XSF file
#
if test -f $XCRYSDEN_TOPDIR/bin/$PWI2XSF ; then
    $XCRYSDEN_TOPDIR/bin/$PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
elif test -f $XCRYSDEN_LIB_BINDIR/$PWI2XSF ; then
    $XCRYSDEN_LIB_BINDIR/$PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
else
    $PWI2XSF < pw.$$ | tee pwi2xsf.xsf_out
fi
#rm -f pw.$$

if [ "$r" -eq 1 ]; then
    if [ -f nuclei.charges ]; then rm nuclei.charges; fi
fi

for file in pw.$$
do
    if test -f $file; then rm -f $file; fi
done

exit 0
