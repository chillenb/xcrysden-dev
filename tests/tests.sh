#
# this file is sourced !!!
#
here=$XCRYSDEN_TOPDIR/tests
base=$XCRYSDEN_TOPDIR/examples
if test -z $XCRYSDEN_SCRATCH; then
    XCRYSDEN_SCRATCH=/tmp
fi
message_file=$XCRYSDEN_SCRATCH/xc_tests.$$

# ------------------------------------------------------------------------
# STRUCTURES: molecules & crystals examples
# ------------------------------------------------------------------------

structures() {
    echo "1. Executing MOLECULES examples ..." >> $message_file
    xcrysden -s $here/molecules.tcl

    echo "
2. Executing CRYSTALS examples ..." >> $message_file
    xcrysden -s $here/crystals.tcl
}


# ------------------------------------------------------------------------
# WIEN examples
# ------------------------------------------------------------------------

wien() {
    if test -d $HOME/test/WIEN97/FermiSurface/tic ; then 
	echo "
3. Executing WIEN examples ..." >> $message_file
	echo "   3.1 Executing Renderdensity-WIEN examples ..." >> $message_file
	xcrysden --wien_renderdensity $HOME/test/WIEN97/FermiSurface/tic
	echo "   3.2 Executing Density-WIEN examples ..." >> $message_file
	xcrysden --wien_density       $HOME/test/WIEN97/FermiSurface/tic
	echo "   3.3 Executing Kpath-WIEN examples ..." >> $message_file
	xcrysden --wien_kpath         $HOME/test/WIEN97/FermiSurface/tic
	echo "   3.4 Executing Fermisurface-WIEN examples ..." >> $message_file
	xcrysden --wien_fermis        $HOME/test/WIEN97/FermiSurface/tic
    fi
}


# ------------------------------------------------------------------------
# CRYSTAL inp|f9 examples
# ------------------------------------------------------------------------

crystal() {
    echo "
4. Executing CRYSTALxx examples ..." >> $message_file
    echo "   4.1 Please test CRYSTALxx AdvGeom options ..." >> $message_file
    xcrysden --crystal_inp $base/CRYSTALxx_input_files/Pt_fcc.r1
    if test -f $HOME/test/crystal14/urea/urea.f9 ; then 
	echo "   4.2 Please test CRYSTALxx Properties options ..." >> $message_file
	xcrysden --crystal_f9  $HOME/test/crystal14/urea/urea.f9
    fi
    if test -f $HOME/test/crystal14/urea_UHF/urea_UHF.f9 ; then
	echo "   4.3 Please test CRYSTALxx Properties options (SPIN-POLARISED) ..." >> $message_file
	xcrysden --crystal_f9  $HOME/test/crystal14/urea_UHF/urea_UHF.f9
    fi
}


# ------------------------------------------------------------------------
# PWscf I/O examples
# ------------------------------------------------------------------------
pwscf() {
    echo "
5. Executing PWSCF examples ..." >> $message_file
    xcrysden -s $here/pwscf.tcl
}


# ------------------------------------------------------------------------
# Scripting examples
# ------------------------------------------------------------------------
scripting() {
    cd $base/Scripting

    echo "
5. Executing Scripting examples ..." >> $message_file
    
    for file in *.tcl
      do
      echo "   5.x Executing Scripting example: $file ..." >> $message_file
      
      xcrysden -s $file
    done
}
