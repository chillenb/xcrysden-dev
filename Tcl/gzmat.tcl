#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/gzmat.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ------------------------------------------------------------------------
# loads structure from Gaussian input file (requires babel)
# ------------------------------------------------------------------------
proc gzmat {filedir {can .mesa} {update 0}} {
    global xcMisc env system

    if { ! [info exists xcMisc(babel)] } {
	if [winfo exists .title] { 
	    destroy .title 
	}
	ErrorDialog "--gzmat option requires the definition of xcMisc(babel) in the ~/.xcrysden/custom-definitions"
	exit
    } else {
	set env(BABEL) $xcMisc(babel)
    }

    # OpenBabel seems not to need BABEL_DIR variable !!!

    if { [info exists xcMisc(babel_dir)] } {
	set env(BABEL_DIR) $xcMisc(babel_dir)
    }

    set head [file rootname [file tail $filedir]]

    #xcCatchExecReturnRedirectStdErr
    set cmd "exec sh $system(TOPDIR)/scripts/gzmat2xsf.sh $filedir > $system(SCRDIR)/$head.xsf 2> /dev/null"
    xcDebug -stderr "Executing: $cmd"
    if { [catch $cmd] } {
	destroyWelcome
	#ErrorDialog "an error occured while executing BABEL program."
        ErrorDialog "an error occured while parsing Gaussian Z-matrix input file"
        if { ![xcIsActive render] } {
            CloseCase
        }            
	return
    }

    xsfOpenParsed $system(SCRDIR)/$head.xsf $can $update
}
	

proc gzmat_menu {openwhat} {
    global system
    set filedir [tk_getOpenFile -defaultextension .xsf \
		     -filetypes { 
			 { {All Files}  {.*}  }
		     } -initialdir $system(PWD) -title $openwhat]
    if { $filedir == "" } {
	return
    }
    reloadRegister $filedir GZMAT
    gzmat $filedir
}
