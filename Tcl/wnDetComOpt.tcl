#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnDetComOpt.tcl                                  #
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnDetComOpt dir {
    global wn

    set wn(complex)        0
    set wn(spin_polarized) 0

    set fh [file tail $dir]
    if { [file exists $dir/$fh.in1c] } { 
	if {[file size $dir/$fh.in1c] > 0 } { 
	    set wn(complex) 1 
	}
    }
    if { [file exists $dir/$fh.clmdn] } {
	if {[file size $dir/$fh.clmdn] > 0 } { 
	    set wn(spin_polarized) 1
	}
    }
}
