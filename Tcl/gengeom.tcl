#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/gengeom.tcl                                      #
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc GetGengM3Arg {mode3 {mode4 {}}} {    
    global geng

    if { $mode3 == "ANGS" } {
	set m31 1
    } elseif { $mode3 == "BOHR" } {
	set m31 2
    }

    if { $mode4 == "95" } {
	set m32 1
    } elseif { $mode4 == "98" || $mode4 == "03" || $mode4 == "06" || $mode4 == "09" || $mode4 == "14" } {
	set m32 2
    } else {
	set m32 1
    }

    return ${m31}${m32}
}


proc GenGeom {m1 m2 m3 igroup nx ny nz output} {
    global system xcMisc periodic working_XSF_file

    if { [xcIsActive c95] || [xcIsActive external34] } { 
	# for CRYSTALXX or EXTERNAL34 read from unit 34
	xcDebug -debug "exec $system(BINDIR)/gengeom \
		$m1 $m2 $m3 $igroup $nx $ny $nz $output"	
	xcCatchExecReturn $system(BINDIR)/gengeom $m1 $m2 $m3 $igroup $nx $ny $nz $output
    } elseif [xcIsActive wien] {
	# for WIENXX read from xc_str2xcr.$$
	xcDebug -debug "exec $system(BINDIR)/gengeom \
		$m1 $m2 $m3 $igroup $nx $ny $nz $output \
		$system(SCRDIR)/xc_str2xcr.$system(PID)"
	xcCatchExecReturn $system(BINDIR)/gengeom $m1 $m2 $m3 $igroup \
		$nx $ny $nz $output $system(SCRDIR)/xc_str2xcr.$system(PID)
    } elseif [xcIsActive external] {
	# for EXTERNAL read from $xcMisc(external_xsf_name)
	xcDebug -debug "exec $system(BINDIR)/gengeom \
		$m1 $m2 $m3 $igroup $nx $ny $nz $output \
		$xcMisc(external_xsf_name)"	
	xcCatchExecReturn $system(BINDIR)/gengeom $m1 $m2 $m3 $igroup \
		$nx $ny $nz $output $xcMisc(external_xsf_name)
    } elseif { [info exists working_XSF_file] } {
	# it must be some XSF struct file
	xcDebug -debug "exec $system(BINDIR)/gengeom \
		$m1 $m2 $m3 $igroup $nx $ny $nz $output $working_XSF_file"
	xcCatchExecReturn $system(BINDIR)/gengeom $m1 $m2 $m3 $igroup \
		$nx $ny $nz $output $working_XSF_file
    } else {
	ErrorDialog "don't know what to do in GenGeom (this is a bug)"
	return
    }

    #
    # query dimensionality of the sistem
    #
    set channel [open $output r]
    set dim 0
    while { ! [eof $channel] } {
	set line [gets $channel]
	if { [string match *DIM* $line] } {
	    set dim [lindex [gets $channel] 0]
	    break
	}
    }
    close $channel
    set periodic(dim) $dim
}

