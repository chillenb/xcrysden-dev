#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnReadStruct.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnReadStruct filehead {
    global system wn geng

    cd $wn(dir)
    if [catch {exec $system(FORDIR)/str2xcr $filehead}] {
	tk_dialog [WidgetName] ERROR \
		"ERROR while executing \"str2xcr\" program" \
		error 0 OK
    }
    # now copy $filehead.xcr to $system(SCRDIR)/xc_str2xcr.$$
    exec mv ${filehead}.xcr $system(SCRDIR)/xc_str2xcr.$system(PID)
    set dirname [file tail $filehead]
    exec mkdir $system(SCRDIR)/$dirname
    exec cp ${filehead}.struct $system(SCRDIR)/$dirname

    cd $system(SCRDIR)

    #
    # WIEN97 struct file is in BOHRs, thatwhy xc_str2xcr.$$ is in BOHRs
    #
    set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
    xcAppendState wien
    set periodic(igroup) 1
    xcDebug "Going to GenGeom"
    GenGeom $geng(M1_PRIM) $geng(M2_CELL) $geng(M3_ARGUMENT) \
	    $periodic(igroup) \
	    1 1 1 $system(SCRDIR)/xc_struc.$system(PID)

    xcDebug "Going to xc_readXSF"
    ################################################
    xc_readXSF $system(SCRDIR)/xc_struc.$system(PID)
    ################################################
}
