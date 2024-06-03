#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xsf2Open.tcl
# ------                                                                    #
# Copyright (c) 2017 by Anton Kokalj                                        #
#############################################################################

proc xsf2Open {filedir {can .mesa} {update 0}} {
    global xsf2 system geng xcMisc periodic radio sInfo xsfAnim working_XSF_file

    if { ! [info exists xcMisc(xsf2_manipulator)] } {
	ErrorDialog "Cannot open XSF.v2 files: xsf2_manipulator does not exist"
	CloseCase
	return
    }

    set filedir [gunzipXSF $filedir]
    UpdateWMTitle $filedir

    # so far XSF2 does not support ANIMSTEPS
    set fID  [open $filedir]
    set l1st [gets $fID]
    close $fID
    if { [string match "*ANIMSTEP*" $l1st] } {
	ErrorDialog "An Error occured, while reading XSF.v2 file $filedir.\n\nXSF.v2 currently does not support ANIMSTEPS keyword"
        if { ! $update } {
            CloseCase
        }
	return
    }
    
    # convert XSF2 to XSF
    if { [xsf2_ $filedir filexsf] } {
        if { ! $update } {
            CloseCase
        }
	return
    }
    
    if { ! [info exists xcMisc(reduce_to)] } {
	set xcMisc(reduce_to) {}
    }
    eval {exec $system(BINDIR)/xsf2xsf $filexsf $filexsf.raw} $xcMisc(reduce_to)
    set filexsf $filexsf.raw

    #
    # the working XSF file
    set working_XSF_file $filexsf

    if { $update } {
	GetDimGroupXSF periodic(dim) periodic(igroup) $filexsf
	puts stderr "*** update XSF:  dim = $periodic(dim)"
	if { $periodic(dim) < 3 && $radio(cellmode) == "conv" } {
	    set radio(cellmode) prim
	}
	if { $periodic(dim) > 1 } {
	    CellMode 1
	} else {
	    UpdateStruct .mesa $filexsf
	}
	xcUpdateState
    } else {
	#
	# now open the file
	#
	ResetDispModes
	if { [catch {xc_openstr xcr $filexsf $can PL}] } {
	    ErrorDialog "An Error occured, while reading XSF file $filexsf"
	    CloseCase
	    return
	}
	
	#
	# set the proper state
	#
	Get_sInfoArray
	DisplayDefaultMode
	xcAppendState render
	xcUpdateState
	set periodic(dim) $sInfo(dim)
	
	#
	# if the structure is periodic update it according to the "state"
	#
	if { $periodic(dim) > 0 } {
	    set geng(M3_ARGUMENT) [GetGengM3Arg ANGS]
	    
	    # the state should be set as well
	    cd $system(SCRDIR)
	    GenGeom $geng(M1_INFO) $geng(M2_CELL) $geng(M3_ARGUMENT) \
		1 1 1 1 xc_gengeom.$system(PID)

	    GetDimGroupXSF periodic(dim) periodic(igroup) $system(SCRDIR)/xc_gengeom.$system(PID)	    
	    CellMode 1
	    
	    Get_sInfoArray
	    xcUpdateState
	}
	
	CrysFrames 
    }
    return        
}

proc xsf2_ {xsf2 xsfVar} {
    global xcMisc
    upvar $xsfVar xsf

    if { ! [info exists xcMisc(xsf2_manipulator)] } {
	error "variable xcMisc(xsf2_manipulator) does not exist"
    }

    set xsf [file rootname $xsf2].xsf
    if { [string match $xsf2 $xsf] } {
	set xsf $xsf2.xsf
    }
    
    return [xcCatchExec $xcMisc(xsf2_manipulator) $xsf2 $xsf]
}
