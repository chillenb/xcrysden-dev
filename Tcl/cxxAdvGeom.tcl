#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/cxxAdvGeom.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc cxxAdvGeom.manualOption {{retry 0}} {
    global system AdvGeom

    set file $system(SCRDIR)/CRYSTAL.option

    set helpText {#
# Please add an "advanced geometrical" option manually. Refer
# to CRYSTAL Manual if you want to get a help on 
# "advanced geometrical option"
#
# Start at the begining of next line:
    }

    if { $retry == 0 } {
	WriteFile $file $helpText w
    } else {
	set old  [ReadFile -nonewline $file]
	set text [format "%s%s\n" $helpText $old]
	WriteFile $file $text w
    }
    update
    xcEditFile $file -foreground

    #
    # drop the comments from file
    #
    set content {}
    foreach line [split [ReadFile -nonewline $file] \n] {
	if { [string match *\#* $line] == 0 } {
	    if { $line != {} } {
		append option [format "%s\n" $line]
	    }
	}
    }

    xcDebug -stderr "manualOption:"
    xcDebug -stderr "-------------"
    xcDebug -stderr $option

    #
    # register the option
    #
    set n [xcAdvGeomState new]
    set AdvGeom($n,option) $option
    
    #
    # check if option is OK
    #
    set cxxInput [MakeInput]
    set inp xc_inp.$system(PID)
    set out xc_out.$system(PID)
    WriteFile $inp $cxxInput w

    set status [cxxAdvGeom.testINPUT $inp $out]
    
    if { $status == "noretry" } {
	xcAdvGeomState delete
	return 0
    } elseif { $status == "retry" } {
	xcAdvGeomState delete
	set status [cxxAdvGeom.manualOption retry]
	return $status
    } else {
	#
	# option is OK -> update the display
	#
	GenCommUndoRedo "Add an Option Manually"
	CalStru
	xcUpdateState

	return 1
    }
}



proc cxxAdvGeom.viewScript {} {
    xcDisplayVarText [MakeInput] {CRYSTAL Input Script}
}


proc cxxAdvGeom.manualEdit {{retry 0}} {
    global system cxx AdvGeom

    if { $retry == 0 } {
	set input [MakeInput]
    } else {
	set input $cxx(tmpInputContent)
    }

    #
    # edit a file
    # 
    cd $system(SCRDIR)
    set inp   xc_inp.$system(PID)
    set out   xc_out.$system(PID)
    WriteFile $inp $input w
    update
    xcEditFile $inp -foreground
    
    # handle correctly the EXTPRT/COORPRT/STOP keywords
    set cxx(tmpInputContent) [cxxHandleEXTPRT [ReadFile -nonewline $inp]]
    WriteFile $inp $cxx(tmpInputContent) w

    #
    # test a new file
    #
    set status [cxxAdvGeom.testINPUT $inp $out]
    
    if { $status == "noretry" } {
	return 0
    } elseif { $status == "retry" } {
	set status [cxxAdvGeom.manualEdit retry]
	return $status
    } else {
	#
	# option is OK -> proceed
	#
	
	# register the option
	set n [xcAdvGeomState new]
	set AdvGeom($n,edit) $cxx(tmpInputContent)
	
	GenCommUndoRedo "Edit Manually"
	CalStru
	xcUpdateState
	return 1
    }
}


proc  cxxAdvGeom.testINPUT {inp out} {
    global system
    
    set catchCode [catch {exec $system(c95_integrals) < $inp > $out} errMsg]
    if { $errMsg == {} } {
	set errMsg "CRYSTAL module: $system(c95_integrals) exited with and 
exit status 0, but the \"ERROR ****\" string 
exists in the output"
    }
    if { $errMsg == "FORTRAN STOP" && $catchCode } {
	# printing the "FORTRAN STOP" to stderr caused the error.
	set catchCode 0
    }

    set content [ReadFile $out]
    if { [string match "*ERROR \*\*\*\**" $content] || $catchCode > 0 } {
	# option is BAD - and error occured!!!
	
	set id [tk_dialog [WidgetName] ERROR \
		    "An ERROR occur while executing CRYSTAL module: $system(c95_integrals)" error 0 OK ErrorInfo {View Crystal Ouput}]
	if { $id == 1 } {
	    set t [xcDisplayVarText $errMsg {Error Info}]
	    tkwait window $t
	} elseif { $id == 2 } {
	    set cxxOutput [ReadFile $out]
	    set t [xcDisplayVarText $cxxOutput {Crystal Ouput}]
	    tkwait window $t
	}
	
	# ask user "Do you want top retry?
	
	set id [tk_dialog [WidgetName] QUESTION \
		    "Do You want to retry ?" question 1 No Yes]
	if { $id == 1 } {
	    return retry
	} else {
	    return noretry
	}
    } else {
	return 1
    }
}


proc cxxHandleEXTPRT {geoInput} {
    # maybe EXTPRT/COORPRT/STOP are already specified, but there may
    # be some additional geometry manipulation after EXTPRT keyword,
    # so the only safe thing is to throw that out and specifying it
    # again
    xcDebug -debug "bug-fixing: geoInput: $geoInput"
    set swap $geoInput
    regsub -all EXTPRT|COORPRT|STOP $swap {} geoInput
    append geoInput "\nEXTPRT\nCOORPRT\nSTOP\n"

    # skip empty lines
    set geoInput [xcSkipEmptyLines $geoInput]
    xcDebug -debug "bug-fixing: geoInput_new: $geoInput"
    return $geoInput
}