#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/unmapWin.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# If one wants to hide toplevel with "Hide", then Toplevel must be 
# registered with the following command
# 
proc xcRegisterUnmapWindow {w bframe what {args {}}} {
    global unmapWin
    # w ........ unmap toplevel to register
    # bframe ... button's frame, i.e. which frame will hold the button
    # what ..... some kind of identifier
    set image  {}
    set bitmap {}
    set text   {}
    set textimage {}
    set i 0    
    foreach option $args {
	incr i
	if { $i%2 } {
	    set tag $option
	} else {
	    switch -- $tag {
		"-image"  {set image  $option}
		"-bitmap" {set bitmap $option}
		"-text"   {set text   $option}
		"-textimage" {set textimage $option}
		default { tk_dialog .mb_error Error \
			"ERROR: Bad xRegisterUnmapWindow's button configure option $tag" \
			error 0 OK
		return }
	    }
	}
    }

    if { $i%2 } {
	tk_dialog .mb_error1 Error "ERROR: You called xcMenuEntry with an odd number of args !" \
		error 0 OK
	return 0
    }

    lappend unmapWin($bframe,registered) $w

    if ![info exists unmapWin($bframe,column)] {
	set unmapWin($bframe,column) 0
	set unmapWin($bframe,count)  0
    } else {
	incr unmapWin($bframe,column)
	incr unmapWin($bframe,count)
    }

    if { $textimage == {} } {
	set unmapWin(but,$w,$what) \
		[button $bframe.unmap$unmapWin($bframe,count) -bd 1 \
		-text $text -bitmap $bitmap -image $image \
		-highlightthickness 0 \
		-command [list xcUnmapWindow map $w $w $bframe $what]]
    } else {
	set text  [lindex $textimage 0]
	set image [lindex $textimage 1]
	set unmapWin(but,$w,$what) \
		[xcHideButton $bframe.unmap$unmapWin($bframe,count) \
		$image right -text $text \
		-command [list xcUnmapWindow map $w $w $bframe $what]]
    }

    set unmapWin($w,$what,column) $unmapWin($bframe,column)
    # this should be changed:
    set unmapWin($w,$what,row)    0
}

#
# the "Hide" toplevel window must be have special <Map> <Unmap> bindings;
# here is the procedure for bindings
#
proc xcUnmapWindow {event type w bframe what} {
    global unmapWin
    # event ..... what event (map/unmap)
    # type, w ... to check if event was triggered by toplevel window
    # bframe .... button's frame, i.e. which frame will hold the button
    # what   .... some kind of identifier

    if { $type != $w } {
	# event wasn't triggered by toplevel window
	return
    }
    #xcDebug "Event:: $event $type"

    #
    # check if $w was registered
    #
    if ![info exists unmapWin($bframe,registered)] {
	return
    }
    set registered 0
    foreach win $unmapWin($bframe,registered) {
	if { $win == $w } {
	    set registered 1
	    break
	}
    }
    if !$registered {
	return
    }

    if { $event == "unmap" } {
	#xcDebug "Action: unmap"
	if ![winfo ismapped $bframe] {
	    eval pack $bframe $unmapWin(packinfo,$bframe)
	} elseif ![info exists unmapWin(packinfo,$bframe)] {
	    set unmapWin(packinfo,$bframe) [pack info $bframe]
	}
	grid configure $unmapWin(but,$w,$what) \
		-column $unmapWin($w,$what,column) -row $unmapWin($w,$what,row)
    } elseif { $event == "map" } {
	#xcDebug "Action: map"
	grid forget $unmapWin(but,$w,$what)
	xcDebug "HowManyUnmapWin:::: [HowManyUnmapWin $bframe]; bframe=$bframe"
	if { [HowManyUnmapWin $bframe] == 0 } {
	    pack forget $bframe 
	}
	wm deiconify $w
    }
}

#
# this is command that the "Hide" button calls
#
proc HideWin {t what} {
    global unmapWin
    wm withdraw $t
    xcUnmapWindow unmap $t $t $unmapWin(frame,main) $what
}


#------------------------------------------------------------------------------
# from here on are procedures needed for Hide technology
#
proc HowManyUnmapWin bframe {
    global unmapWin
    set count 0
    foreach elem [array names unmapWin but,*] {
	if [winfo ismapped $unmapWin($elem)] { incr count }
    }
    return $count
}


proc UnmapCleanAll {} {
    global unmapWin

    foreach elem [array names unmapWin but,*] {
	xcDebug "unmapWin Button:::       $elem"
	if [winfo exists $unmapWin($elem)] {
	    destroy $unmapWin($elem)
	}
    }
    # unset all unmapWin elements, but the bframes and packinfo,*
    foreach pattern {
	*,registered
	but,*
	*,column
	*,row
	*,count
    } {
	foreach elem [array names unmapWin $pattern] {
	    unset unmapWin($elem)
	}
    }
}
