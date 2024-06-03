#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/dispModes.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc Lighting {var} {
    global dispmode xcFonts xcMisc light
    # var -- can be On ro Off
    #        On  - 3D displayModes
    #        Off - 2D displaymodes
    if { ![info exists dispmode(mode3D_ModeFrame)] } {
	set dispmode(mode3D_ModeFrame) {}
    }

    if { $var == "Off" || $var == "On" } {
        set light $var     
    }   

    if { $var == "Off" } {
	# if 3D widgets exists delete them
	if {[winfo exists .ctrl.c.f.fr3.canv.mode3D]} { 
	    destroy .ctrl.c.f.fr3.canv.mode3D
	}
	if {[winfo exists $dispmode(mode3D_ModeFrame)]} {
	    #xcDisableAll $dispmode(mode3D_ModeFrame)
	    destroy $dispmode(mode3D_ModeFrame)
	}
	pack forget .ctrl.c.f.fr3.2
	    
	# if 2D widgets already exists, just return silently
	if {[winfo exists .ctrl.c.f.fr3.canv.mode2D]} { return }

	set f [frame .ctrl.c.f.fr3.canv.mode2D -bd 0]
	set ff [frame $f.2d]
	.ctrl.c.f.fr3.canv create window 0 0 \
	    -tags frame -anchor nw \
	    -width [winfo width .ctrl.c.f.fr3.canv] \
	    -height [winfo height .ctrl.c.f.fr3.canv] -window $f
	set b1 [radiobutton $ff.wire -text "WireFrame" -image wireframes_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value WF \
		    -command [list Display2D WF]]
	set b2 [radiobutton $ff.point -text "PointLine" -image pointlines_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value PL \
		    -command [list Display2D PL]]
	set b3 [radiobutton $ff.pipe -text "Pipe&Ball" -image pipeballs_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value PB \
		    -command [list Display2D PB]]
	set b4 [radiobutton $ff.ball1 -text "BallStick-1" -image ballsticks1_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value BS1 \
		    -command [list Display2D BS1]] 
	set b5 [radiobutton $ff.ball2 -text "BallStick-2" -image ballsticks2_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value BS2 \
		    -command [list Display2D BS2]]
	set b6 [radiobutton $ff.fill  -text "SpaceFill" -image spacefills_2d \
		    -indicatoron 0 -variable dispmode(mode2D_name) -selectcolor \#44ff44 -value SF \
		    -command [list Display2D SF]]
	#set i15 [expr round(15 * $xcMisc(resolution_ratio1))]
	#pack $b1 $b2 $b3 $b4 $b5 $b6 -side top -padx $i15 -pady 1 -fill x	
	grid $b1 -row 0 -column 0 -pady 1 -padx 1 -sticky nswe
	grid $b2 -row 0 -column 1 -pady 1 -padx 1 -sticky nswe
	grid $b3 -row 1 -column 0 -pady 1 -padx 1 -sticky nswe
	grid $b4 -row 1 -column 1 -pady 1 -padx 1 -sticky nswe
	grid $b5 -row 2 -column 0 -pady 1 -padx 1 -sticky nswe
	grid $b6 -row 2 -column 1 -pady 1 -padx 1 -sticky nswe

	foreach {wid text} [list \
				$b1 "WireFrame display-mode" \
				$b2 "PointLine display-mode" \
				$b3 "Pipe&Ball display-mode" \
				$b4 "BallStick-1 display-mode" \
				$b5 "BallStick-2 display-mode" \
				$b6 "SpaceFill display-mode"] {
	    DynamicHelp::register $wid balloon $text
	}
	pack $ff -side top
	#new
	set _h [expr [winfo reqheight $b1] * 3 + 10]
	.ctrl.c.f.fr3.canv configure -height $_h
	.ctrl.c.f.fr3.canv itemconfigure frame -height  $_h
	#/
	
	bind $b4 <Visibility> ConfigControlWid 
	From3Dto2D    	 
    } elseif { $var == "On" } {
	# if 2D widgets exists delete them
	if [winfo exists .ctrl.c.f.fr3.canv.mode2D] { 
	    destroy .ctrl.c.f.fr3.canv.mode2D 
	}
	# if 3D widgets already exists, just return silently
	if {[winfo exists .ctrl.c.f.fr3.canv.mode3D]} { return }
	
	eval pack .ctrl.c.f.fr3.2 $dispmode(mode3D_f2_packinfo) -after .ctrl.c.f.fr3.1

	set f [frame .ctrl.c.f.fr3.canv.mode3D -bd 0]
	.ctrl.c.f.fr3.canv create window 0 0 \
	    -tags frame -anchor nw \
	    -width [winfo width .ctrl.c.f.fr3.canv] \
	    -height [winfo height .ctrl.c.f.fr3.canv] -window $f 
	
	if ![info exists dispmode(mode3D)] {
	    set dispmode(mode3D) "Preset"
	}
	global radio_but_cmd_frame
	#if { ![winfo exists $dispmode(mode3D_ModeFrame)] } {
	#    ;
	#} else {
	#    xcEnableAll $dispmode(mode3D_ModeFrame)
	#}
	set rbs [RadioButCmd .ctrl.c.f.fr3.2 "Mode:" dispmode(mode3D) \
		     Mode3D top left 1 1 5 "Preset" "Logic"]
	foreach r [lrange $rbs 1 end] {
	    $r config -indicatoron 0 -selectcolor \#4444ff -font $xcFonts(small)
	    pack configure $r -padx 1 -pady 3 -ipadx 5 -ipady 1 -fill x -expand 1
	}
	set dispmode(mode3D_ModeFrame) $radio_but_cmd_frame
	#set rb1 [lindex $rbs 1]
	#set rb2 [lindex $rbs 2]
	#$rb1 configure -font $xcFonts(small)
	#$rb2 configure -font $xcFonts(small)

	# display widgets for $dispmode(mode3D)
	Mode3D $dispmode(mode3D)	
	# what displayMode3D to render
	From2Dto3D
    } else {
	puts stderr "invalid value \"$var\" submitted to proc DispMode, must be one of On, Off"
	flush stderr
        return
    }

    # add the additional created widgets to the control-toolbar scrolling
    scrollControlToolboxCmd .ctrl.c
    return
}


#==============================================================
# instead of xc_displayMode2D this proc should be used !!!!!!!!
#            ^^^^^^^^^^^^^^^^
proc Display2D {mode} {
    global mode2D dispmode

    set mode2D(WF)  Off
    set mode2D(PL)  Off
    set mode2D(PB)  Off
    set mode2D(BS1) Off
    set mode2D(BS2) Off
    set mode2D(SF)  Off

    if { $mode == "WF" } {
	set mode2D(WF) On
	set dispmode(mode2D_name) WF
	xc_displayMode2D .mesa WF
    } elseif { $mode == "PL" } {
	set mode2D(PL) On
	set dispmode(mode2D_name) PL
	xc_displayMode2D .mesa PL
    } elseif { $mode == "PB" } {
	set mode2D(PB) On
	set dispmode(mode2D_name) PB
	xc_displayMode2D .mesa PB
    } elseif { $mode == "BS1" } {
	set mode2D(BS1) On
	set dispmode(mode2D_name) BS1
	xc_displayMode2D .mesa BS1
    } elseif { $mode == "BS2" } {
	set mode2D(BS2) On
	set dispmode(mode2D_name) BS2
	xc_displayMode2D .mesa BS2
    } elseif { $mode == "SF" } {
	set mode2D(SF) On
	set dispmode(mode2D_name) SF
	xc_displayMode2D .mesa SF
    } else {
	puts stderr "invalid value \"$mode\" submitted to proc Display2D, must be one of WF, PL, PB, BS1, BS2, SF"
	flush stderr
    }
    return
}


proc Mode3D {mode} {
    global mode3D xcMisc

    set f .ctrl.c.f.fr3.canv.mode3D
    
    if { $mode == "Preset" } {
	# if "Boolean/Logic" buttons CheckButtons exists, delete them
	if { [winfo exists $f.bool] } { 
	    destroy $f.bool
	}
	if { [winfo exists $f.over] } {
	    return
	}
	set f [frame $f.over]
	set b1 [radiobutton $f.over1 -text "Stick" -image sticks_3d \
		    -indicatoron 0 -variable dispmode(mode3D_name) -selectcolor \#44ff44 -value S \
		    -command [list DisplayOver3D S]]
	set b2 [radiobutton $f.over1a -text "Pipe&Ball" -image pipeballs_3d \
		    -indicatoron 0 -variable dispmode(mode3D_name) -selectcolor \#44ff44 -value PB \
		    -command [list DisplayOver3D PB]]
	set b3 [radiobutton $f.over2 -text "BallStick" -image ballsticks_3d \
		    -indicatoron 0 -variable dispmode(mode3D_name) -selectcolor \#44ff44 -value BS \
		    -command [list DisplayOver3D BS]]
	set b4 [radiobutton $f.over3 -text "SpaceFill" -image spacefills_3d \
		    -indicatoron 0 -variable dispmode(mode3D_name) -selectcolor \#44ff44 -value SF \
		    -command [list DisplayOver3D SF]] 
	#set i15 [expr round(15 * $xcMisc(resolution_ratio1))]
	#pack $b1 $b2 $b3 $b4 -side top -padx $i15 -pady 1 -fill x -ipadx 0 -ipady 0
	grid $b1 -row 0 -column 0 -pady 1 -padx 1 -sticky nswe
	grid $b2 -row 0 -column 1 -pady 1 -padx 1 -sticky nswe
	grid $b3 -row 1 -column 0 -pady 1 -padx 1 -sticky nswe
	grid $b4 -row 1 -column 1 -pady 1 -padx 1 -sticky nswe

	pack $f -side top
	foreach {wid text} [list \
				$b1 "Stick display-mode" \
				$b2 "Pipe&Ball display-mode" \
				$b3 "BallStick display-mode" \
				$b4 "SpaceFill display-mode"] {
	    DynamicHelp::register $wid balloon $text
	}

	bind $b4 <Visibility> ConfigControlWid 
	#new
	set _h [expr [winfo reqheight $b1] * 2 + 10]
	.ctrl.c.f.fr3.canv configure -height $_h
	.ctrl.c.f.fr3.canv itemconfigure frame -height  $_h
	#/

    } elseif { $mode == "Logic" } {
	# Logic stands for Boolean
	if { [winfo exists $f.over] } { 
	    destroy $f.over
	}
	# create frame where CheckButtons will be displayed
	set ff [frame $f.bool]
	pack $ff -side top
	if { ![info exists mode3D(pipe)]   } { set mode3D(pipe)   "Off" }
	if { ![info exists mode3D(sticks)] } { set mode3D(sticks) "Off" }
	if { ![info exists mode3D(balls)]  } { set mode3D(balls)  "Off" }
	if { ![info exists mode3D(space)]  } { set mode3D(space)  "Off" }
	
	set f1 [frame $ff.1]
	set f2 [frame $ff.2]
	set f3 [frame $ff.3]
	checkbutton $f1.cb -variable mode3D(sticks) \
	    -onvalue On -offvalue Off -text "STICK" \
	    -command [list DispBool3D mode3D(sticks)] -anchor sw 

	checkbutton $f2.cb -variable mode3D(balls) \
	    -onvalue On -offvalue Off -text "BALL" \
	    -command [list DispBool3D mode3D(balls)] -anchor sw 

	checkbutton $f3.cb -variable mode3D(space) \
	    -onvalue On -offvalue Off -text "SPACEFILL" \
	    -command [list DispBool3D mode3D(space)] -anchor sw 

	pack $f1.cb $f2.cb $f3.cb -side top -fill both \
	    -padx 0 -pady 0 -ipadx 0 -ipady 0
	pack $f1 $f2 $f3 -side top -fill x -expand 1
	bind $f3 <Visibility> ConfigControlWid 
	#new
	set _h [expr [winfo reqheight $f1.cb] * 3 + 10]
	.ctrl.c.f.fr3.canv configure -height $_h
	.ctrl.c.f.fr3.canv itemconfigure frame -height  $_h
	#/

    } else {
	ErrorDialog "invalid value \"$mode\" submitted to proc Mode3D, must be one of Preset, Boolean"
    }

    # add the additional created widgets to the control-toolbar scrolling
    scrollControlToolboxCmd .ctrl.c
    return
}


# ============================================================================
# istead of "xc_displayMode3D over" this proc should be used
#            ^^^^^^^^^^^^^^^^^^^^^
proc DisplayOver3D {mode} {
    global mode3D dispmode

    set mode3D(sticks) Off
    set mode3D(balls)  Off
    set mode3D(space)  Off
    
    if { $mode == "S" } {
	set mode3D(sticks) On
	set dispmode(mode3D_name) S
	xc_displayMode3D .mesa over S
    } elseif { $mode == "PB" } {	
	set mode3D(pipe)   On
	set mode3D(sticks) On
	set mode3D(balls)  On
	set dispmode(mode3D_name) PB
	xc_displayMode3D .mesa over PB
    } elseif { $mode == "BS" } {
	set mode3D(pipe)   Off
	set mode3D(sticks) On
	set mode3D(balls)  On
	set dispmode(mode3D_name) BS
	xc_displayMode3D .mesa over BS
    } elseif { $mode == "SF" } {
	set mode3D(space) On
	set dispmode(mode3D_name) SF
	xc_displayMode3D .mesa over SF
    } else {
	ErrorDialog "invalid value \"$mode\" submitted to proc DisplayOver3D, must be one of S, PB, BS, SF"
    }
    return
}


# ============================================================================
# istead of "xc_displayMode3D boolean" this proc should be used
#            ^^^^^^^^^^^^^^^^^^^^^^^^
proc DispBool3D {var} {
    global mode3D
    
    if { $var == "mode3D(sticks)" } {
	set mode3D(space) Off
	Invert3D mode3D(sticks)
	xc_displayMode3D .mesa boolean S
    } elseif { $var == "mode3D(balls)" } {
	set mode3D(space) Off
	Invert3D mode3D(balls)
	xc_displayMode3D .mesa boolean B
    } elseif { $var == "mode3D(space)" } {
	set mode3D(sticks) Off
	set mode3D(balls)  Off
	Invert3D mode3D(space)
	xc_displayMode3D .mesa boolean SF
    }
    From3Dto3D
    return
}

proc From3Dto3D {} {
    global mode3D dispmode
    
    if { $mode3D(sticks) == "On" && $mode3D(balls) == "Off" } { 
	set dispmode(mode3D_name) S
    } elseif { $mode3D(balls) == "On" && $mode3D(pipe) =="Off" } {
	set dispmode(mode3D_name) BS
    } elseif { $mode3D(balls) == "On" && $mode3D(pipe) == "On" } {
	set dispmode(mode3D_name) PB
    } else {
	set dispmode(mode3D_name) SF
    }
    return
}

# this proc is for DispBool3D
proc Invert3D {var} {
    upvar #0 $var value
    if { $value == "On" } {
	set $var Off
    } else {
	set $var On
    }
}
    


# this proc map from 2D to 3D and display in 3D
proc From2Dto3D {} {
    global mode2D mode3D dispmode
    
    # mapping
    # (WF & PL)   --> Sticks (B)
    # (BS1 & BS2) --> BallSticks (BS)

    if ![info exists dispmode(mode3D)] { set dispmode(mode3D) Preset }

    if { $mode2D(WF) == "On" } {
	set mode3D(sticks) On 
	set dispmode(mode3D_name) S
	xc_displayMode3D .mesa over S
    } elseif { $mode2D(PL) == "On" || $mode2D(PB) == "On" } {
	set mode3D(pipe)   On
	set mode3D(sticks) On
	set mode3D(balls)  On
	set dispmode(mode3D_name) PB
	xc_displayMode3D .mesa over PB
    } elseif { $mode2D(BS1) == "On" || $mode2D(BS2) == "On" } {
	set mode3D(pipe)   Off
	set mode3D(sticks) On
	set mode3D(balls)  On
	set dispmode(mode3D_name) BS
	xc_displayMode3D .mesa over BS
    } else {
	set mode3D(sticks) Off
	set mode3D(balls)  Off
	set mode3D(space)  On
	set dispmode(mode3D_name) SF
	xc_displayMode3D .mesa over SF
    }

    set mode2D(WF)  Off
    set mode2D(PL)  Off
    set mode2D(PB)  Off
    set mode2D(BS1) Off
    set mode2D(BS2) Off
    set mode2D(SF)  Off

    # enable menubuttons for shade-model & draw-style
    set dispmode(style)   3D
    xcUpdateState
    return
}
	

# this proc map from 3D to 2D and display in 2D
proc From3Dto2D {} {
    global mode2D mode3D dispmode
    # mapping:
    # (S || ~B) --> WF
    # (B) --> BS2
    # (SF) --> BS2
    
    # if dispmode(mode3D) do not exists, just return silenlty, because
    # we were not yet on 3D modes and thus going from 3D to 2D is absurd
    if ![info exists dispmode(mode3D)] { 
	puts stdout "3D NOT EXIST"
	return }
	
    if { $mode3D(sticks) == "On" && $mode3D(balls) == "Off" } { 
	set mode2D(WF)   On 
	set dispmode(mode2D_name) WF
	xc_displayMode2D .mesa WF
    } elseif { $mode3D(balls) == "On" && $mode3D(pipe) =="Off" } {
	set mode2D(BS2)  On
	set dispmode(mode2D_name) BS2
	xc_displayMode2D .mesa BS2
    } elseif { $mode3D(balls) == "On" && $mode3D(pipe) == "On" } {
	set mode2D(PB)   On
	set dispmode(mode2D_name) PB
	xc_displayMode2D .mesa PB
    } else {
	set mode2D(SF)   On
	set dispmode(mode2D_name) SF
	xc_displayMode2D .mesa SF
    }
    set mode3D(pipe)   Off
    set mode3D(sticks) Off
    set mode3D(balls)  Off
    set mode3D(space)  Off

    # enable menubuttons for shade-model & draw-style
    set dispmode(style)   2D
    xcUpdateState
    return
}


# this proc deletes all child-widgets of parent
proc DeleteWid {parent} {
    
    set child [pack slaves $parent]
    puts stdout "CHILD = $child\n"
    foreach widget $child {
	puts stdout "WIDGET = $widget\n"
	destroy $widget
    }
}


# reset mode2D and mode3D arrays
proc ResetDispModes {} {
    global mode2D mode3D light

    set light Off
    set mode2D(WF)  Off
    set mode2D(PL)  On
    set mode2D(BS1) Off
    set mode2D(BS2) Off
    set dispmode(mode2D_name) PL


    set mode3D(sticks) Off
    set mode3D(balls)  Off
    set mode3D(space)  Off
    if { [info exists dispmode(mode3D)] } { unset dispmode(mode3D) }
}


proc ConfigControlWid {} {
    global widsize

    set width  [expr [winfo width .ctrl] - $widsize(scrW)]
    set height [winfo height .ctrl]
    set h      [winfo reqheight .ctrl.c.f]
    .ctrl.c config -width $width -height $height \
	    -scrollregion "0 0 $width $h" 
}


proc DisplayDefaultMode {} {
    global light

    set light On
    Lighting On
    Mode3D Preset
    DisplayOver3D BS
}


proc DisplayMode3D {} {
    global light
    set light On
    Lighting On
}


proc DisplayMode2D {} {
    global light
    set light Off
    Lighting Off
}
