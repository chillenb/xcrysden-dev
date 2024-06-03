#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/planeselect.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc IsoPlaneGrid w {
    global isoplane

    #
    # initialize "xc_isospacesel" command
    #
    xc_isospacesel .mesa init

    set fp [frame $w.fp]
    set fp1 [frame $w.fp.1 -relief groove -bd 2]
    set l1  [label $fp1.l  \
	    -text "Select Plane for Property Evaluation" -relief flat]

    #
    # PARALLELOGRAM SELECTION
    #
    proc IsoGrid_ParaPlaneSel t {
	global isoplane xcColors select
	
	IsoPlaneGrid_State $t

	set select(done) 0
	set selw .igpps
	if [winfo exists .igpps] {
	    # we already have this toplevel displayed
	    return
	}
	PreSel $selw .mesa "Select Parallelogram" "Select a three-atoms spanning parallelogram. Please click on three atoms" ParalleSel 3
	tkwait window $selw
	if !$select(done) {
	    set isoplane(plane_sel) {}
	    IsoPlaneGrid_State $t
	    return
	}
	#
	# here I will have to set the four corners
	# isoplane(X/Y/Z,#)
	for {set i 1} {$i <= 3} {incr i} {
	    set isoplane(1,$i)  $select(X$i)
	    set isoplane(2,$i)  $select(Y$i)
	    set isoplane(3,$i)  $select(Z$i)
	}
	#  1---4
	#  |   |
	#  2---3
	#
	# point4 == point3 + point1 - point2
	foreach i {1 2 3} {
	    set isoplane($i,4) [expr $isoplane($i,3) + \
		    $isoplane($i,1) - $isoplane($i,2)]
	}
	#
	# END of IsoGrid_ParaPlaneSel
    }

    if ![info exists isoplane(plane_sel)] {
	set isoplane(plane_sel) {}
    }
    set r1  [radiobutton $fp1.r1 \
	    -text "Three-atoms spanning parallelogram selection" \
	    -variable isoplane(plane_sel) \
	    -value "para" \
	    -width 41 \
	    -command [list IsoGrid_ParaPlaneSel $w] \
	    -anchor w]
    pack $fp  -side top -fill both
    pack $fp1 -side top -fill both -padx 2 -pady 5
    
    grid $l1 -row 0 -column 0 -columnspan 3 -sticky we
    grid $r1 -row 1  -column 0 -columnspan 3 -sticky we

    set isoplane(AB_margin) 0.0
    set isoplane(CD_margin) 0.0
    set isoplane(AD_margin) 0.0
    set isoplane(BC_margin) 0.0
    set f [frame $fp1.f]
    grid $f -row 2 -column 1 -rowspan 4 -sticky wens
    FillEntries $f {{AB margin:} {CD margin:} {AD margin:} {BC margin:}} \
	    {isoplane(AB_margin) isoplane(CD_margin) isoplane(AD_margin) \
	    isoplane(BC_margin)} 10 10 
    set c11 [checkbutton $fp1.c11 \
	    -text "Rectangular parallelogram" \
	    -variable isoplane(rectangular) \
	    -command [list IsoPlaneGrid_Update $w] \
	    -anchor w]
    grid $c11 -row 6 -column 1 -columnspan 3 -rowspan 2 -sticky wen

    #
    # CENTERED SELECTION
    #
    proc IsoGrid_CenterPlaneSel t {
	global isoplane xcColors  isogrid_center
	
	IsoPlaneGrid_State $t

	if [winfo exists .igcps] {
	    # we already have this toplevel displayed
	    return
	}
	set isogrid_center {}
	set tp [xcToplevel .igcps "Atom centered selection" \
		"Atom centered selection" . -0 0 1]
	set f [frame $tp.f -class StressText]
	pack $f -padx 5 -pady 5
	set m1  [message $f.l  \
		-text "Atom centered selection is composed of two tasks:\n     1.) select an atom center\n     2.) Select a plane direction" \
		-aspect 500 \
		-justify left \
		-relief groove -bd 2]
	pack $m1 -side top

	proc IsoGrid_SelCenter {} {
	    global select isoplane xcColors
	    set select(done) 0
	    set selw [WidgetName]
	    PreSel $selw .mesa "Select an Atom Center" "For selecting an atom center click one atom" AtomSel 1
	    tkwait window $selw
	    if !$select(done) {
		return
	    }
	    #
	    # here get the data
	    #
	    set isoplane(c,X) $select(X1)
	    set isoplane(c,Y) $select(Y1)
	    set isoplane(c,Z) $select(Z1)

	    set isoplane(but_center) "Select an atom center\n(done)"
	    $isoplane(b1) config -background $xcColors(normal_bg)
	    if { $isoplane(but_center) == "Select an atom center\n(done)" && \
		   $isoplane(but_nml) == "Select a plane direction\n(done)" } {
		$isoplane(ok) config -state normal
	    }
	}

	proc IsoGrid_SelDirection {} {
	    global select isoplane xcColors
	    set select(done) 0
	    set selw [WidgetName]
	    PreSel $selw .mesa "Select a Plane Direction" "For selecting a plane direction click three atoms that determine the desired direction" \
		    ParalleSel 3 
	    tkwait window $selw
	    if !$select(done) {
		return
	    }

	    # 1
	    # |
	    # 2--3; vec1 = p1 - p2, vec2 = p3 - p2	   
	    # nml = vec1 x vec2; vec = v(X/Y/Z,#)
	    foreach i {X Y Z} {
		set v($i,1) [expr $select(${i}1) - $select(${i}2)]
		set v($i,2) [expr $select(${i}3) - $select(${i}2)]
	    }
	    set isoplane(nml,1) [expr $v(Y,1) * $v(Z,2) - $v(Y,2) * $v(Z,1)]
	    set isoplane(nml,2) [expr $v(Z,1) * $v(X,2) - $v(Z,2) * $v(X,1)]
	    set isoplane(nml,3) [expr $v(X,1) * $v(Y,2) - $v(X,2) * $v(Y,1)]

	    
	    set isoplane(but_nml) "Select a plane direction\n(done)"
	    $isoplane(b2) config -background $xcColors(normal_bg)
	    if { $isoplane(but_center) == "Select an atom center\n(done)" && \
		   $isoplane(but_nml) == "Select a plane direction\n(done)" } {
		$isoplane(ok) config -state normal
	    }
	}

	set isoplane(but_center) "Select an atom center\n(to be done)"
	set isoplane(but_nml)    "Select a plane direction\n(to be done)"
	set isoplane(b1) [button $f.b1 -textvariable isoplane(but_center) \
		-bg "#ffffdd" -command IsoGrid_SelCenter]
	set isoplane(b2) [button $f.b2 -textvariable isoplane(but_nml) \
		-bg "#ffffdd" -command IsoGrid_SelDirection]
	pack $isoplane(b1) $isoplane(b2) -side left -padx 5 -pady 5

	proc IsoGrid_CenterCan {t t1} {
	    global isogrid_center isoplane
	    set isoplane(plane_sel) {}
	    IsoPlaneGrid_State $t1
	    unset isogrid_center
	    CancelProc $t
	}
	proc IsoGrid_CenterOK {t w} {
	    global isogrid_center
	    # check if both steps (center & direction) have been determined
	    destroy $t
	    IsoPlaneGrid_Update $w
	    unset isogrid_center
	}

	set can [button $tp.can -text Cancel \
		-command [list IsoGrid_CenterCan $tp $t]]
	set ok [DefaultButton $tp.ok -text OK \
		-command [list IsoGrid_CenterOK $tp $t]]
	pack $can $ok -side right -padx 10 -pady 5
	set isoplane(ok) $ok.b
	$isoplane(ok) config -state disabled
	#
	# END of IsoGrid_CenterPlaneSel
    }

    set r2  [radiobutton $fp1.r2 \
	    -text "Atom centered selection" \
	    -variable isoplane(plane_sel) \
	    -value "center" \
	    -width 41 \
	    -command [list IsoGrid_CenterPlaneSel $w] \
	    -anchor w]
    if ![info exists isoplane(r)] {
	set isoplane(r) 1.0
    }
    if ![info exists isoplane(rot)] {
	set isoplane(rot) 0.0
    }
    grid $r2  -row 8 -column 0 -columnspan 3 -sticky we
    set f [frame $fp1.ff]
    grid $f -row 9 -column 1 -rowspan 2 -sticky wens
    FillEntries $f {radius: rotate:} {isoplane(r) isoplane(rot)} 10 10     
    set rot [button $fp1.rot -text "rotate" \
	    -command [list IsoPlaneGrid_Update $w 1]]
    grid $rot -row 10 -column 2 -sticky s

    frame $fp1.fupd -relief raised -bd 1
    set upd [button $fp1.upd -text "Update Display" -bd 1 \
	    -command [list IsoPlaneGrid_Update $w]]
    grid $fp1.fupd -row 11 -column 0 -columnspan 3 -sticky we -padx 3 -pady 3
    pack $upd -in $fp1.fupd -pady 3
    #
    # get the correct state
    #
    IsoPlaneGrid_State $w
}


proc IsoPlaneGrid_State t {
    global isoplane xcColors
    
    if { $isoplane(plane_sel) == "para" } {
	foreach i {1 2 3 4} {
	    $t.fp.1.f.f1.$i.entry$i config -state normal -relief sunken
	    $t.fp.1.f.f1.$i.lab$i config -foreground $xcColors(enabled_fg)
	}
	$t.fp.1.c11 config -state normal
	$t.fp.1.ff.f1.1.entry1 config -state disabled -relief flat
	$t.fp.1.ff.f1.1.lab1 config -foreground $xcColors(disabled_fg)
	$t.fp.1.ff.f1.2.entry2 config -state disabled -relief flat
	$t.fp.1.ff.f1.2.lab2 config -foreground $xcColors(disabled_fg)
	$t.fp.1.rot config -state disabled
	$t.fp.1.upd config -state normal
    } elseif { $isoplane(plane_sel) == "center" } {
	foreach i {1 2 3 4} {
	    $t.fp.1.f.f1.$i.entry$i config -state disabled -relief flat
	    $t.fp.1.f.f1.$i.lab$i config -foreground $xcColors(disabled_fg)
	}
	$t.fp.1.c11 config -state disabled
	$t.fp.1.ff.f1.1.entry1 config -state normal -relief sunken
	$t.fp.1.ff.f1.1.lab1 config -foreground $xcColors(enabled_fg)
	$t.fp.1.ff.f1.2.entry2 config -state normal -relief sunken
	$t.fp.1.ff.f1.2.lab2 config -foreground $xcColors(enabled_fg)
	$t.fp.1.rot config -state normal
	$t.fp.1.upd config -state normal
    } else {
	foreach i {1 2 3 4} {
	    $t.fp.1.f.f1.$i.entry$i config -state disabled -relief flat
	    $t.fp.1.f.f1.$i.lab$i config -foreground $xcColors(disabled_fg)
	}
	$t.fp.1.c11 config -state disabled
	$t.fp.1.ff.f1.1.entry1 config -state disabled -relief flat
	$t.fp.1.ff.f1.1.lab1 config -foreground $xcColors(disabled_fg)
	$t.fp.1.ff.f1.2.entry2 config -state disabled -relief flat
	$t.fp.1.ff.f1.2.lab2 config -foreground $xcColors(disabled_fg)
	$t.fp.1.rot config -state disabled
	$t.fp.1.upd config -state disabled
    }
}


proc IsoPlaneGrid_Update {w {rot 0}} {
    global isoplane

    if { $isoplane(plane_sel) == "para" } {
	set varlist { {isoplane(AB_margin) real} {isoplane(CD_margin) real} \
		{isoplane(AD_margin) real} {isoplane(BC_margin) real} }
	set foclist [list \
		$w.fp.1.f.f1.1.entry1 $w.fp.1.f.f1.2.entry2 \
		$w.fp.1.f.f1.3.entry3 $w.fp.1.f.f1.4.entry4]	
	if ![check_var $varlist $foclist] {
	    return;
	}
	set isoplane(points) [xc_isospacesel .mesa paralle2D \
		-points [list \
		$isoplane(1,1) $isoplane(2,1) $isoplane(3,1) \
		$isoplane(1,2) $isoplane(2,2) $isoplane(3,2) \
		$isoplane(1,3) $isoplane(2,3) $isoplane(3,3) \
		$isoplane(1,4) $isoplane(2,4) $isoplane(3,4)] \
		-margins [list \
		$isoplane(AB_margin) $isoplane(CD_margin) \
		$isoplane(AD_margin) $isoplane(BC_margin)] \
		-rectangular $isoplane(rectangular)]
	# set the main OK button to normal state
	$isoplane(OK) config -state normal
	if { [winfo exists $isoplane(editButton)] } {
	    $isoplane(editButton) config -state normal
	}
    } elseif { $isoplane(plane_sel) == "center" } {
	set varlist { {isoplane(r) real} {isoplane(rot) real} }
	set foclist [list $w.fp.1.ff.f1.1.entry1 $w.fp.1.ff.f1.2.entry2]
	if ![check_var $varlist $foclist] {
	    return;
	}
	
	set rotate 0.0
	if $rot {
	    # rotate button was pressed
	    set rotate $isoplane(rot)
	}

	set isoplane(points) [xc_isospacesel .mesa center2D \
		-center [list  $isoplane(c,X) \
		$isoplane(c,Y) $isoplane(c,Z)] \
		-normal [list  $isoplane(nml,1) \
		$isoplane(nml,2) $isoplane(nml,3)] \
		-size $isoplane(r) \
		-rotate $rotate]
	# set the main OK button to normal state
	$isoplane(OK) config -state normal
	if { [winfo exists $isoplane(editButton)] } {
	    $isoplane(editButton) config -state normal
	}
    }
}
