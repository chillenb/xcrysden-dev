#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/hbonds.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc Hbonds {togl} {
    global check 

    if { $check(Hbonds) } {
	xc_hbonds $togl on
	$togl render
	xcSwapBuffers
    } else {
	xc_hbonds $togl off
	$togl render
	xcSwapBuffers
    }
}


proc HbondsSetting {togl} {
    global check Hbonds

    set t .hbond[lindex [split $togl .] end]

    if { [winfo exists $t] } {
	return
    }

    _HbondsSetting:Init $togl


    # from here on: widget managing
    xcToplevel $t "H-bonds: Settings" "H-bonds" . -0 0 1
    set f1 [frame $t.f1]
    set f2 [frame $t.f2]
    set f3 [frame $t.f3]
    pack $f1 $f2 $f3 -side top -expand 1 -fill x -padx 3m -pady 3m

    # frame-1: simple entries

    FillEntries $f1 {
	"H-like atom list:"
	"O-like atom list:"
	"Minimum H-bond length:"
	"Maximum H-bond length:"
	"Minimum H-bond angle:"
	"H-bond line width:"
	"H-bond line pattern:"
	"H-bond pattern size:"
    } {
	Hbonds(H_like_list) 
	Hbonds(O_like_list) 
	Hbonds(length_min)  
	Hbonds(length_max)  
	Hbonds(angle_min)   
	Hbonds(line_width)   
	Hbonds(line_pattern)   
	Hbonds(line_patternsize)   
    } 23 20 top left

    # frame-2: H-bond color
    
    xcModifyColor $f2 "H-bonds color:" [rgb_f2h $Hbonds(color)] \
	groove left left 100 100 100 5 20 $Hbonds(colorID)

    # frame-3: Close/Update buttons 

    set b1 [button $f3.close  -text "Close"  -command [list DestroyWid $t]]
    set b2 [button $f3.update -text "Update" -command [list HbondsSetting:update $togl]]
    
    pack $b1 $b2 -side left -padx 3m -pady 3m -expand 1
}


proc HbondsSetting:update {togl} {
    global Hbonds
 
    # check the values of the HBonds array here ...
    # ...insert...
    #--

    if { [info exists Hbonds(colorID)] } {
	set Hbonds(color) [xcModifyColorGet $Hbonds(colorID) float RGB]
    }

    xc_hbonds $togl set \
	-H_like_list  $Hbonds(H_like_list) \
	-O_like_list  $Hbonds(O_like_list) \
	-color        $Hbonds(color) \
	-length_min   $Hbonds(length_min) \
	-length_max   $Hbonds(length_max) \
	-angle_min    $Hbonds(angle_min) \
	-line_width   $Hbonds(line_width) \
	-line_pattern $Hbonds(line_pattern) \
	-line_patternsize $Hbonds(line_patternsize)

    $togl render
    xcSwapBuffers
}


proc _HbondsSetting:Init {togl} {
    global Hbonds
    foreach option {
	H_like_list
	O_like_list
	color
	length_min
	length_max
	angle_min
	line_width
	line_pattern
	line_patternsize
    } {
	if { ![info exists Hbonds($option)] } {
	    set Hbonds($option) [xc_hbonds $togl get -$option]
	}
    }

    if { ! [info exists Hbonds(colorID)] } {
	set Hbonds(colorID) [xcModifyColorID]
    }
}
