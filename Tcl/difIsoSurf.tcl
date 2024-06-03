#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/difIsoSurf.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc DiffIsoSurf_Widget {{dim 3D}} {
    global prop dif_isosurf
    
    # if window already exist destroy it and create new one
    catch {destroy .isodif}
    set t [xcToplevel .isodif "Difference Maps" "Diff. Maps" . 100 50]

    set top   [frame $t.top       -relief raised -bd 2]
    frame $t.mid
    set left  [frame $t.mid.left  -relief raised -bd 2]
    set right [frame $t.mid.right -relief raised -bd 2]
    set bot   [frame $t.bot       -relief raised -bd 2]

    pack $top $t.mid $bot -side top -expand 1 -fill both
    pack $left $right -side left -expand 1 -fill both

    set l  [label $top.l -text "Difference Map == Map A  -  Map B"]
    set ll [label $left.ll  -text "Map A:" -relief ridge -bd 2]
    set lr [label $right.lr -text "Map B:" -relief ridge -bd 2]
    pack $l -padx 10 -pady 10 -expand 1
    pack $ll $lr -side top -padx 10 -pady 5 -ipadx 20 -ipady 5

    ########################################
    # make everything that goes in LEFT FRAME
    frame $left.f1 -relief groove -bd 2
    set cb [checkbutton $left.f1.cb -text "" -anchor w]
    $cb invoke
    $cb config -state disabled

    set l [label $left.f1.l -textvariable dif_isosurf(file_textvar_A) \
	    -relief flat -anchor w]
    pack $cb $l -side top -fill x -padx 5 -pady 5
    frame $left.f2 -relief groove -bd 2

    xcMenuEntry $left.f2 "Specify Property:" 30  \
	    dif_isosurf(prop_A) $prop(dif_prop3D_list) \
	    -labelwidth 16
    pack $left.f1 $left.f2 -padx 2 -pady 10 -ipady 3 -fill x

    RadioButVarCmd $left "Density matrix to use:" \
	    dif_isosurf(denmat_A) DiffIsoSurf_LeftCmd \
	    top top 1 1 \
	    "SCF density matrix" \
	    "Density matrix as superposition of atomic densities"

    checkbutton $left.f.cb \
	    -text "Give new basis set/electron configuration" \
	    -variable prop(PATO_newbasis_A) \
	    -onvalue 1 \
	    -offvalue 0 \
	    -anchor e
    $left.f.cb config -state disabled

    pack $left.f.cb
    DiffIsoSurf_LeftCmd $dif_isosurf(denmat_A)

    ########################################
    # make everything that goes in RIGHT FRAME
    frame $right.f1 -relief groove -bd 2

    set cb [checkbutton $right.f1.cb \
	    -text "Load Crystal95's unit 9 from another file" \
	    -variable dif_isosurf(load_another_unit9) \
	    -onvalue 1 -offvalue 0 \
	    -command LoadAnotherUnit9 \
	    -anchor w]

    set l [label $right.f1.l -textvariable dif_isosurf(file_textvar_B) \
	    -relief flat -anchor w]
    pack $cb $l -side top -fill x -padx 5 -pady 5
    frame $right.f2 -relief groove -bd 2

    xcMenuEntry $right.f2 "Specify Property:" 30  \
	    dif_isosurf(prop_B) $prop(dif_prop3D_list) \
	    -labelwidth 16
    pack $right.f1 $right.f2 -side top -padx 2 -pady 10 -ipady 3 -fill x

    RadioButVarCmd $right "Density matrix to use:" \
	    dif_isosurf(denmat_B) DiffIsoSurf_RightCmd \
	    top top 1 1 \
	    "SCF density matrix" \
	    "Density matrix as superposition of atomic densities"

    checkbutton $right.f.cb \
	    -text "Give new basis set/electron configuration" \
	    -variable prop(PATO_newbasis_B) \
	    -onvalue 1 \
	    -offvalue 0 \
	    -anchor e
    $right.f.cb config -state disabled
    pack $right.f.cb
    DiffIsoSurf_RightCmd $dif_isosurf(denmat_B)

    ########################################
    # bottom frame
    set can [button $bot.can -text "Cancel" -command [list CancelProc $t]]
    set ok  [button $bot.ok  -text "OK" \
	    -command [list DiffIsoSurf_WidgetOK $t $dim]]
    pack $can $ok -side left -expand 1 -padx 10 -pady 10
}



proc DiffIsoSurf_LeftCmd item {
    global prop

    if { $item == "SCF density matrix" } {
	#.isodif.mid.left.f.cb configure -state disabled
	set prop(PATO_newbasis_A) 0
    } else {
	#.isodif.mid.left.f.cb configure -state normal
    }
}


proc DiffIsoSurf_RightCmd item {
    global prop 

    if { $item == "SCF density matrix" } {
	#.isodif.mid.right.f.cb configure -state disabled
	set prop(PATO_newbasis_B) 0
    } else {
	#.isodif.mid.right.f.cb configure -state normal
    }
}


proc LoadAnotherUnit9 {} {
    global dif_isosurf fileselect

    # if we turn off the "load another Crystal95's unit 9"
    if { $dif_isosurf(load_another_unit9) == 0 } { 
	set dif_isosurf(file_textvar) "loaded UNIT 9:\nsame as in Map A"
	return 
    } 

    fileselect "Load another Crystal95's unit 9"

    # maybe CANCEL button was pressed
    if { $fileselect(path) == "" } {
	set dif_isosurf(file_textvar) "loaded UNIT 9:\nsame as in Map A"
	set dif_isosurf(load_another_unit9) 0
	return
    }
    
    set dif_isosurf(file_textvar_B) "loaded Unit 9:\n$fileselect(path)"
    set dif_isosurf(unit9_B) $fileselect(path)
}


proc DiffIsoSurf_WidgetOK {t {dim 3D}} {
    global prop dif_isosurf system
    
    set match_A 0
    set match_B 0
    foreach item $prop(dif_prop3D_list) {
	if { $item == [string toupper $dif_isosurf(prop_A)] } {
	    set match_A 1
	}
	if { $item == [string toupper $dif_isosurf(prop_B)] } {
	    set match_B 1
	}
    }

    # replace {,} with " in prop(dif_prop3D_lis) and write to $p_list
    regsub -all -- \{ $prop(dif_prop3D_list) \" p_list
    regsub -all -- \} $p_list \" p_list
    
    if { $match_A == 0 } {
	tk_dialog [WidgetName] ERROR "ERROR: unknown property_A \"$dif_isosurf(prop_A)\", should be of $p_list" error 0 OK
	focus .isodif.mid.left.f2.e
	return
    }
    if { $match_B == 0 } {
	tk_dialog [WidgetName] ERROR "ERROR: unknown property_B \"$dif_isosurf(prop_A)\", should be of $p_list" error 0 OK
	focus .isodif.mid.right.f2.e
	return
    }
	
    # we come so far, so everything is OK
    set dif_isosurf(prop_A) [string toupper $dif_isosurf(prop_A)]
    set dif_isosurf(prop_B) [string toupper $dif_isosurf(prop_B)]

    if [winfo exist $t] {destroy $t}

    # now map from $dif_isosurf(prop_X) to {ECHG ECHD POTM}
    set propA [C95Name2Com $dif_isosurf(prop_A)]
    set propB [C95Name2Com $dif_isosurf(prop_B)]
    
    # now space-selection must be specefied 
    SetIsoGridSpace "Difference 3D Map - Grid Specifications" "Grid" \
	    [list $propA $propB] $dim
}



proc C95Name2Com name {
    # thiw proc map from property-name to C95's propertiy-command

    if { $name == "CHARGE DENSITY" } {
	return ECHG\n0
    } elseif { $name == "CHARGE DENSITY GRADIENT" } {
	return ECHG\n1
    } elseif { $name == "ELECTROSTATIC POTENTIAL" } {
	return POTM
    }
    return 0
}
