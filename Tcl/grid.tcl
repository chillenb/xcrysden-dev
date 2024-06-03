#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/grid.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

###############################################################################
# some properties can be rendered as ALPHA/BETA/ALPHA+BETA/ALPHA-BETA if
# type_of_run == UHF. Such properties are ECHD & ECHG
# 
# for ECHD & ECHG first record is ALPHA+BETA and second one is: ALPHA-BETA
# so ALPHA = 0.5 * (record1 + record2)
#    BETA  = 0.5 * (record1 - record2)
###############################################################################

###############################################################################
# isosurf(res_type)........resolution type for isosurface, can be points, 
#                          angstroms, bohrs
# isosurf(resol_poi).......resolution of grid for isosurfaces in points
# isosurf(resol_ang).......resolution of grid for isosurfaces in angstroms
# isosurf(mb_angs/bohr)....are resolution specified in Angstoms/Bohrs
# isosurf(3Dinterpl_degree)...degree of tri cubic spline interpolation
# isosurf(space_sel).......how space is selected (whole_cell/manually)
# isosurf(Y_Sel)...........how Y-space is selected (centered/min_max/offset)
# isosurf(Y_size)..........size of Y-space (for $isosurf(Y_sel) == centered)
# isosurf(Y_min)...........minimun Y (for) $isosurf(Y_sel) == min_max
# isosurf(Y_max)...........maximum_Y
# isosurf(Y_offset)........Y offset for $isosurf(Y_sel) == "offset"
# isosurf(Z_Sel)
# isosurf(Z_size)
# isosurf(Z_min)
# isosurf(Z_max)
# isosurf(Z_offset)
##########
# isosurf(minvalue)........minimal grid value 
# isosurf(maxvalue)........maximum grid value
# isosurf(rangevalue)......[expr $isosurf(maxvalue) - $isosurf(minvalue)]
# isosurf(isolevel)........isovalue to render
# isosurf(type_of_surf)....type of isosurface; solid/wire
# isosurf(shade_model).....smooth/flat
# isosurf(twoside_lighting)...off/on
# isosurf(transparency).......off/on
##########
# isosurf(spin)............what spin to take for UHF type_of_run
##########
# isosurf(isovalue_entry)..this is name of entry-widget for specifying
#                          isolevel
##############################################################################

proc SetIsoGridSpace {title iconname comlist {dim 3D}} {
    global periodic isosurf varlist foclist isogrid_space prop \
	    dif_isosurf spin_entry isoplane

    set varlist {}
    set foclist {}

    set t [xcToplevel .isoget $title $iconname . 100 50 1]

    set f1 [frame $t.f1 -relief raised -bd 2]
    set f2 [frame $t.f2 -relief raised -bd 2]
    set f3 [frame $t.f3 -relief raised -bd 2]

    #this is the name of OK button at the bottom of the window
    set isoplane(OK)         $f3.ok
    set isoplane(editButton) $f3.edit

    set com1            [lindex $comlist 0]    
    set spinmatch       0
    set prop(spin_case) 0
    if { ! [info exists prop(spin_prop_list)] } {
	# for WIEN2k this variable doesn't exists
	set prop(spin_prop_list) NONEXISTENT
    }
    if { [string match "Difference*" $title] } {
	# we have difference map
	set com2 [lindex $comlist 1]
	set dif_isosurf(dif_map) 1
	# do we have to deal with spin-dependent properties
	foreach _com $prop(spin_prop_list) {
	    if { [string match "*$com1*" $_com] || \
		     [string match "*$com2*" $_com] } {
		set spinmatch 1
		break
	    }
	}
    } else {
	# we don't have difference map
	foreach _com $prop(spin_prop_list) {
	    if { [string match "*$com1*" $_com] } {
		set spinmatch 1
		break
	    }
	}
	set dif_isosurf(dif_map) 0
    }
    
    ###########################################################
    # for UHF there is a need to specify what spin to take, 
    # but just if we are dealing with spin dependent properties
    if { $prop(type_of_run) == "UHF" && $spinmatch } {
	#######################
	# we have "SPIN_CASE" #
	#######################
	set prop(spin_case) 1
	#######################
	set f0 [frame $t.f0 -relief raised -bd 2]
	frame $f0.1 -relief groove -bd 2
	set spin_entry [xcMenuEntry $f0.1 "What SPIN to take:" 30  \
		isosurf(spin) {ALPHA BETA ALPHA+BETA ALPHA-BETA} \
		-labelwidth 16]
	pack $f0 -side top -fill both
	pack $f0.1 -padx 2 -pady 5 -ipady 3 -fill x
    }

    pack $f1 $f2 $f3 -side top -fill both
    ########################################
    # SPECIFY GRID RESOLUTION
    set f11 [frame $f1.f1 -relief groove -bd 2]
    set l11 [label $f11.l  -text "Specify Grid Resolution" -relief flat]

    proc SetIsoGridSpace_com1 {} {
	focus .isoget.f1.f1.e1 
	.isoget.f1.f1.e1  configure -state normal   -relief sunken
	.isoget.f1.f1.e2  configure -state disabled -relief flat
	.isoget.f1.f1.mb1 configure -state disabled
    }
    proc SetIsoGridSpace_com2 {} { 
	focus .isoget.f1.f1.e2
	.isoget.f1.f1.e1  configure -state disabled -relief flat
	.isoget.f1.f1.e2  configure -state normal   -relief sunken
	.isoget.f1.f1.mb1 configure -state normal
    }

    set f111 [frame $f11.1 -relief flat]
    set r11 [radiobutton $f11.r1 \
	    -text "Number of points along one grid-segment:" \
	    -variable isosurf(res_type) \
	    -value "points" \
	    -width 41 \
	    -command SetIsoGridSpace_com1 \
	    -anchor w]
    set e11 [entry $f11.e1 \
	    -textvariable isosurf(resol_poi) \
	    -width 10]
    lappend varlist {isosurf(resol_poi) posint}
    lappend foclist $f11.e1
    
    set f112 [frame $f11.2 -relief flat]
    set r12 [radiobutton $f11.r2 \
	    -text "Resolution in Angstroms/Bohrs:" \
	    -variable isosurf(res_type) \
	    -value "angstroms" \
	    -width 41 \
	    -command SetIsoGridSpace_com2 \
	    -anchor w] 
    set e12 [entry $f11.e2 \
	    -textvariable isosurf(resol_ang) \
	    -width 10]
    lappend varlist {isosurf(resol_ang) posreal}
    lappend foclist $f11.e2

    set mb11 [menubutton $f11.mb1 \
	    -textvariable isosurf(mb_angs/bohr) \
	    -menu $f11.mb1.menu \
	    -indicatoron 1 \
	    -relief raised \
	    -width 9 \
	    -anchor w]
    set menu11 [menu $mb11.menu -relief raised -tearoff 0]
    $menu11 add command -label "Angstroms" \
	    -command [list set isosurf(mb_angs/bohr) "Angstroms"]
    $menu11 add command -label "Bohrs" \
	    -command [list set isosurf(mb_angs/bohr) "Bohrs"]

    if ![info exists isosurf(3Dinterpl_degree)] {
	set isosurf(3Dinterpl_degree) 1
	
    }
    if { $dim == "3D" } {
	set label "Degree of triCubic Spline Interpolation:"
    } elseif { $dim == "2D" } {
	set label "Degree of biCubic Spline Interpolation:"
    }
    set f113 [frame $f11.3 -relief flat]
    scale $f11.sc -from 1 -to 4 -length 100 \
	    -variable isosurf(3Dinterpl_degree) -orient horizontal \
	    -label $label \
	    -tickinterval 1 -resolution 1 \
	    -showvalue true 
    trace variable isosurf(3Dinterpl_degree) w xcTrace
    ##########
    # get init states
    if { $isosurf(res_type) == "points" } {
	SetIsoGridSpace_com1
    } else {
	SetIsoGridSpace_com2
    }

    pack $f11  -side top -fill both -padx 2 -pady 5
    pack $l11  -side top -expand 1

    pack $f111 -side top -padx 5 -fill x 
    pack $r11  -side left -in $f111
    pack $e11  -side left -in $f111 -fill x -expand 1

    pack $f112 -side top -padx 5 -fill x 
    pack $r12  -side left -in $f112 -pady 5
    pack $e12  -side left -in $f112 -pady 5 -fill x -expand 1
    pack $mb11 -side left -in $f112 -pady 5

    pack $f113 -side top -padx 5 -fill x 
    pack $f11.sc -side left -in $f113 -pady 5 -padx 5 -fill x -expand 1
    
    ########################################
    # SPACE SELECTIONS
    ########################################
    if { $dim == "3D" } {
	proc SetIsoGridSpace_com3 {} {
	    xcEnableAll {.isoget.f2.f1.mar .isoget.f2.f1.1 .isoget.f2.f1.2}
	    set p0 .isoget.f2.f1.mar
	    set p1 .isoget.f2.f1.1
	    set p2 .isoget.f2.f1.2
	    set p3 .isoget.f2.f1.3
	    foreach ent { e1 e2 e3 e4 e5 e6 } { 
		catch {$p0.$ent configure -relief sunken}
		catch {$p1.$ent configure -relief sunken}
		catch {$p2.$ent configure -relief sunken}
		catch {$p3.$ent configure -relief sunken}
	    }
	}
	set f21 [frame $f2.f1 -relief groove -bd 2]
	set l21 [label $f21.l  -text "Select Space for Isosurface Evaluation" \
		-relief flat]

	pack $f21  -side top -fill both -padx 2 -pady 5
	pack $l21  -side top -expand 1
	
	if { $periodic(dim) == 3 } {
	    ###########
	    # CRYSTAL #
	    ###########
	    set r21 [radiobutton $f21.r1 \
		    -text "Whole unit cell" \
		    -variable isosurf(space_sel) \
		    -value "whole_cell" \
		    -command SetIsoGridSpace_com3 \
		    -width 41 \
		    -anchor w]
	    pack $r21  -side top -fill x -expand 1 -padx 5
	    IsoCellSel 
	} elseif { $periodic(dim) < 3 && $periodic(dim) > 0 } {
	    ###################
	    #  SLAB & POLYMER #	
	    ###################	
	    if { $periodic(dim) == 2 } {
		set r21text "Whole unit cell in XY plane"
		set ii 1
	    }
	    if { $periodic(dim) == 1 } {
		set r21text "Whole unit cell along X line"
		set ii 2
	    }
	    	    
	    set r21 [radiobutton $f21.r1 \
		    -text $r21text  \
		    -variable isosurf(space_sel) \
		    -value "whole_cell" \
		    -width 41 \
		    -command SetIsoGridSpace_com3 \
		    -anchor w]
	    pack $r21 -side top -fill x -expand 1 -padx 5
	    IsoCellSel
	    
	    ##########
	    # get init states
	    if { $isosurf(space_sel) == "whole_cell" } {
		SetIsoGridSpace_com3
	    }
	    
	    if { $periodic(dim) == 1 } {
		####################
		# determine Y-space
		set f211  [frame $f21.1 -relief flat]
		set f211a [frame $f211.a -relief flat -width 40]
		set f211b [frame $f211.b -relief flat]
		
		set l211b0 [label $f211.l0 -text "Determine Y-space:" \
			-relief flat -anchor w]
		
		proc SetIsoGridSpace_com4 {} {
		    focus .isoget.f2.f1.1.e1
		    .isoget.f2.f1.1.e1 configure -state normal   -relief sunken
		    .isoget.f2.f1.1.e2 configure -state disabled -relief flat
		    .isoget.f2.f1.1.e3 configure -state disabled -relief flat 
		    .isoget.f2.f1.1.e4 configure -state disabled -relief flat
		}
		proc SetIsoGridSpace_com5 {} {
		    focus .isoget.f2.f1.1.e2
		    .isoget.f2.f1.1.e1 configure -state disabled -relief flat
		    .isoget.f2.f1.1.e2 configure -state normal   -relief sunken
		    .isoget.f2.f1.1.e3 configure -state normal   -relief sunken
		    .isoget.f2.f1.1.e4 configure -state disabled -relief flat
		}
		proc SetIsoGridSpace_com6 {} {		
		    focus .isoget.f2.f1.1.e4
		    .isoget.f2.f1.1.e1 configure -state disabled -relief flat
		    .isoget.f2.f1.1.e2 configure -state disabled -relief flat
		    .isoget.f2.f1.1.e3 configure -state disabled -relief flat
		    .isoget.f2.f1.1.e4 configure -state normal   -relief sunken
		}
		
		set f211c1 [frame $f211.c1 -relief flat]
		set r211c1 [radiobutton $f211.r1 \
			-text "set center of Y at Y=0.0;" \
			-variable isosurf(Y_Sel) \
			-value "centered" \
			-command SetIsoGridSpace_com4 \
			-anchor w]
		set l211c1 [label $f211.l1 -text "Y-size:" \
			-width 10 \
			-anchor e]
		set e211c1 [entry $f211.e1 \
			-textvariable isosurf(Y_size) \
			-width 20]
		lappend varlist {isosurf(Y_size) real}
		lappend foclist $f211.e1
		
		set f211c2 [frame $f211.c2 -relief flat -width 150]
		set r211c2 [radiobutton $f211.r2 \
			-text "min Y:" \
			-variable isosurf(Y_Sel) \
			-value "min_max" \
			-command SetIsoGridSpace_com5 \
			-anchor w ]
		set e211c2 [entry $f211.e2 \
			-textvariable isosurf(Y_min) \
			-width 20]
		lappend varlist {isosurf(Y_min) real}
		lappend foclist $f211.e2
		
		set l211c3 [label $f211.l2 -text " max Y:" \
			-width 10 \
			-anchor e]
		set e211c3 [entry $f211.e3 \
			-textvariable isosurf(Y_max) \
			-width 20]
		lappend varlist {isosurf(Y_max) real}
		lappend foclist $f211.e3
		
		set f211c4 [frame $f211.c4 -relief flat]
		set r211c4 [radiobutton $f211.r4 \
			-text "whole structure per cell with Y-offset;" \
			-variable isosurf(Y_sel) \
			-value "offset" \
			-anchor w \
			-command SetIsoGridSpace_com6]
		set l211c4 [label $f211.l4 -text "Y-offset:" \
			-anchor e  \
			-width 10]
		set e211c4 [entry $f211.e4 \
			-textvariable isosurf(Y_offset) \
			-width 20]
		lappend varlist {isosurf(Y_offset) real}	
		lappend foclist $e211c4
		
		pack $f211 -side top -fill x -expand 1 -padx 5
		pack $f211a -side left
		pack $f211b -side left -expand 1 -fill x
		pack $l211b0 -in $f211b -side top -expand 1 -fill x
		pack $f211c1 $f211c2 $f211c4 -in $f211b -side top -fill x
		pack $r211c1 -in $f211c1 -side left -fill x
		pack $l211c1 -in $f211c1 -side left -fill x -expand 1
		pack $e211c1 -in $f211c1 -side left  
		pack $r211c2 $e211c2 -in $f211c2 -side left -fill x
		pack $l211c3 -in $f211c2 -side left -fill x -expand 1
		pack $e211c3 -in $f211c2 -side left 
		pack $r211c4 $l211c4 -in $f211c4 -side left -fill x
		pack $e211c4 -in $f211c4 -side left
		
		##########
		# get init states
		if { $isosurf(Y_Sel) == "centered" } { SetIsoGridSpace_com4 }
		if { $isosurf(Y_Sel) == "min_max" }  { SetIsoGridSpace_com5 }
		if { $isosurf(Y_Sel) == "offset" }   { SetIsoGridSpace_com6 }
	    }
	    
	    ####################
	    # determine Z-space
	    set f211  [frame $f21.2 -relief flat]
	    set f211a [frame $f211.a -relief flat -width 40]
	    set f211b [frame $f211.b -relief flat]
	    
	    proc SetIsoGridSpace_com7 {} {
		focus .isoget.f2.f1.2.e1
		.isoget.f2.f1.2.e1 configure -state normal   -relief sunken
		.isoget.f2.f1.2.e2 configure -state disabled -relief flat
		.isoget.f2.f1.2.e3 configure -state disabled -relief flat
		.isoget.f2.f1.2.e4 configure -state disabled -relief flat
	    }
	    proc SetIsoGridSpace_com8 {} {
		focus .isoget.f2.f1.2.e2
		.isoget.f2.f1.2.e1 configure -state disabled -relief flat
		.isoget.f2.f1.2.e2 configure -state normal   -relief sunken
		.isoget.f2.f1.2.e3 configure -state normal   -relief sunken
		.isoget.f2.f1.2.e4 configure -state disabled -relief flat
	    }
	    proc SetIsoGridSpace_com9 {} {
		focus .isoget.f2.f1.2.e4
		.isoget.f2.f1.2.e1 configure -state disabled -relief flat
		.isoget.f2.f1.2.e2 configure -state disabled -relief flat
		.isoget.f2.f1.2.e3 configure -state disabled -relief flat
		.isoget.f2.f1.2.e4 configure -state normal   -relief sunken
	    }
	    
	    
	    set text "Determine Z-space:"
	    set l211b0 [label $f211.l0 -text "Determine Z-space:" \
		    -relief flat -anchor w]
	    
	    set f211c1 [frame $f211.c1 -relief flat]
	    set r211c1 [radiobutton $f211.r1 \
		    -text "set center of Z at Z=0.0;" \
		    -variable isosurf(Z_Sel) \
		    -value "centered" \
		    -anchor w \
		    -command SetIsoGridSpace_com7]
	    set l211c1 [label $f211.l1 -text "Z-size:" \
		    -width 10 \
		    -anchor e]
	    set e211c1 [entry $f211.e1 \
		    -textvariable isosurf(Z_size) \
		    -width 20]
	    lappend varlist {isosurf(Z_size) real}
	    lappend foclist $f211.e1
	    
	    set f211c2 [frame $f211.c2 -relief flat -width 150]
	    set r211c2 [radiobutton $f211.r2 \
		    -text "min Z:" \
		    -variable isosurf(Z_Sel) \
		    -value "min_max" \
		    -anchor w \
		    -command SetIsoGridSpace_com8]
	    set e211c2 [entry $f211.e2 \
		    -relief sunken \
		    -bd 2 \
		    -textvariable isosurf(Z_min) \
		    -width 20]
	    lappend varlist {isosurf(Z_min) real}
	    lappend foclist $f211.e2
	    
	    set l211c3 [label $f211.l2 -text " max Z:" \
		    -width 10 \
		    -anchor e]
	    set e211c3 [entry $f211.e3 \
		    -textvariable isosurf(Z_max) \
		    -width 20]
	    lappend varlist {isosurf(Z_max) real}
	    lappend foclist $f211.e3
	    
	    set f211c4 [frame $f211.c4 -relief flat]
	    set r211c4 [radiobutton $f211.r4 \
		    -text "whole structure per cell with Z-offset;" \
		    -variable isosurf(Z_Sel) \
		    -value "offset" \
		    -anchor w \
		    -command SetIsoGridSpace_com9]
	    set l211c4 [label $f211.l4 -text "Z-offset:" \
		    -anchor e \
		    -width 10]
	    set e211c4 [entry $f211.e4 \
		    -textvariable isosurf(Z_offset) \
		    -width 20]
	    lappend varlist {isosurf(Z_offset) real}	
	    lappend foclist $e211c4
	    
	    pack $f211 -side top -fill x -expand 1 -padx 5
	    pack $f211a -side left
	    pack $f211b -side left -expand 1 -fill x
	    pack $l211b0 -in $f211b -side top -expand 1 -fill x
	    pack $f211c1 $f211c2 $f211c4 -in $f211b -side top -fill x
	    pack $r211c1 -in $f211c1 -side left -fill x
	    pack $l211c1 -in $f211c1 -side left -fill x -expand 1
	    pack $e211c1 -in $f211c1 -side left  
	    pack $r211c2 $e211c2 -in $f211c2 -side left -fill x
	    pack $l211c3 -in $f211c2 -side left -fill x -expand 1
	    pack $e211c3 -in $f211c2 -side left 
	    pack $r211c4 $l211c4 -in $f211c4 -side left -fill x
	    pack $e211c4 -in $f211c4 -side left
	    
	    ##########
	    # get init states
	    if { $isosurf(Z_Sel) == "centered" } { SetIsoGridSpace_com7 }
	    if { $isosurf(Z_Sel) == "min_max" }  { SetIsoGridSpace_com8 }
	    if { $isosurf(Z_Sel) == "offset" }   { SetIsoGridSpace_com9 }
	} elseif { $periodic(dim) == 0 } {
	    ############
	    # MOLECULE #
	    ############
	    proc SetIsoGridSpace_com10 {} {
		foreach ent {.isoget.f2.f1.1.e1 .isoget.f2.f1.2.e2 \
			.isoget.f2.f1.3.e3} { 
		    $ent configure -state normal -relief sunken 
		}
	    }
	    
	    set r21 [radiobutton $f21.r1 \
		    -text "Whole structure with offsets" \
		    -variable isosurf(space_sel) \
		    -value "whole_cell" \
		    -width 41 \
		    -anchor w \
		    -command SetIsoGridSpace_com10]
	    
	    pack $r21  -side top -fill x -expand 1 -padx 5
	    
	    frame $f21.1 -relief flat
	    frame $f21.2 -relief flat
	    frame $f21.3 -relief flat
	    label $f21.1.l1 -text "          X-offset: " -anchor w -width 18
	    label $f21.2.l2 -text "          Y-offset: " -anchor w -width 18
	    label $f21.3.l3 -text "          Z-offset: " -anchor w -width 18
	    entry $f21.1.e1 -textvariable isosurf(X_offset) -width 20
	    entry $f21.2.e2 -textvariable isosurf(Y_offset) -width 20
	    entry $f21.3.e3 -textvariable isosurf(Z_offset) -width 20
	    append varlist { {isosurf(X_offset) real} \
		    {isosurf(Y_offset) real} \
		    {isosurf(Z_offset) real} }
	    lappend foclist $f21.1.e1 $f21.2.e2 $f21.3.e3
	    pack $f21.1 $f21.2 $f21.3 -side top -padx 5 -fill x
	    pack $f21.1.l1 $f21.2.l2 $f21.3.l3 -side left
	    pack $f21.1.e1 $f21.2.e2 $f21.3.e3 -side left	
	    
	    IsoCellSel
	    ##########
	    # get init states
	    if { $isosurf(space_sel) == "whole_cell" } { 
		SetIsoGridSpace_com10 
	    }
	}
	
	proc SetIsoGridSpace_com11 {} {
	    xcDisableAll {.isoget.f2.f1.mar \
		    .isoget.f2.f1.1 .isoget.f2.f1.2 .isoget.f2.f1.3}
	    
	    set p0 .isoget.f2.f1.mar
	    set p1 .isoget.f2.f1.1
	    set p2 .isoget.f2.f1.2
	    set p3 .isoget.f2.f1.3
	    foreach ent { e1 e2 e3 e4 e5 e6 } { 
		catch {$p0.$ent configure -relief flat }
		catch {$p1.$ent configure -relief flat }
		catch {$p2.$ent configure -relief flat }
		catch {$p3.$ent configure -relief flat }
	    }
	}
	
	set r22 [radiobutton $f21.r2 \
		-text "Manually selected space" \
		-variable isosurf(space_sel) \
		-value "manually" \
		-width 41 \
		-command SetIsoGridSpace_com11 \
		-anchor w]
	### t.k.: temporarily
	$r22 config -state disabled
	###
	pack $r22  -side top -fill x -expand 1 -padx 5
    } elseif { $dim == "2D" } {
	IsoPlaneGrid $f2
    }

    # warning: $ok should be equal to isoplane(OK); look 5 lines above !!!
    set ok  [button $isoplane(OK) -text "Submit" \
		 -command [list SetIsoGridSpaceOK $t $comlist $dim]]
    if { [xcIsActive properties] } { 
	set edit [button $isoplane(editButton) -text "Edit Script & Submit" \
		      -command [list SetIsoGridSpaceOK $t $comlist $dim edit]]
    }

    if { $isoplane(OK) == "2D" } {
	$ok   config -state disabled
	if { [info exists edit] } {
	    config -state disabled
	}
    }
    proc SetIsoGridSpace_Cancel {t dim} {
	global isoplane

	xc_isospacesel .mesa clear

	if [xcIsActive wien] {
	    # if cancel button was pressed warn the user that 
	    # application will exit
	    exit_pr
	} else {
	    CancelProc $t
	}
    }
    set can [button $f3.can -text "Cancel" \
	    -command [list SetIsoGridSpace_Cancel $t $dim]]
    
    pack $can  -side left -expand 1 -padx 10 -pady 10
    if { [info exists edit] } {
	pack $edit -side left -expand 1 -padx 10 -pady 10
    }
    pack $ok   -side left -expand 1 -padx 10 -pady 10

    ##########
    # get init states
    if { $isosurf(space_sel) == "manually" } { SetIsoGridSpace_com11 }
}


proc SetIsoGridSpaceOK {t comlist {dim 3D} {edit {}}} {
    global varlist foclist isosurf periodic prop radio vec system \
	    dif_isosurf isostack isosign isodata isofiles mody \
	    spin_entry 	isoplane
    # radio(cellmode).......which mode for crystal cell; prim/conv

    set f0 $t.f0
    set f2 $t.f2

    # if UHF, check that spin was specified
    if { [winfo exists $f0] } {
	if { ! [info exists isosurf(spin)] } {
	    set isosurf(spin) ""
	}
	if { $isosurf(spin) == "" } {
	    dialog .number2 ERROR "ERROR !\nYou forgot to specify the spin" error 0 OK
	    return
	}
    }

    if { $edit == {} } {
	set prop(editScript) 0
    } else {
	set prop(editScript) 1
    }

    # check if spin was entered correctly
    if { $prop(spin_case) } {
	if { ! [CheckSpin $spin_entry] } {
	    return
	}
    }

    set ie {}
    if { $isosurf(res_type) == "points" } { 
	lappend ie 1 
    } else { 
	lappend ie 2
    }

    if { $dim == "3D" } {
	if { $isosurf(space_sel) == "whole_cell" } {
	    if { $periodic(dim) == 0 } {
		append ie { 3 4 5}
	    }
	    if { $periodic(dim) == 1 } { 
		append ie { 3 6}
		if { $isosurf(Y_Sel) == "centered" } { 
		    lappend ie 9 
		} elseif { $isosurf(Y_sel) == "min_max" } {
		    append ie { 10 11}
		} else {
		    lappend ie 12
		}
		if { $isosurf(Z_Sel) == "centered" } { 
		    lappend ie 13 
		} elseif { $isosurf(Z_Sel) == "min_max" } {
		    append ie { 14 15}
		} else {
		    lappend ie 16
		}
	    }    
	    if { $periodic(dim) == 2 } { 
		append ie { 3 4  6 7}
		if { $isosurf(Z_Sel) == "centered" } {
		    lappend ie 9
		} elseif { $isosurf(Z_Sel) == "min_max" } {
		    append ie { 10 11}
		} else {
		    lappend ie 12
		}
	    }
	}
        
    } elseif { $dim == "2D" } {
	#
	# was the isoplane correctly selected
	if { $isoplane(plane_sel) == {} } { return }
	#xc_isospacesel .mesa clear
    }

    foreach i $ie {
	incr i -1	
	lappend new_varlist [lindex $varlist $i]
	lappend new_foclist [lindex $foclist $i]
    }
        
    if { ! [check_var $new_varlist $new_foclist] } {
	return
    }

    if { $dim == "2D" && ! [info exists isosurf(origVec)] } {
	# the "Update-Display" button was not pressed for 2D; update a display
	IsoPlaneGrid_Update $f2
    }    
	
    if { $comlist == "UPDATE_DISPLAY" } {
	# we should just update display and return
	IsoSpaceUpdate $periodic(dim)
	set isosurf(origVec) [xc_isospacesel .mesa cell3D \
		-ctype $radio(cellmode) \
		-margins [list $isosurf(A_margin) $isosurf(B_margin) \
		$isosurf(C_margin) $isosurf(A*_margin) \
		$isosurf(B*_margin) $isosurf(C*_margin)]]
	return 0
    }

    if { [winfo exists $t] } { destroy $t }

    ###################################################################
    # OK now everything need to be prepared for isosurface evaluation #
    ###################################################################
    # initialise origin of "isosurf cell"
    for {set i 0} {$i < 3} {incr i} {
	set isosurf(origin,$i) 0.0
    }
    	    

    ######################
    # for difference-maps we have two command in $comlist, but 
    # NPX,NPY & NPZ must be equal for both, so take the first command 
    # and determine NPX, NPZ, NPY
    set command [lindex $comlist 0]

    if { $dim == "3D" } {
	# t.k::
	# next two lines are not tested for CRYSTALXX (if & xc_isosp...)
	if { ! [info exists isosurf(origVec)] } {
	    IsoSpaceUpdate $periodic(dim)
	    set isosurf(origVec) [xc_isospacesel .mesa cell3D \
		    -ctype $radio(cellmode) \
		    -margins [list $isosurf(A_margin) $isosurf(B_margin) \
		    $isosurf(C_margin) $isosurf(A*_margin) \
		    $isosurf(B*_margin) $isosurf(C*_margin)]]
	}
	xc_isospacesel .mesa clear
	if { $comlist == "WIEN" } {
	    global wn
	    # --wien_density option;
	    # "whole_cell" & "manual" option !!!!
	    wnMakeIn5_2D3D $isosurf(origVec) 3D
	    set wn(done) 1
	    return
	} else {
	    # 
	    # this is for CRYSTALXX
	    #

	    # delete old unit-25
	    if { [file exists $prop(dir)/$prop(file)25] } {
		file delete $prop(dir)/$prop(file)25
	    }
	    
	    if { $isosurf(space_sel) == "whole_cell" } {
		GetCageVecOrig $isosurf(origVec)
		GetNPY $command
		GetNumberOfPoints $command 3
	    }
	}
    } elseif { $dim == "2D" } {
	xc_isospacesel .mesa clear

	# order of points in isoplane(points)::
	# 0--3
	# |  |
	# 1--2; p1 will be set to origin, vec0 = p2 - p1, vec1 = p0 - p1 
	
	#
	# CRYSTAL or WIEN ????
	#
	if { $comlist == "WIEN" } {
	    global wn
	    # --wien_density option;
	    wnMakeIn5_2D3D $isoplane(points) 2D
	    set wn(done) 1
	    return
	} else {
	    # CRYSTALXX
	    set p $isoplane(points)
	    xcDebug "\npoints:: $isoplane(points)\n"
	    set isosurf(origin,0) [Angs2Bohr [lindex $p 3]]
	    set isosurf(origin,1) [Angs2Bohr [lindex $p 4]]
	    set isosurf(origin,2) [Angs2Bohr [lindex $p 5]]
	    
	    set vec(0,0) [Angs2Bohr [expr [lindex $p 6] - [lindex $p 3]]]
	    set vec(0,1) [Angs2Bohr [expr [lindex $p 7] - [lindex $p 4]]]
	    set vec(0,2) [Angs2Bohr [expr [lindex $p 8] - [lindex $p 5]]]
	    
	    set vec(1,0) [Angs2Bohr [expr [lindex $p 0] - [lindex $p 3]]]
	    set vec(1,1) [Angs2Bohr [expr [lindex $p 1] - [lindex $p 4]]]
	    set vec(1,2) [Angs2Bohr [expr [lindex $p 2] - [lindex $p 5]]]
	    
	    set vec(2,0) 0.0
	    set vec(2,1) 0.0
	    set vec(2,2) 0.0
	    
	    GetNPY $command
	    set prop(NPZ) 1
	}
    }
    ########################################
    # now caclulate isosurfaces
    set n_com [llength $comlist]
    set com1  [lindex $comlist 0]
    set com2  [lindex $comlist 1]
    xcDebug "\nCOMLIST>> $com1 $com2"

    ############# CALCULATE map_A
    update
    update idletask
    if { ! [IsoCalc 1 $com1] } {
	return
    }
    update
    #############
    set prop(unit25_1) $system(SCRDIR)/$prop(file)25.A    
    file rename -force $system(SCRDIR)/$prop(file)25 $prop(unit25_1)

    # If diff. map calculate map B
    # For diff. map we may use one unit9 or load each map from separate unit9
    # For diff. maps we always make two units 25 
    if { $n_com > 1 } { 
	set dirB $system(SCRDIR)
	if { $dif_isosurf(load_another_unit9) == 1 } {
	    file copy -force $dif_isosurf(unit9_B) $system(SCRDIR_1)/$prop(file)9
	    set dirB $system(SCRDIR_1)
	}
	############### CALCULATE map_B
	if { ! [IsoCalc 2 $com2 $dirB/xc_output.$system(PID) $dirB] } {
	    return
	}
	###############
	set prop(unit25_2)    $system(SCRDIR)/$prop(file)25.B
	file rename -force $dirB/$prop(file)25 $prop(unit25_2)
    }

    ###################################################
    # now we have to determine xc_isostack parameters #
    ###################################################
    
    #######################
    # ISOFILES parameters #
    #######################
    set isofiles "$prop(unit25_1) "
    if { $n_com > 1 } { 
	append isofiles "$prop(unit25_2)"
    }

    #######################
    # ISOSTACK parameters #
    #######################
    set isostack {}
    if { $dif_isosurf(dif_map) == 0 } {
	#######################
	# NOT DIFFERENCE MAPS #
	#######################
	append isostack "xc_isostack 3; "
	append isostack "xc_isostack 0 2 1; "
	append isostack "xc_isostack 0 1 $prop(NPZ); "
    } elseif { $dif_isosurf(dif_map) == 1 } {
	###################
	# DIFFERENCE MAPS #
	###################
	append isostack "xc_isostack 4; "
	append isostack "xc_isostack 0 3 2; "
	append isostack "xc_isostack 0 2 1; "
	append isostack "xc_isostack 1 2 1; "
	append isostack "xc_isostack 0 1 $prop(NPZ); "
	append isostack "xc_isostack 1 1 $prop(NPZ); "
    }

    #########################
    # set isosign & isodata #
    #########################
    SetXC_Iso $dim

    ########################################
    #   PREPARE ISOSURFACE/PROPERTYPLANE   #
    ########################################
    xcPrepareIsosurf $com1 $dim

    #################################################################
    # now make isosurf_struct global variable, which will be used for
    # isoControl tracking system of "isoControl" changes
    Set_UpdateIsosurf_Struct

    #######################
    # display IsoControls #
    if { $dim == "3D" } {
	set prop(datagridDim) 3
	IsoControl 
    } elseif { $dim == "2D" } {
	set prop(datagridDim) 2
	IsoControl2D
    }
}



proc CheckSpin {foc} {
    global isosurf 

    set s_list "\"ALPHA\" \"BETA\" \"ALPHA-BETA\" \"ALPHA+BETA\""
    if ![string match "*$isosurf(spin)*" $s_list] {
	tk_dialog [WidgetName] ERROR "ERROR: unknown spin \"$isosurf(spin)\", should be of $s_list" error 0 OK
	focus $foc
	return 0
    }
    return 1
}



proc IsoCellSel {} {
    uplevel 1 {
	frame $f21.mar -relief flat
	pack $f21.mar -expand 1
	
	array set isosurf {A_margin 0.0 B_margin 0.0 C_margin 0.0 \
		A*_margin 0.0 B*_margin 0.0 C*_margin 0.0}
	
	label $f21.mar.l1 -text " A margin: " -anchor e -width 12
	label $f21.mar.l2 -text " B margin: " -anchor e -width 12
	label $f21.mar.l3 -text " C margin: " -anchor e -width 12
	label $f21.mar.l4 -text "A* margin: " -anchor e -width 12
	label $f21.mar.l5 -text "B* margin: " -anchor e -width 12
	label $f21.mar.l6 -text "C* margin: " -anchor e -width 12
	
	entry $f21.mar.e1 -textvariable isosurf(A_margin) -width 7
	entry $f21.mar.e2 -textvariable isosurf(B_margin) -width 7
	entry $f21.mar.e3 -textvariable isosurf(C_margin) -width 7
	entry $f21.mar.e4 -textvariable isosurf(A*_margin) -width 7
	entry $f21.mar.e5 -textvariable isosurf(B*_margin) -width 7
	entry $f21.mar.e6 -textvariable isosurf(C*_margin) -width 7
	append varlist { {isosurf(A_margin) real} \
		{isosurf(B_margin) real} \
		{isosurf(C_margin) real} \
		{isosurf(A*_margin) real} \
		{isosurf(B*_margin) real} \
		{isosurf(C*_margin) real} }
	lappend foclist $f21.mar.e1 $f21.mar.e2 $f21.mar.e3 $f21.mar.e4 $f21.mar.e5 $f21.mar.e6
	
	set list {}
	if { $periodic(dim) == 1 } {
	    set list {1 4}
	} elseif { $periodic(dim) == 2 } {
	    set list {1 2  4 5}
	} elseif { $periodic(dim) == 3 } {
	    set list {1 2 3 4 5 6}
	}

	set j 0
	foreach i $list {
	    set ii [expr ($i-1) * 2 - $j*6]
	    grid $f21.mar.l$i -row $j -column $ii -pady 3
	    grid $f21.mar.e$i -row $j -column [expr $ii + 1] -pady 3
	    if { $i == $periodic(dim) } { incr j }
	}
	
	frame $f21.mar.2 -relief raised -bd 1
	set upd [button $f21.mar.2.upd -text "Update Display" -bd 1 \
		-command [list SetIsoGridSpaceOK $t UPDATE_DISPLAY 3D]]
	if { $periodic(dim) < 1 } {
	    grid $f21.mar.2 -row 2 -column 0 \
		    -columnspan 2 -sticky ew -pady 5
	} else {
	    grid $f21.mar.2 -row 2 -column [expr $periodic(dim) - 1] \
		    -columnspan 2 -sticky ew -pady 5
	}
	pack $upd -pady 3
    }
}



proc IsoSpaceUpdate dim {
    global mody isosurf

    if { $dim < 3 && $dim > 0 } {
	if { $isosurf(Z_Sel) == "centered" } {
	    # isosurf(Z_size)
	    set z2 [expr $isosurf(Z_size) / 2.0]
	    set isosurf(C_margin)  $z2
	    set isosurf(C*_margin) $z2
	} elseif { $isosurf(Z_Sel) == "min_max" } {
	    set isosurf(C_margin)  [expr -1.0 * $isosurf(Z_min)]
	    set isosurf(C*_margin) $isosurf(Z_max)
	} elseif { $isosurf(Z_Sel) == "offset" } {
	    set isosurf(C_margin) [expr -1.0 * \
		    [xc_getvalue $mody(GET_SS_MINZ)] + $isosurf(Z_offset)]
	    set isosurf(C*_margin) [expr \
		    [xc_getvalue $mody(GET_SS_MAXZ)] + $isosurf(Z_offset)]
	}
    }

    if { $dim == 1 } {
	if { $isosurf(Y_Sel) == "centered" } {
	    # isosurf(Y_size)
	    set y2 [expr $isosurf(Y_size) / 2.0]
	    set isosurf(B_margin)  $y2
	    set isosurf(B*_margin) $y2
	} elseif { $isosurf(Y_Sel) == "min_max" } {
	    set isosurf(B_margin)  [expr -1.0 * $isosurf(Y_min)]
	    set isosurf(B*_margin) $isosurf(Y_max)
	} elseif { $isosurf(Y_Sel) == "offset" } {
	    set isosurf(B_margin) [expr -1.0 * \
		    [xc_getvalue $mody(GET_SS_MINY)] + $isosurf(Y_offset)]
	    set isosurf(B*_margin) [expr \
		    [xc_getvalue $mody(GET_SS_MAXY)] + $isosurf(Y_offset)]
	}
    }

    if { $dim == 0 } {
	set isosurf(C_margin) [expr -1.0 * \
		[xc_getvalue $mody(GET_SS_MINZ)] + $isosurf(Z_offset)]
	set isosurf(C*_margin) [expr \
		[xc_getvalue $mody(GET_SS_MAXZ)] + $isosurf(Z_offset)]
	set isosurf(B_margin) [expr -1.0 * \
		[xc_getvalue $mody(GET_SS_MINY)] + $isosurf(Y_offset)]
	set isosurf(B*_margin) [expr \
		[xc_getvalue $mody(GET_SS_MAXY)] + $isosurf(Y_offset)]
	set isosurf(A_margin) [expr -1.0 * \
		[xc_getvalue $mody(GET_SS_MINX)] + $isosurf(X_offset)]
	set isosurf(A*_margin) [expr \
		[xc_getvalue $mody(GET_SS_MAXX)] + $isosurf(X_offset)]	
    }
}
