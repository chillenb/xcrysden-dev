#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/colors.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc ColorSchemeState {f1 f2 f3} {
    global colSh
    xcEnableAll  $f1
    xcDisableAll $f2
    xcDisableAll $f3
    xcEnableOne [list $colSh(r1) $colSh(r2) $colSh(r3)]
}


proc ColorSchemeUpdate w {
    global colSh

    puts stderr "------------------------------> ColorSchemeUpdate"
    
    switch -- $colSh(scheme) {
	atomic { xc_colorscheme structure $w -default }
	slab   { 
	    if { $colSh(slab_fractional) } {
		set min $colSh(slabrange_min)
		set max $colSh(slabrange_max)
	    } else {
		global mody
		if { [string match ?x $colSh(slab_dir)] } {
		    set size2 [xc_getvalue $mody(GET_AT_MAXX)]
		} elseif { [string match ?y $colSh(slab_dir)] } {
		    set size2 [xc_getvalue $mody(GET_AT_MAXY)]
		} else {
		    set size2 [xc_getvalue $mody(GET_AT_MAXZ)]
		}
		
		set min [expr {0.5 * ($colSh(slabrange_min) + $size2) / $size2}]
		set max [expr {0.5 * ($colSh(slabrange_max) + $size2) / $size2}]
		
	    }
	    puts stderr "------------------------------> min,max = $min, $max"
	    xc_colorscheme structure $w -slab $colSh(slab_dir) \
		-slabrange   [list $min $max] \
		-colorscheme $colSh(slab_colbas) \
		-colortype   $colSh(slab_coltyp) -alpha $colSh(slab_alpha) 
	}
	dist   { xc_colorscheme structure $w \
		     -nn [list $colSh(dist_x) $colSh(dist_y) $colSh(dist_z)] \
		     -colorscheme $colSh(dist_colbas) \
		     -colortype $colSh(dist_coltyp) -alpha $colSh(dist_alpha) \
		     -r $colSh(dist_r)
	}
    }
}


proc ColorSchemeClose w {
    ColorSchemeUpdate $w
    ColorSchemeCancel
}


proc ColorSchemeCancel {} {
    global colSh
    set t [winfo toplevel $colSh(f)]
    CancelProc $t
    unset colSh(f)
}


proc ColorSchemeCoor {} {
    global colSh SelHoleCL

    if [info exists SelHoleCL(transl)] {unset SelHoleCL(transl)}
    PreSel .addatom .mesa "Select coordinates" \
	    "Select coordinates by clicking on desired number of atoms.\n\
	    Selected coorinates are geometrical center of selected atoms." \
	    SelCentreHoleCL 15; # 15 is maximum allowed number of 
    tkwait variable SelHoleCL(done)
    set colSh(dist_x) [lindex $SelHoleCL(centre) 0]
    set colSh(dist_y) [lindex $SelHoleCL(centre) 1]
    set colSh(dist_z) [lindex $SelHoleCL(centre) 2]
}


proc ColorScheme {} {
    global colSh mody

    if { [info exists colSh(f)] } {
	if { [winfo exists $colSh(f)] } {
	    return
	}
    }
    set f [xcUpdateWindow -cancelcom ColorSchemeCancel \
    	    -closecom "ColorSchemeClose .mesa" \
	    -updatecom "ColorSchemeUpdate .mesa"]
    set colSh(f) $f
    
    #
    # container frames
    #
    set f1 [frame $f.1 -class GrooveFrame]
    set f2 [frame $f.2 -class GrooveFrame]
    set f3 [frame $f.3 -class GrooveFrame]

    #--------------------------------------------------
    # Atomic Colors
    #--------------------------------------------------
    set colSh(r1) [radiobutton $f1.r1 \
	    -text "Colors according to atomic numbers" \
    	    -variable colSh(scheme) -value atomic -anchor w \
    	    -command [list ColorSchemeState $f1 $f2 $f3]]
    pack $f1 $f2 $f3 $colSh(r1) -side top -padx 10 -pady 10 -fill both
    
    
    #--------------------------------------------------
    # SLAB colors
    #--------------------------------------------------

    if { ! [info exists colSh(slab_min_from)] } {
	# X
	set max [xc_getvalue $mody(GET_AT_MAXX)]
	
	# Y
	set _max [xc_getvalue $mody(GET_AT_MAXY)]
	set max [expr {$_max < $max ? $_max : $max}]
	
	# Z
	set _max [xc_getvalue $mody(GET_AT_MAXZ)]
	set max [expr {$_max < $max ? $_max : $max}]
	
	set colSh(slab_min_from) [expr -1.5 * $max]
	set colSh(slab_min_to)   [expr +1.5 * $max]
	set colSh(slab_max_from) $colSh(slab_min_from)
	set colSh(slab_max_to)   $colSh(slab_min_to)
	
	set colSh(slab_absrange_min_from) $colSh(slab_min_from)
	set colSh(slab_absrange_min_to)   $colSh(slab_min_to)  
	set colSh(slab_absrange_max_from) $colSh(slab_min_from)
	set colSh(slab_absrange_max_to)   $colSh(slab_min_to)  
	
	set colSh(slab_absrange_min) $colSh(slab_absrange_min_from)
	set colSh(slab_absrange_max) $colSh(slab_absrange_max_to)
    
	if { $colSh(slab_fractional) } {
	    set colSh(slab_min_from) 0.0
	    set colSh(slab_min_to)   $colSh(slabrange_max)
	    set colSh(slab_max_from) $colSh(slabrange_min)
	    set colSh(slab_max_to)   1.0
	}
    }

    frame $f2.left
    pack $f2.left -fill both -side left
    set colSh(r2) [radiobutton $f2.r2 -text "Slab colors" \
    	    -variable colSh(scheme) -value slab -anchor w \
    	    -command [list ColorSchemeState $f2 $f1 $f3]]

    set mb1 [xcMenuButton $f2 -labeltext "Direction:" -labelwidth 12 \
	    -textvariable colSh(slab_dir) \
	    -menu {
    	+x {set colSh(slab_dir) +x}
    	-x {set colSh(slab_dir) -x}
    	+y {set colSh(slab_dir) +y}
    	-y {set colSh(slab_dir) -y}
    	+z {set colSh(slab_dir) +z}    
    	-z {set colSh(slab_dir) -z}
    }]
	
    set mb2 [xcMenuButton $f2 -labeltext "Color basis:" -labelwidth 12 \
    	     -textvariable colSh(slab_colbas) \
    	     -menu {
    	monochrome  {set colSh(slab_colbas) monochrome}
    	rgb        {set colSh(slab_colbas) rgb}
    	rainbow    {set colSh(slab_colbas) rainbow}
    	geographic {set colSh(slab_colbas) geographic}
    }]			      
    set mb3 [xcMenuButton $f2 -labeltext "Color type:" -labelwidth 12 \
    	     -textvariable colSh(slab_coltyp) \
    	     -menu {
    	override  {set colSh(slab_coltyp) override}
    	combined  {set colSh(slab_coltyp) combined}
    }]
    frame $f2.f
    FillEntries $f2.f "Alpha:" colSh(slab_alpha) 12 10 
    
    set cb [checkbutton $f2.cb -text "Specify slab-ranges in fractional units" \
		-onvalue 1 -offvalue 0 -variable colSh(slab_fractional)]
    set colSh(scMin) [scale $f2.scmin \
			  -from $colSh(slab_min_from) -to $colSh(slab_min_to) \
			  -variable colSh(slabrange_min) \
			  -orient horizontal \
			  -label "Slab-range -- minimum:" \
			  -digits 3 \
			  -resolution 0.01 \
			  -tickinterval 0.2 \
			  -showvalue true \
			  -width 10]
    set colSh(scMax) [scale $f2.scmax \
			  -from $colSh(slab_max_from) -to $colSh(slab_max_to) \
			  -variable colSh(slabrange_max) \
			  -orient horizontal \
			  -label "Slab-range -- maximum:" \
			  -digits 3 \
			  -resolution 0.01 \
			  -tickinterval 0.2 \
			  -showvalue true \
			  -width 10]
    trace variable colSh(slab_fractional) w xcTrace
    trace variable colSh(slabrange_max)   w xcTrace
    trace variable colSh(slabrange_min)   w xcTrace
    
    #frame $f2.f1
    #frame $f2.f2
    #FillEntries $f2.f1 {"Fixed minimum:"} colSh(slab_fixmin) 14 8 top top
    #FillEntries $f2.f2 {"Fixed maximum:"} colSh(slab_fixmax) 14 8 top top
    
    grid $colSh(r2) -in $f2.left -column 0 -row 0 -sticky w -padx 10 -pady 5
    grid $mb1  -in $f2.left -column 0 -row 1 -sticky w  -padx 10 -pady 5
    grid $mb2  -in $f2.left -column 1 -row 1 -sticky e  -padx 10 -pady 5
    grid $mb3  -in $f2.left -column 0 -row 2 -sticky w  -padx 10 -pady 5
    grid $f2.f -in $f2.left -column 1 -row 2 -sticky we -padx 10 -pady 5
    grid $cb   -in $f2.left -column 0 -row 3 -sticky w \
	-columnspan 2 -padx 10 -pady 5
    grid $colSh(scMin) -in $f2.left -column 0 -row 4 -sticky we \
	-columnspan 2 -padx 10 -pady 5
    grid $colSh(scMax) -in $f2.left -column 0 -row 5 -sticky we \
	-columnspan 2 -padx 10 -pady 5
    #grid $f2.f1 -in $f2.left -column 1 -row 3 -sticky we \
    #	-columnspan 1 -padx 10 -pady 5
    #grid $f2.f2 -in $f2.left -column 1 -row 4 -sticky we \
    #	-columnspan 1 -padx 10 -pady 5
    
    #--------------------------------------------------
    # distance colors
    #--------------------------------------------------
    set colSh(r3) [radiobutton $f3.r3 -text "Distance colors" \
    	    -variable colSh(scheme) -value dist -anchor w \
    	    -command [list ColorSchemeState $f3 $f1 $f2]]
    frame $f3.f1
    label $f3.f1.l -text Coordinates:
    pack $f3.f1.l -side left
    FillEntries $f3.f1 {X: Y: Z:} \
    	     {colSh(dist_x) colSh(dist_y) colSh(dist_z)} 2 8 left  left
    button $f3.f1.b -text "Select Coor" -command ColorSchemeCoor
    pack $f3.f1.b -side left
    set mb1 [xcMenuButton $f3 -labeltext "Color basis:" -labelwidth 12 \
    	     -textvariable colSh(dist_colbas) \
    	     -menu {
    	monochrome  {set colSh(dist_colbas) monochrome}
    	rgb        {set colSh(dist_colbas) rgb}
    	rainbow    {set colSh(dist_colbas) rainbow}
    	geographic {set colSh(dist_colbas) geographic}
    }]
    set mb2 [xcMenuButton $f3 -labeltext "Color type:" -labelwidth 12 \
    	     -textvariable colSh(dist_coltyp) \
    	     -menu {
    	override  {set colSh(dist_coltyp) override}
    	combined  {set colSh(dist_coltyp) combined}
    }]
    frame $f3.f2
    FillEntries $f3.f2 "Alpha:" colSh(dist_alpha) 12 10

    set colSh(scR) [scale $f3.scr \
			-from 0 -to 1 \
			-variable colSh(dist_r) \
			-orient horizontal \
			-label "Distance-range:" \
			-digits 3 \
			-resolution 0.01 \
			-tickinterval 0.2 \
			-showvalue true \
			-width 10]

    grid $colSh(r3) -column 0 -row 0 -sticky w -padx 10 -pady 5
    grid $f3.f1 -column 0 -row 1 -columnspan 2 -sticky w -padx 10 -pady 5
    grid $mb1   -column 0 -row 2 -sticky w -padx 10 -pady 5
    grid $mb2   -column 1 -row 2 -sticky e -padx 10 -pady 5
    grid $f3.f2 -column 0 -row 3 -sticky w -padx 10 -pady 5
    grid $f3.scr -column 0 -row 4 -columnspan 3 -sticky ew -padx 10 -pady 5

    switch -- $colSh(scheme) {
        atomic {ColorSchemeState $f1 $f2 $f3}
        slab   {ColorSchemeState $f2 $f1 $f3}
        dist   {ColorSchemeState $f3 $f1 $f2}
    }
}


proc ColorScheme:fromAbsToFrac {} {
    global colSh mody

    if { [string match ?x $colSh(slab_dir)] } {
	set size2 [xc_getvalue $mody(GET_AT_MAXX)]
    } elseif { [string match ?y $colSh(slab_dir)] } {
	set size2 [xc_getvalue $mody(GET_AT_MAXY)]
    } else {
	set size2 [xc_getvalue $mody(GET_AT_MAXZ)]
    }
    
    set _min [expr {0.5 * ($colSh(slabrange_min) + $size2) / $size2}]
    set _max [expr {0.5 * ($colSh(slabrange_max) + $size2) / $size2}]
    
    if { $_min < 0.0 } { set _min 0.0 }
    if { $_min > 1.0 } { set _min 1.0 }
    if { $_max < 0.0 } { set _max 0.0 }
    if { $_max > 1.0 } { set _max 1.0 }
    
    set colSh(slabrange_min) $_min
    set colSh(slabrange_max) $_max
}


proc ColorScheme:fromFracToAbs {} {
    global colSh mody

    if { [string match ?x $colSh(slab_dir)] } {
	set size2 [xc_getvalue $mody(GET_AT_MAXX)]
    } elseif { [string match ?y $colSh(slab_dir)] } {
	set size2 [xc_getvalue $mody(GET_AT_MAXY)]
    } else {
	set size2 [xc_getvalue $mody(GET_AT_MAXZ)]
    }

    set min [expr 2.0 * $colSh(slabrange_min) * $size2 - $size2]
    set max [expr 2.0 * $colSh(slabrange_max) * $size2 - $size2]

    puts stderr "min,max:: $min,$max"

    #trace vdelete colSh(slabrange_min) w xcTrace
    #trace vdelete colSh(slabrange_max) w xcTrace
    #$colSh(scMin) set $min
    #$colSh(scMax) set $max
    set colSh(slabrange_min) $min
    set colSh(slabrange_max) $max
    puts stderr "colSh(slabrange_min): $colSh(slabrange_min), $min"
    puts stderr "colSh(slabrange_max): $colSh(slabrange_max), $max"
    #trace variable colSh(slabrange_min) w xcTrace
    #trace variable colSh(slabrange_max) w xcTrace    
}
