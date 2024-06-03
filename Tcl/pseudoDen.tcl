#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/pseudoDen.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc PseudoDensity {} {
    global pDen check

    # incoming ...
    # THIS ...
    #if { ! $check(pseudoDens) && ![info exists pDen(nsurface)] } {
    #	return
    #}
    #/

    if { $check(pseudoDens) } {
	if { ![info exists pDen(nsurface)] } {
	    PseudoDenDefaults 0
	    set pDen(nsurface) 1
	    set pDen(0,ident) [xc_molsurfreg .mesa]	    
	    
	    SetWatchCursor
	    update
	    xc_molsurf .mesa \
		-ident        $pDen(0,ident) \
		-type         $pDen(0,type) \
		-radius       $pDen(0,radius) \
		-level        $pDen(0,level) \
		-cutoff       $pDen(0,cutoff) \
		-colorscheme  $pDen(0,colorscheme) \
		-drawstyle    $pDen(0,drawstyle) \
		-transparent  $pDen(0,transparent) \
		-shademodel   $pDen(0,shademodel) \
		-monocolor    $pDen(0,monocolor) \
		-surfacetype  $pDen(0,surfacetype) \
		-resolution   $pDen(0,resolution) \
		-smoothsteps  $pDen(0,smoothsteps) \
		-smoothweight $pDen(0,smoothweight) \
		-tessellation $pDen(0,tessellation) \
		-normals      $pDen(0,normals)
	    set pDen(0,new) 0
	    ResetCursor
	    return;
	}
    }

    if { ![info exists pDen(nsurface)] } { return }

    for {set i 0} {$i < $pDen(nsurface)} {incr i} {
	xcDebug -debug "xc_molsurfconfig: $check(pseudoDens)"
	xc_molsurfconfig .mesa \
	    -ident  $pDen($i,ident) \
	    -render $check(pseudoDens)
	
    }
}


proc PseudoDenDefaults ind {
    global pDen

    set pDen($ind,render)       1
    set pDen($ind,type)         $pDen(type)         
    set pDen($ind,radius)       $pDen(radius)       
    set pDen($ind,level)        [expr $pDen(level) - $ind * 0.2]
    set pDen($ind,cutoff)       $pDen(cutoff)       
    set pDen($ind,colorscheme)  $pDen(colorscheme)  
    set pDen($ind,drawstyle)    $pDen(drawstyle)    
    set pDen($ind,surfacetype)  $pDen(surfacetype)  
    set pDen($ind,resolution)   $pDen(resolution)   
    set pDen($ind,smoothsteps)  $pDen(smoothsteps)  
    set pDen($ind,smoothweight) $pDen(smoothweight) 
    set pDen($ind,transparent)  $pDen(transparent)
    set pDen($ind,shademodel)   $pDen(shademodel)
    set pDen($ind,tessellation) $pDen(tessellation)
    set pDen($ind,normals)      $pDen(normals)

    set pDen($ind,type_old)        $pDen($ind,type)
    set pDen($ind,radius_old)      $pDen($ind,radius)
    set pDen($ind,level_old)       $pDen($ind,level)  
    set pDen($ind,cutoff_old)      $pDen($ind,cutoff) 
    set pDen($ind,surfacetype_old) $pDen($ind,surfacetype)
    set pDen($ind,resolution_old)  $pDen($ind,resolution) 	

    if { $ind == 0 } {
	set pDen($ind,monocolor)    $pDen(monocolor)
    } else {
	set i [expr $ind - 1]
	set r [expr [lindex $pDen($i,monocolor) 0] + 0.2]
	if { $r > 1.0 } { set r [expr $r - 1] }

	set g [expr [lindex $pDen($i,monocolor) 1] + 0.2]
	if { $g > 1.0 } { set g [expr $g - 1] }

	set b [expr [lindex $pDen($i,monocolor) 2] + 0.2]
	if { $b > 1.0 } { set b [expr $b - 1] }

	trace vdelete pDen($ind,monocolor) w xcTrace	
	set pDen($ind,monocolor) [list $r $g $b]
	trace variable pDen($ind,monocolor) w xcTrace
    }

    set pDen($ind,t_type)         $pDen(t_type)         
    set pDen($ind,t_radius)       $pDen(t_radius)       
    set pDen($ind,t_colorscheme)  $pDen(t_colorscheme)  
    set pDen($ind,t_drawstyle)    $pDen(t_drawstyle)    
    set pDen($ind,t_surfacetype)  $pDen(t_surfacetype)  
    set pDen($ind,t_shademodel)   $pDen(t_shademodel)
    set pDen($ind,t_tessellation) $pDen(t_tessellation)
    set pDen($ind,t_normals)      $pDen(t_normals)

    set pDen($ind,new) 1
}

proc PseudoDenGenNew t {
    global pDen

    set i $pDen(nsurface)
    PseudoDenDefaults $i
    incr pDen(nsurface)

    $pDen(nb) insert $i n$i -text "Surf #$i"
    set f [$pDen(nb) getframe n$i]
    PseudoDenNoteBook $i $f $t
    $pDen(nb) raise n$i
    set pDen($i,ident) [xc_molsurfreg .mesa]	    
    
    #SetWatchCursor
    #update
    #xc_molsurf .mesa \
    #	     -ident        $pDen($i,ident) \
    #	     -type         $pDen($i,type) \
    #	     -radius       $pDen($i,radius) \
    #	     -level        $pDen($i,level) \
    #	     -cutoff       $pDen($i,cutoff) \
    #	     -colorscheme  $pDen($i,colorscheme) \
    #	     -transparent  $pDen($i,transparent) \
    #	     -shademodel   $pDen($i,shademodel) \
    #	     -monocolor    $pDen($i,monocolor) \
    #	     -drawstyle    $pDen($i,drawstyle) \
    #	     -surfacetype  $pDen($i,surfacetype) \
    #	     -resolution   $pDen($i,resolution) \
    #	     -smoothsteps  $pDen($i,smoothsteps) \
    #	     -smoothweight $pDen($i,smoothweight)
    #ResetCursor
}

proc PseudoDenConfigWid {} {
    global pDen check

    set t [xcToplevel [WidgetName] "Molecular Surface Settings\
	    " "Molecular Surface" . 0 0 1]

    button $t.b1 -text "Generate new molecular-surface" \
	    -command [list PseudoDenGenNew $t]
    pack $t.b1 -side top -expand 1 -fill x -padx 5m -pady 2m

    set nb $t.n
    set pDen(nb) $nb
    #tixNoteBook $nb
    # now we use BWidgets
    NoteBook $nb
    pack $nb

    set ind 0
    if { ![info exists pDen(nsurface)] } { 
	set pDen(nsurface) 0
	PseudoDenGenNew $t
	set ind 1
    }
    for {set i $ind} {$i < $pDen(nsurface)} {incr i} {
	$nb insert $i n$i -text "Surf #$i"
	set f [$nb getframe n$i]
	PseudoDenNoteBook $i $f $t	
	$nb raise n$i
    }
}

proc PseudoDenNoteBook {ind f t} {
    global pDen
    #		 -ident        <identifier>
    #		 -type         gauss|exp|unigauss|uniexp 
    #		 -radius       cov|VdW 
    #		 -level        <isoleve> 
    #		     -cutoff       <cutoff> 
    #		 -colorscheme  atomic|monochrome 
    #		 -drawstyle    solid|wire|dot 
    #		     -surfacetype  molsurf|gap
    #		     -resolution   <x> 
    #		     -smoothsteps  <steps> 
    #		     -smoothweight <weight>
    

    set f1 [frame $f.1]
    set f2 [frame $f.2 -relief ridge -bd 2]
    pack $f1 $f2 -side top -padx 2m -pady 2m -fill x

    checkbutton $f1.cb0 -text "Display surface:" -variable pDen($ind,render) \
	    -command [list PseudoDenUpdateSurf $ind]
    
    pack $f1.cb0 -side top -pady 2

    set mb1 [xcMenuButton $f1 \
	    -labeltext "Function type:" -labelwidth 20 \
	    -textvariable pDen($ind,t_type) \
	    -menu [list \
	    GAUSSIAN               [list set pDen($ind,t_type) GAUSSIAN] \
	    EXPONENTIAL            [list set pDen($ind,t_type) EXPONENTIAL] \
	    {constant GAUSSIAN}    [list set pDen($ind,t_type) {constant GAUSSIAN}] \
	    {constant EXPONENTIAL} [list set pDen($ind,t_type) {constant EXPONENTIAL}] \
	    {distance FUNCTION}    [list set pDen($ind,t_type) {distance FUNCTION}]] ]

    set mb2 [xcMenuButton $f1 \
	    -labeltext "Radius type:" -labelwidth 20 \
	    -textvariable pDen($ind,t_radius) \
	    -menu [list \
	    {Covalent radii}      [list set pDen($ind,t_radius) {Covalent radii}] \
	    {Van der Waals radii} [list set pDen($ind,t_radius) {Van der Waals radii}]]
    ] 

    set mb3 [xcMenuButton $f1 \
	    -labeltext "Color Scheme:" -labelwidth 20 \
	    -textvariable pDen($ind,t_colorscheme) \
	    -menu [list \
	    {Atomic colors}      [list set pDen($ind,t_colorscheme) {Atomic colors}] \
	    {Monocolor}          [list set pDen($ind,t_colorscheme) {Monocolor}]]
    ]
    
    set mb4 [xcMenuButton $f1 \
	    -labeltext "Draw Style:" -labelwidth 20 \
	    -textvariable pDen($ind,t_drawstyle) \
	    -menu [list \
	    Solid      [list set pDen($ind,t_drawstyle) Solid] \
	    Wire       [list set pDen($ind,t_drawstyle) Wire] \
	    Dot        [list set pDen($ind,t_drawstyle) Dot]]
    ]
    
    set mb5 [xcMenuButton $f1 \
	    -labeltext "Shade Model:" -labelwidth 20 \
	    -textvariable pDen($ind,t_shademodel) \
	    -menu [list \
	    {Smooth}      [list set pDen($ind,t_shademodel) {Smooth}] \
	    {Flat}        [list set pDen($ind,t_shademodel) {Flat}]]
    ]

    set mb51 [xcMenuButton $f1 \
		  -labeltext "Surface tessellation type:" -labelwidth 20 \
		  -textvariable pDen($ind,t_tessellation) \
		  -menu [list \
			     {Cube}        [list set pDen($ind,t_tessellation) {Cube}] \
			     {Tetrahedron} [list set pDen($ind,t_tessellation) {Tetrahedron}]]
	     ]
    
    set mb52 [xcMenuButton $f1 \
		  -labeltext "Surface normals type:" -labelwidth 20 \
		  -textvariable pDen($ind,t_normals) \
		  -menu [list \
			     {Gradient}  [list set pDen($ind,t_normals) {Gradient}] \
			     {Triangle}  [list set pDen($ind,t_normals) {Triangle}]]
	     ]
    
    set mb6 [xcMenuButton $f1 \
		 -labeltext "Surface Type:" -labelwidth 20 \
		 -textvariable pDen($ind,t_surfacetype) \
		 -menu [list \
			    {Pseudo density}      [list set pDen($ind,t_surfacetype) {Pseudo density}] \
			    {Gap analysis}        [list set pDen($ind,t_surfacetype) {Gap analysis}]]
	    ]
    
    ###########################################################
    # TEMPORAL: until the C-code for GAP-ANALYSIS for slabs-polymers-molecules
    #          is written
    global periodic
    if { $periodic(dim) < 3 } {
	$mb6.mb.menu entryconfig {Gap analysis} -state disabled
    }
    ###########################################################

    pack $mb1 $mb2 $mb3 $mb4 $mb5 $mb51 $mb52 $mb6 -side top -expand 1 -fill x -padx 2m -pady 2

    global fillEntries
    set e [FillEntries $f1 [list \
	    "Isovalue (between \[0.0,2.0\]):" "Cutoff:" \
	    "Resolution (in Angs.):" "N. of smoothing steps:" \
	    "Smoothing weight:"] \
	    [list pDen($ind,level) pDen($ind,cutoff) pDen($ind,resolution) \
	    pDen($ind,smoothsteps) pDen($ind,smoothweight)] \
	    30 7]
    set foclist $fillEntries
    set varlist [list pDen($ind,level) \
	    pDen($ind,cutoff) pDen($ind,resolution) \
	    pDen($ind,smoothsteps) pDen($ind,smoothweight)]

    set f11 [frame [WidgetName $f1]]
    set b  [button $f11.b -text "Set Surface Monocolor" \
	    -command [list PseudoDenChangeCol $ind]] 
    set l  [label $f11.l -text "Current monocolor:"]
    set pDen($ind,colbut) [button $f11.colb \
	    -bg [rgb_f2h $pDen($ind,monocolor)] -width 5 \
            -state disabled -relief sunken -bd 1]
    trace variable pDen($ind,monocolor) w xcTrace

    pack $f11 -side top -pady 2m -fill x
    pack $b $l $pDen($ind,colbut) -side left -pady 2 -fill x

    checkbutton $f1.cb -text "Surface transparency:" \
	    -variable pDen($ind,transparent)
    pack $f1.cb -side top -pady 2

    button $f2.1 -text "Close"  -command [list CancelProc $t]
    button $f2.2 -text "Update" \
	    -command [list PseudoDenUpdate $ind $varlist $foclist]
    pack $f2.1 $f2.2 -side left -expand 1 -padx 2m -pady 3
}


proc PseudoDenUpdate {ind varlist foclist} {
    global pDen

    set all(type) {
	{GAUSSIAN gauss}
	{EXPONENTIAL exp}
	{{constant GAUSSIAN} unigauss}
	{{constant EXPONENTIAL} uniexp}
	{{distance FUNCTION} distf}
    }
    set all(radius) {
	{{Covalent radii} cov}
	{{Van der Waals radii} VdW}
    }
    set all(colorscheme) {
	{{Atomic colors} atomic}
	{Monocolor monochrome}
	{{Distance colors} distance}
	{{Slab colors} slab}
    }
    set all(drawstyle) {
	{Solid solid}
	{Wire wire}
	{Dot dot}
    }
    set all(shademodel) {
	{Smooth smooth}
	{Flat   flat}
    }
    set all(surfacetype) {
	{{Pseudo density} molsurf}
	{{Gap analysis} gap}
    }

    set all(tessellation) {
	{Cube cube}
	{Tetrahedron tetrahedron}
    }
    set all(normals) {
	{Gradient gradient}
	{Triangle triangle}
    }
    
    foreach var {type radius colorscheme drawstyle shademodel surfacetype tessellation normals} {
	set match 0
	foreach item $all($var) {
	    set text  [lindex $item 0]
	    set value [lindex $item 1]
	    if { $text == $pDen($ind,t_$var) } {
		set pDen($ind,$var) $value
		set match 1
	    }
	}
	if { !$match } {
	    tk_dialog [WidgetName] Error "PseudoDenUpdate Error" error 0 OK
	}
    }

    if { ![check_var $varlist $foclist] } {
	return
    }
    PseudoDenUpdateSurf $ind
}


proc PseudoDen:_xc_molsurf {ind} {
    global pDen
    xc_molsurf .mesa \
	-ident        $pDen($ind,ident) \
	-type         $pDen($ind,type) \
	-radius       $pDen($ind,radius) \
	-level        $pDen($ind,level) \
	-cutoff       $pDen($ind,cutoff) \
	-colorscheme  $pDen($ind,colorscheme) \
	-transparent  $pDen($ind,transparent) \
	-shademodel   $pDen($ind,shademodel) \
	-monocolor    $pDen($ind,monocolor) \
	-drawstyle    $pDen($ind,drawstyle) \
	-surfacetype  $pDen($ind,surfacetype) \
	-resolution   $pDen($ind,resolution) \
	-smoothsteps  $pDen($ind,smoothsteps) \
	-smoothweight $pDen($ind,smoothweight) \
	-tesselation  $pDen($ind,tessellation) \
	-normals      $pDen($ind,normals)
}


proc PseudoDenUpdateSurf ind {
    global pDen check

    if {  ( $pDen($ind,type_old)        != $pDen($ind,type) || \
	    $pDen($ind,radius_old)      != $pDen($ind,radius) || \
	    $pDen($ind,cutoff_old)      != $pDen($ind,cutoff) || \
	    $pDen($ind,surfacetype_old) != $pDen($ind,surfacetype) || \
	    $pDen($ind,resolution_old)  != $pDen($ind,resolution) ) \
	    || \
	    $pDen($ind,new) == 1 } {
	set pDen($ind,new) 0
	PseudoDen:_xc_molsurf $ind
    } else {
	xc_molsurfconfig .mesa \
	    -ident        $pDen($ind,ident) \
	    -render       $pDen($ind,render) \
	    -level        $pDen($ind,level) \
	    -colorscheme  $pDen($ind,colorscheme) \
	    -transparent  $pDen($ind,transparent) \
	    -shademodel   $pDen($ind,shademodel) \
	    -monocolor    $pDen($ind,monocolor) \
	    -drawstyle    $pDen($ind,drawstyle) \
	    -smoothsteps  $pDen($ind,smoothsteps) \
	    -smoothweight $pDen($ind,smoothweight) \
	    -tessellation $pDen($ind,tessellation) \
	    -normals      $pDen($ind,normals)	      
    }
    
    set pDen($ind,type_old)        $pDen($ind,type)
    set pDen($ind,radius_old)      $pDen($ind,radius)
    set pDen($ind,level_old)       $pDen($ind,level)  
    set pDen($ind,cutoff_old)      $pDen($ind,cutoff) 
    set pDen($ind,surfacetype_old) $pDen($ind,surfacetype)
    set pDen($ind,resolution_old)  $pDen($ind,resolution) 	

    set check(pseudoDens) 1
}

proc PseudoDenChangeCol ind {
    global pDen

    set c1 [d2h [expr [lindex $pDen($ind,monocolor) 0] * 255]]
    set c2 [d2h [expr [lindex $pDen($ind,monocolor) 1] * 255]]
    set c3 [d2h [expr [lindex $pDen($ind,monocolor) 2] * 255]]

    set t [xcToplevel [WidgetName] "Set Canvas Background" "SetBg" . 0 0 1]
    xcModifyColor $t "Set Surface Monocolor" #${c1}${c2}${c3} groove \
	    left left 100 100 70 5 20

    proc PseudoDenChangeColOK {type t ind} {
	global mody_col mody pDen
	
	if { $type == "OK" } {
	    set cID [xcModifyColorGetID]
	    set pDen($ind,monocolor) [list $mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue) 1.0]
	    set pDen($ind,colbut) [rgb_f2h [list $mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue)]]
	}
	destroy $t
    }

    set ok  [DefaultButton [WidgetName $t] -text "OK" \
	    -command [list PseudoDenChangeColOK OK     $t $ind]]
    set can [button [WidgetName $t] -text "Cancel" \
	    -command [list PseudoDenChangeColOK Cancel $t $ind]]
    pack $ok $can -padx 10 -pady 10 -expand 1
}
