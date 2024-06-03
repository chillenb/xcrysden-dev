#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/FS_Main.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# NOTE: XSF file and (band)XSF file must alreadu be read !!!
#
proc FS_GoFermi {{spin {}}} {
    global fs xcMisc system

    if { ! [info exists fs(counter)] } {
	set fs(counter) 0
    } else {
	incr fs(counter) 
    }

    set fs($spin,togl_w) [expr int(750 * $xcMisc(resolution_ratio1))]
    set fs($spin,togl_h) $fs($spin,togl_w) 

    xcDebug -debug "FS_GoFermi> xcMisc(resolution_ratio1) = $xcMisc(resolution_ratio1)"

    # NOTE:
    #-------
    # prevent the mass that can be done by going several times trough
    # the "Render Fermi Surface" process for WIEN
    set t .fs${spin}
    if { [winfo exists $t] } {
	return
    }
    toplevel $t
    if { $spin == "dn" } {
    	wm geometry $t +0-0
    } else {
    	wm geometry $t -0+0
    }

    global exit_viewer_win
    set exit_viewer_win $t
    bind $t <Destroy> "exit_viewer $t"

    if { ! [info exists xcMisc(titlefile)] } {
        set xcMisc(titlefile) $fs(titlefile)
    }
    if { $spin != {} } {
	wm title $t "*** XCrySDen - Fermi Surface ($spin spin): [file tail $xcMisc(titlefile)]"
    } else {
	wm title $t "*** XCrySDen - Fermi Surface: [file tail $xcMisc(titlefile)]"
    }
    wm iconname $t "Fermi Surface"    
    wm iconbitmap . @$system(BMPDIR)/xcrysden.xbm
    
    set nb [NoteBook $t.nb -width $fs($spin,togl_w) -height $fs($spin,togl_h)]
    pack $nb -expand 1 -fill both
    set fs($spin,nb) $nb

    set fs($spin,bandlist) ""
    set fs($spin,togllist) ""

    set _first_band -1
    putsFlush stderr "NBANDS = $fs($spin,nbands)"
    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	if { $fs($spin,$i,band_selected) } {
	    if { $_first_band == -1 } {
		set _first_band $i
	    }
	    # initialize variables
	    FS_InitVar $i $spin

	    $nb insert $i band$i -text "Band #$i" \
		    -createcmd 	[list FS_RenderSurface $i $spin]

	    #
	    # page container frame
	    #
	    set f    [$nb getframe band$i]
	    set togl $f.togl$i
	    lappend fs($spin,bandlist) $i
	    lappend fs($spin,togllist) $togl

	    #
	    # toolbox frame
	    #
	    set ft [frame $f.container -relief raised -bd 1]
	    pack $ft -side top -expand 0 -fill x -padx 0m -pady 0m 
	    FS_Toolbox $ft $togl $spin $i
	    set fs($spin,$i,show_toolbox_frame)       1
	    set fs($spin,$i,toolbox_frame)            $ft
	    set fs($spin,$i,toolbox_frame_pack)       [pack info $ft]
	    set fs($spin,$i,toolbox_frame_packbefore) $togl

	    #
	    # Togl
	    #
            global toglOpt
            
	    set fs($spin,$i,togl) \
                [togl $togl \
                     -ident  $togl \
                     -rgba           true  \
                     -rgba           $toglOpt(rgba)          \
                     -redsize        $toglOpt(redsize)       \
                     -greensize      $toglOpt(greensize)     \
                     -bluesize       $toglOpt(bluesize)      \
                     -double         $toglOpt(double)        \
                     -depth          $toglOpt(depth)         \
                     -depthsize      $toglOpt(depthsize)     \
                     -accum          $toglOpt(accum)         \
                     -accumredsize   $toglOpt(accumredsize)  \
                     -accumgreensize $toglOpt(accumgreensize) \
                     -accumbluesize  $toglOpt(accumbluesize)  \
                     -accumalphasize $toglOpt(accumalphasize) \
                     -alpha          $toglOpt(alpha)         \
                     -alphasize      $toglOpt(alphasize)     \
                     -stencil        $toglOpt(stencil)       \
                     -stencilsize    $toglOpt(stencilsize)   \
                     -auxbuffers     $toglOpt(auxbuffers)    \
                     -overlay        $toglOpt(overlay)       \
                     -stereo         $toglOpt(stereo)        \
                     -time           $toglOpt(time)          \
                     -create    togl_create \
                     -display   togl_display \
                     -reshape   togl_reshape \
                     -destroy   togl_destroy \
                     -timer     togl_timer ]
            
	    pack $togl -fill both -expand 1
            update
            
	    # take care of togl's background
	    FS_UserBackground $togl
	    
	    #
	    bind $fs($spin,$i,togl) <B1-Motion>        { xc_B1motion %W %x %y }
	    bind $fs($spin,$i,togl) <B2-Motion>        { xc_B2motion %W %x %y }
	    bind $fs($spin,$i,togl) <B1-ButtonRelease> { xc_Brelease %W B1; MouseZoomBrelease %W }
	    bind $fs($spin,$i,togl) <B2-ButtonRelease> { xc_Brelease %W B2 }
	    bind $fs($spin,$i,togl) <Button-3>         [list FS_PopupMenu %W %X %Y $i $spin]
	    bind $fs($spin,$i,togl) <Shift-B1-Motion>         {  MouseZoom %W %X %Y }
	    bind $fs($spin,$i,togl) <Shift-B1-ButtonRelease>  {  MouseZoomBrelease %W }
	    
	    global tcl_platform
	    if { $tcl_platform(platform) == "unix" } {
		bind $fs($spin,$i,togl) <Button-4>  {  MouseWheelZoom %W +}
		bind $fs($spin,$i,togl) <Button-5>  {  MouseWheelZoom %W -}
	    } else {
		bind $t <MouseWheel> [list WindowsMouseWheel $fs($spin,$i,togl) %D ]
	    }

	    bind $t <Control-q>     exit_pr
	    bind $t <Control-p>     [list FSbind_printTogl $spin]
	    bind $t <Control-Alt-p> FSbind_printSetup

	    bind $t <S>         [list FSbind_SetSurfColor  $spin]
	    bind $t <C>         [list FSbind_SetCellColor  $spin]
	    bind $t <L>         [list FSbind_glLight       $spin]
	    bind $t <D>         [list FSbind_ModDepthCuing $spin]
	    bind $t <A>         [list FSbind_ModAntiAlias  $spin]

	    bind $t t [list FS_ToggleMenuCheckbutton transparent $spin FS_Config    ]
	    bind $t c [list FS_ToggleMenuCheckbutton displaycell $spin FS_fsConfig  ]
	    bind $t p [list FS_ToggleMenuCheckbutton cropbz      $spin FS_fsConfig  ]
	    bind $t d [list FS_ToggleMenuCheckbutton depthcuing  $spin FS_DepthCuing]
	    bind $t a [list FS_ToggleMenuCheckbutton antialias   $spin FS_AntiAlias ]

	    #
	    # here is some setting optimized for rendering Fermi Surfaces
	    #
	    global mody
	    #xc_setGLparam lightmodel -disable_light 1
	    xc_newvalue $togl $mody(SET_FOG_DENSITY)      2.0
	    xc_newvalue $togl $mody(SET_FOG_ORT_START_F)  0.0
	    xc_newvalue $togl $mody(SET_FOG_ORT_END_F)    0.8
	    xc_newvalue $togl $mody(SET_ANTIALIAS_DEGREE) 2
	    xc_newvalue $togl $mody(SET_ANTIALIAS_OFFSET) 0.9
	    FS_DepthCuing $i $spin

	    set fs($spin,$i,ident) [cry_surfreg $fs($spin,$i,togl)]	    
	    cry_dispfunc $fs($spin,$i,togl) fermisurface

	    #
	    # small toolbox frame ontop of Togl
	    #
	    set small_toolbox [frame $togl.f -relief raised -bd 1 -class mea]
	    place $small_toolbox -x 0 -y 0
	    set fs($spin,$i,show_small_toolbox_frame) 1
	    set fs($spin,$i,small_toolbox_frame) $small_toolbox
	    set fs($spin,$i,toolbox_frame_place) [place info $small_toolbox]

	    set separator_1 [frame $small_toolbox.s1 -height 2 -relief raised -bd 1]
	    set separator_2 [frame $small_toolbox.s2 -height 2 -relief raised -bd 1]

	    set bz_b   [radiobutton $small_toolbox.bz   -image fs_bz   -highlightthickness 1 \
			    -variable fs($spin,$i,celltype) -value bz -indicatoron 0 \
			    -selectcolor \#ff4444 -highlightbackground \#000000 \
			    -command [list FSbutton_SmallToolbox bz $i $spin]]

	    set para_b [radiobutton $small_toolbox.para -image fs_cell -highlightthickness 1 \
			    -variable fs($spin,$i,celltype) -value para -indicatoron 0 \
			    -selectcolor \#ff4444 -highlightbackground \#000000 \
			    -command [list FSbutton_SmallToolbox para $i $spin]]

	    set nocrop_b [checkbutton $small_toolbox.nocrop -image fs_nocrop -highlightthickness 1 \
			      -selectcolor \#44ff44 -highlightbackground \#000000 \
			      -variable fs($spin,$i,nocropbz) -command [list FSbutton_SmallToolbox nocrop $i $spin] -indicatoron 0]	

	    foreach button {nocell wirecell solidcell solidwirecell} {		
		set ${button}_b [radiobutton $small_toolbox.$button \
				     -image fs_$button -highlightthickness 1 \
				     -selectcolor \#4444ff -highlightbackground \#000000 \
				     -variable fs($spin,$i,radiobutton_celldisplaytype) -value $button -indicatoron 0 \
				     -command  [list FSbutton_SmallToolbox $button $i $spin]]
	    }
	    
	    set b_pack_option {-side top -fill x -padx 0 -pady 0 -ipadx 0 -ipady 0}
	    set s_pack_option {-side top -fill x -padx 0 -pady 2 -ipadx 0 -ipady 0}
	    
	    eval pack $bz_b $para_b $b_pack_option
	    eval pack $separator_1  $s_pack_option
	    eval pack $nocrop_b     $b_pack_option
	    eval pack $separator_2  $s_pack_option
	    eval pack $nocell_b $wirecell_b $solidcell_b $solidwirecell_b $b_pack_option

	    global xcFonts
	    foreach {wid text} {
	    	bz            "display Fermi surface in Brillouin zone"
	    	para          "display Fermi surface in reciprocal unit cell"
	    	nocrop        "toggle croping of Fermi surface to Brillouin zone"
	    	nocell        "do not display cell"
	    	wirecell      "display wire cell"
	    	solidcell     "display solid cell" 
	    	solidwirecell "display solid+wire cell"
	    } {
	    	set path $small_toolbox.$wid		
	    	DynamicHelp::register $path balloon $text
	    }
    
	    #
	    # Status frame
	    #
	    set ff [frame $f.f -relief ridge -bd 4]
	    pack $ff -side top -expand 0 -fill x -padx 0m -pady 0m 
	    #set fs($spin,$i,status_f) $ff

	    set fs($spin,$i,show_status_frame) 1
	    set fs($spin,$i,status_frame)      $ff
	    set fs($spin,$i,status_frame_pack) [pack info $ft]


	    if { $spin != {} } {
		set l1  [label $ff.l1 -text "Spin: [string toupper $spin]" \
			-anchor w -relief sunken -bd 1]
		pack $l1 -side left -padx 1m -ipadx 1m -ipady 1m
	    }
	    set l2 [label $ff.l2 -text "FERMI Energy: $fs(Efermi)" \
		    -anchor w -relief sunken -bd 1]
	    set l3 [label $ff.l3 -text "Min Ene: $fs($spin,$i,minE)" \
		    -anchor w -relief sunken -bd 1]
	    set l4 [label $ff.l4 -text "Max Ene: $fs($spin,$i,maxE)" \
		    -anchor w -relief sunken -bd 1]		
	    set fff [frame $ff.f -relief sunken -bd 1]
	    pack $l2 $l3 $l4 \
		    -side left -padx 1m -ipadx 1m -ipady 1m
	    pack $fff \
		    -side right -fill x -padx 1m -ipadx 1m -ipady 1m
	    
	    set l5 [label $fff.l5 -text "Isolevel:"]
	    set e  [entry $fff.e \
                        -width 15 -textvariable fs($spin,$i,isolevel) \
                        -validate key -validatecommand {string is double %P}]
	    pack $e $l5 -side right -padx 0m -pady 0m
	    bind $e <Return> [list FS_Config $i $spin]	
	}
    }
    if { $_first_band != -1 } {
	$nb raise band$_first_band
	update
	set fs($spin,toolbox_frame_height) [winfo height $fs($spin,$_first_band,toolbox_frame)]
	set fs($spin,status_frame_height)  [winfo height $fs($spin,$_first_band,status_frame)]
	xcDebug -debug "fs($spin,toolbox_frame_height) == [winfo height $fs($spin,$_first_band,toolbox_frame)]"
	xcDebug -debug "fs($spin,status_frame_height)  == [winfo height $fs($spin,$_first_band,status_frame)] "
    } else {
	WarningDialog "no band selected !!! Aplication will exit."
	exit 0
    }

    update
    set fs($spin,togl_w)      [winfo width  $fs($spin,$_first_band,togl)]
    set fs($spin,togl_h)      [winfo height $fs($spin,$_first_band,togl)]
    set fs($spin,top_w)       [winfo width $t]
    set fs($spin,top_h)       [winfo height $t]
    set fs($spin,top_togl_dw) [expr $fs($spin,top_w) - $fs($spin,togl_w)]
    set fs($spin,top_togl_dh) [expr $fs($spin,top_h) - $fs($spin,togl_h)]

    xcDebug -debug "fs($spin,togl_w)      == [winfo width  $fs($spin,$_first_band,togl)]"
    xcDebug -debug "fs($spin,togl_h)      == [winfo height $fs($spin,$_first_band,togl)]"
    xcDebug -debug "fs($spin,top_w)       == [winfo width $t]				"
    xcDebug -debug "fs($spin,top_h)       == [winfo height $t]				"
    xcDebug -debug "fs($spin,top_togl_dw) == [expr $fs($spin,top_w) - $fs($spin,togl_w)]"
    xcDebug -debug "fs($spin,top_togl_dh) == [expr $fs($spin,top_h) - $fs($spin,togl_h)]"

    xcDebug -debug "Notebook's width:  [winfo width  $nb]" 
    xcDebug -debug "Notebook's height: [winfo height $nb]" 
 
    #if { $fs($spin,top_h) < $fs($spin,togl_h) } {
    #	set fs($spin,top_h) [expr $fs($spin,togl_h) + $fs($spin,toolbox_frame_height) + $fs($spin,status_frame_height) + 30]
    #	set fs($spin,top_w) [expr $fs($spin,togl_w) + 4]
    #}
    #set w $fs($spin,top_w)
    #set h $fs($spin,top_h)
    #if { $spin == "dn" } {
    #	wm geometry $t ${w}x${h}+0-0
    #} else {
    #	wm geometry $t ${w}x${h}-0+0
    #}

    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	if { $fs($spin,$i,band_selected) } {
	    set fs($spin,$i,togl_w)      [winfo width  $fs($spin,$_first_band,togl)]
	    set fs($spin,$i,togl_h)      [winfo height $fs($spin,$_first_band,togl)]
	    set fs($spin,$i,top_togl_dw) [expr $fs($spin,top_w) - $fs($spin,togl_w)]
	    set fs($spin,$i,top_togl_dh) [expr $fs($spin,top_h) - $fs($spin,togl_h)]
	}
    }

    FS_Multi $nb $spin

    bind $t <Configure> [list FS_ResizeWin %W %w %h $t $spin]
}

proc FS_PopupMenu {W x y i {spin {}} {multiband {}}} {
    global fs

    set togl $fs($spin,$i,togl)

    if { [winfo exists $W.menu] } {
	destroy $W.menu
    }
    set m [menu $W.menu -tearoff 1]
    tk_popup $m $x $y

    #$m add command -label "PopUp Menu" -state disabled
    #$m add separator

    # ------------------------------------------------------------------------
    # Pop-Up menu
    # ------------------------------------------------------------------------

    if { $multiband == "" } {
	$m add command -label "Render Surface" \
	    -command [list FS_RenderSurface $i $spin]

    #$m add separator
    #$m add command -label "Interpolation" \
    #	-command [list FS_Interpolation  "Interpolation for band \# $i:" $togl $spin $i]
    #$m add command -label "Zoom" \
    #	-command [list toglZoom "Zoom for band \# $i:" $togl]
    #
	$m add separator
    }

    # Palette-cascade
    $m add cascade -image colors      -menu $m.colors 
    menu $m.colors -tearoff 1
    ColorMenu $W $m.colors
    
    # File-cascade
    $m add cascade -label "File"      -menu $m.file
    set mfile [menu $m.file -tearoff 1]
    $mfile add command -label "Save Fermi Surface(s) in BXSF format" \
	-command [list FS_SaveBXSF $i $spin multiband]
    $mfile add separator
    $mfile add command -label "Print Setup" -command printSetup -accelerator "Ctlr-Alt-p"
    $mfile add command -label "Print " -command [list printTogl $togl] -accelerator "Ctrl-p"

    if { $multiband == "" } {
	# View-cascade
	$m add cascade -label "View"    -menu $m.view
	menu $m.view -tearoff 1

	# Display-cascade
	$m add cascade -label "Display" -menu $m.dis 
	menu $m.dis -tearoff 1

	# Modify-cascade
	$m add cascade -label "Modify"  -menu $m.mody
	menu $m.mody -tearoff 1

	FS_ViewMenu    $m.view $W $i $spin
	FS_DisplayMenu $m.dis  $W $i $spin
	FS_ModifyMenu  $m.mody $W $i $spin

	#$m add cascade -label "Modify"    -menu $m.mod 
	#menu $m.mod -tearoff 1
	#FS_ModifyMenu $m.dis $W $i $spin
	#$m add cascade -label "Tools"     -menu $m.tools 
	#menu $m.tools -tearoff 0
	#FS_ToolsMenu $m.dis $W $i $spin
    }

    $m add separator
    $m add command -label "Print " -command [list printTogl $togl] -accelerator "Ctrl-p"

    $m add separator

    $m add command -label "Exit" -command exit_pr -accelerator "Ctrl-q"
}

proc FS_ResizeWin {W w h t {spin {}}} {
    global fs

    if { $t != $W } {
	set w [winfo width  $t]
	set h [winfo height $t]
    }

    #xcDebug -debug "FS_ResizeWin> (w,h) == ($w,$h)"

    # update only if size of "toplevel ." has changed
    #if { $w != $fs($spin,top_w) || $h != $fs($spin,top_h) } {
    #	set fs($spin,top_w) $w
    #	set fs($spin,top_h) $h
    #	for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
    #	    if { $fs($spin,$i,band_selected) } {
    #		set fs($spin,$i,togl_w) [expr $w - $fs($spin,$i,top_togl_dw)]
    #		set fs($spin,$i,togl_h) [expr $h - $fs($spin,$i,top_togl_dh)]
    #		$fs($spin,$i,togl) config \
    #		    -width  $fs($spin,$i,togl_w) \
    #		    -height $fs($spin,$i,togl_h)
    #	    }
    #	}
    #}
}

proc FS_InitVar {i {spin {}}} {
    global fs
    
    #
    # set monocolor
    # 
    set rainbow {
	{ 1.0 0.2 0.2 0.5 }
	{ 1.0 1.0 0.2 0.5 }
	{ 0.2 1.0 0.2 0.5 }
	{ 0.2 1.0 1.0 0.5 }
	{ 0.2 0.2 1.0 0.5 }
	{ 1.0 0.2 1.0 0.5 }
    }
    foreach rgb $rainbow {
	set r [expr 1.0 - [lindex $rgb 0]]
	set g [expr 1.0 - [lindex $rgb 1]]
	set b [expr 1.0 - [lindex $rgb 2]]
	set a 0.5
	lappend backrainbow [list $r $g $b $a]
    }

    #
    # hard-coded defaults
    #
    set im [expr $i - 6 * int( $i / 6)]
    set fs($spin,$i,celltype)           bz
    set fs($spin,$i,text_celltype)      "first Brillouin zone"
    set fs($spin,$i,cropbz)             1
    set fs($spin,$i,nocropbz)           0
    set fs($spin,$i,displaycell)        1
    set fs($spin,$i,celldisplaytype)    wire
    set fs($spin,$i,drawstyle)          solid
    set fs($spin,$i,transparent)        0
    set fs($spin,$i,shademodel)         smooth
    set fs($spin,$i,colormodel)         "set front-side color only"
    set fs($spin,$i,monocolor)          [lindex $rainbow $im]
    set fs($spin,$i,backmonocolor)      [lindex $backrainbow $im]
    set fs($spin,$i,smoothsteps)        0
    set fs($spin,$i,smoothweight)       0.2
    set fs($spin,$i,interpolationdegree) {1 1 1}
    set fs($spin,$i,frontface)          CW
    set fs($spin,$i,revertnormals)      0
    
    set fs($spin,$i,wirecellcolor)  {1.00 1.00 1.00 1.00}
    set fs($spin,$i,solidcellcolor) {0.00 0.95 0.95 0.40}
    set fs($spin,$i,antialias)      0
    set fs($spin,$i,depthcuing)     0

    set fs($spin,$i,radiobutton_celldisplaytype) $fs($spin,$i,celldisplaytype)cell
    
    # try to use user specified defaults   
    FS_UserDefaults $i $spin    

    set fs($spin,$i,old_celltype)        $fs($spin,$i,celltype)
    set fs($spin,$i,old_cropbz)          $fs($spin,$i,cropbz)
    set fs($spin,$i,old_displaycell)     $fs($spin,$i,displaycell)
    set fs($spin,$i,old_celldisplaytype) $fs($spin,$i,celldisplaytype)
    set fs($spin,$i,old_drawstyle)       $fs($spin,$i,drawstyle)
    set fs($spin,$i,old_transparent)     $fs($spin,$i,transparent)
    set fs($spin,$i,old_shademodel)      $fs($spin,$i,shademodel)
    set fs($spin,$i,old_monocolor)       $fs($spin,$i,monocolor)
    set fs($spin,$i,old_smoothsteps)     $fs($spin,$i,smoothsteps)
    set fs($spin,$i,old_smoothweight)    $fs($spin,$i,smoothweight)
    set fs($spin,$i,old_interpolationdegree) $fs($spin,$i,interpolationdegree)    
    #set fs($spin,$i,old_frontface)       $fs($spin,$i,frontface)
    #set fs($spin,$i,old_revertnormals)   $fs($spin,$i,revertnormals)
}


proc FS_UserBackground {togl} {	    
    global myParam mody
    # BEWARE: the  myParam(FS_BACKGROUND) needs a special treatment    
    if { [info exists myParam(FS_BACKGROUND)] } {
	if { ! [rgba $myParam(FS_BACKGROUND)] } {
	    error "wrong value \"$myParam(FS_BACKGROUND)\" for myParam(FS_BACKGROUND), should be one of rgba type; correct custom-definition file"
	} else {
	    eval xc_newvalue $togl $mody(L_BACKGROUND) $myParam(FS_BACKGROUND)
	}
    }
}


proc FS_UserDefaults {i {spin {}}} {
    global fs myParam

    # BEWARE: the  myParam(FS_BACKGROUND) needs a special treatment

    foreach {fs_item allowed} {
	FS_CELLTYPE            {bz para}
	FS_CROPBZ              {0 1} 
	FS_CELLDISPLAYTYPE     {none wire solid solidwire}
	FS_DRAWSTYLE           {solid wire dot}
	FS_TRANSPARENT         {0 1}
	FS_SHADEMODEL          {smooth flat}
	FS_INTERPOLATIONDEGREE {@ positiveInteger}
	FS_FRONTFACE           {CW CCW}
	FS_REVERTNORMALS       {0 1}
	FS_WIRECELLCOLOR       {@ rgba}
	FS_SOLIDCELLCOLOR      {@ rgba}
	FS_ANTIALIAS           {0 1}
	FS_DEPTHCUING          {0 1}
    } {
	if { [info exists myParam($fs_item)] } {

	    regsub ^FS_ $fs_item {} _item
	    set item [string tolower $_item]
	
	    if { [lindex $allowed 0] == "@" } {
		
		# not a literal comparison, but a given type is specified
		
		set typeCmd [lindex $allowed 1]

		if { ! [$typeCmd $myParam($fs_item)] } {
		    error "wrong value \"$myParam($fs_item)\" for myParam($fs_item), should be one of $typeCmd type; correct custom-definition file"
		} else {
		    # special treatment for FS_INTERPOLATIONDEGREE
		    if { $fs_item eq "FS_INTERPOLATIONDEGREE" } {
			set fs($spin,$i,$item) [list $myParam($fs_item) $myParam($fs_item) $myParam($fs_item)]
		    } else {
			set fs($spin,$i,$item) $myParam($fs_item)
		    }
		}

	    } else {

		# literal comparison
		
		if { ! [allowedValue $myParam($fs_item) $allowed] } {
		    error "wrong value \"$myParam($fs_item)\" for myParam($fs_item), should be one of: $allowed; correct custom-definition file"
		} else {		    
		    set fs($spin,$i,$item) $myParam($fs_item)		    
		    
		    # handle specialties

		    if { $item eq "cropbz" } {

			set fs($spin,$i,nocropbz)  [expr ! $fs($spin,$i,cropbz)]

		    } elseif { $item eq "celldisplaytype" } {

			switch -- $myParam($fs_item) {
			    none {
				set fs($spin,$i,displaycell)     0
				set fs($spin,$i,celldisplaytype) wire; # just in any case !!!
				set fs($spin,$i,radiobutton_celldisplaytype) nocell
			    }
			    wire - solid - solidwire {
				set fs($spin,$i,displaycell)     1
				set fs($spin,$i,radiobutton_celldisplaytype) $fs($spin,$i,celldisplaytype)cell
			    }
			}
		    }
		}
	    }
	}
    }
}


proc FS_RenderSurface {i {spin {}}} {
    global fs

    if { ! [info exist fs($spin,$i,rendered)] } {
	set fs($spin,$i,rendered) 1
	
	FS:cry_surf $i $spin
	# next lines are a hack-around a "display-bug" to force the display
	set w [lindex [$fs($spin,$i,togl) config -width] end]
	$fs($spin,$i,togl) config -width $w
	$fs($spin,$i,togl) render
	$fs($spin,$i,togl) swapbuffers
	update
    }
}

proc FS:cry_surf {i {spin {}}} {
    global fs

    # band index-identifiers strats from 0 not from 1 in ReadBandGrid !!!
    set iband [expr $i - 1]
    SetWatchCursor
    update

    if { $fs($spin,$i,colormodel) == "set front-side color only" } {

	# monocolor == -monocolor $fs($spin,$i,monocolor)

	cry_surf $fs($spin,$i,togl) \
	    -ident $fs($spin,$i,ident) \
	    -type  fermisurface \
	    -fs    [list \
			-gridindex       $fs($spin,grid_index) \
			-gridsubindex    $fs($spin,grid_subindex) \
			-bandindex       $iband \
			-celltype        $fs($spin,$i,celltype) \
			-cropbz          $fs($spin,$i,cropbz)  \
			-displaycell     $fs($spin,$i,displaycell) \
			-celldisplaytype $fs($spin,$i,celldisplaytype) \
			-interpolationdegree $fs($spin,$i,interpolationdegree) \
			-wirecellcolor   $fs($spin,$i,wirecellcolor) \
			-solidcellcolor  $fs($spin,$i,solidcellcolor)] \
	    -level        $fs($spin,$i,isolevel) \
	    -drawstyle    $fs($spin,$i,drawstyle) \
	    -transparent  $fs($spin,$i,transparent) \
	    -shademodel   $fs($spin,$i,shademodel) \
	    -monocolor    $fs($spin,$i,monocolor) \
	    -smoothsteps  $fs($spin,$i,smoothsteps) \
	    -smoothweight $fs($spin,$i,smoothweight) \
	    -frontface    $fs($spin,$i,frontface) \
	    -revertnormals $fs($spin,$i,revertnormals)
    } else {

	# monocolor == -frontmonocolor $fs($spin,$i,monocolor) \
	#              -backmonocolor $fs($spin,$i,backmonocolor)"

	cry_surf $fs($spin,$i,togl) \
	    -ident $fs($spin,$i,ident) \
	    -type  fermisurface \
	    -fs    [list \
			-gridindex       $fs($spin,grid_index) \
			-gridsubindex    $fs($spin,grid_subindex) \
			-bandindex       $iband \
			-celltype        $fs($spin,$i,celltype) \
			-cropbz          $fs($spin,$i,cropbz)  \
			-displaycell     $fs($spin,$i,displaycell) \
			-celldisplaytype $fs($spin,$i,celldisplaytype) \
			-interpolationdegree $fs($spin,$i,interpolationdegree) \
			-wirecellcolor   $fs($spin,$i,wirecellcolor) \
			-solidcellcolor  $fs($spin,$i,solidcellcolor)] \
	    -level        $fs($spin,$i,isolevel) \
	    -drawstyle    $fs($spin,$i,drawstyle) \
	    -transparent  $fs($spin,$i,transparent) \
	    -shademodel   $fs($spin,$i,shademodel) \
	    -frontmonocolor $fs($spin,$i,monocolor) \
	    -backmonocolor  $fs($spin,$i,backmonocolor) \
	    -smoothsteps  $fs($spin,$i,smoothsteps) \
	    -smoothweight $fs($spin,$i,smoothweight) \
	    -frontface    $fs($spin,$i,frontface) \
	    -revertnormals $fs($spin,$i,revertnormals)
    }

    ResetCursor 
    update
}


proc FS:cry_surfconfig {i {spin {}}} {
    global fs

    if { ! [info exist fs($spin,$i,rendered)] } {
	return
    }

    # band index-identifiers strats from 0 not from 1 in ReadBandGrid !!!
    set iband [expr $i - 1]

    SetWatchCursor
    update

    if { $fs($spin,$i,colormodel) == "set front-side color only" } {

	# monocolor == -monocolor $fs($spin,$i,monocolor) 
	cry_surfconfig $fs($spin,$i,togl) \
	    -ident $fs($spin,$i,ident) \
	    -fs    [list \
			-gridindex       $fs($spin,grid_index) \
			-gridsubindex    $fs($spin,grid_subindex) \
			-bandindex       $iband \
			-celltype        $fs($spin,$i,celltype) \
			-cropbz          $fs($spin,$i,cropbz)  \
			-displaycell     $fs($spin,$i,displaycell) \
			-celldisplaytype $fs($spin,$i,celldisplaytype) \
			-interpolationdegree $fs($spin,$i,interpolationdegree) \
			-wirecellcolor   $fs($spin,$i,wirecellcolor) \
			-solidcellcolor  $fs($spin,$i,solidcellcolor)] \
	    -render       1 \
	    -level        $fs($spin,$i,isolevel) \
	    -drawstyle    $fs($spin,$i,drawstyle) \
	    -transparent  $fs($spin,$i,transparent) \
	    -shademodel   $fs($spin,$i,shademodel) \
	    -monocolor    $fs($spin,$i,monocolor) \
	    -smoothsteps  $fs($spin,$i,smoothsteps) \
	    -smoothweight $fs($spin,$i,smoothweight) \
	    -frontface    $fs($spin,$i,frontface) \
	    -revertnormals $fs($spin,$i,revertnormals)
    } else {
	# monocolor == -frontmonocolor $fs($spin,$i,monocolor) \
	#              -backmonocolor $fs($spin,$i,backmonocolor)"

	cry_surfconfig $fs($spin,$i,togl) \
	    -ident $fs($spin,$i,ident) \
	    -fs    [list \
			-gridindex       $fs($spin,grid_index) \
			-gridsubindex    $fs($spin,grid_subindex) \
			-bandindex       $iband \
			-celltype        $fs($spin,$i,celltype) \
			-cropbz          $fs($spin,$i,cropbz)  \
			-displaycell     $fs($spin,$i,displaycell) \
			-celldisplaytype $fs($spin,$i,celldisplaytype) \
			-interpolationdegree $fs($spin,$i,interpolationdegree) \
			-wirecellcolor   $fs($spin,$i,wirecellcolor) \
			-solidcellcolor  $fs($spin,$i,solidcellcolor)] \
	    -render       1 \
	    -level        $fs($spin,$i,isolevel) \
	    -drawstyle    $fs($spin,$i,drawstyle) \
	    -transparent  $fs($spin,$i,transparent) \
	    -shademodel   $fs($spin,$i,shademodel) \
	    -frontmonocolor $fs($spin,$i,monocolor) \
	    -backmonocolor  $fs($spin,$i,backmonocolor) \
	    -smoothsteps  $fs($spin,$i,smoothsteps) \
	    -smoothweight $fs($spin,$i,smoothweight) \
	    -frontface    $fs($spin,$i,frontface) \
	    -revertnormals $fs($spin,$i,revertnormals)	
    }

    ResetCursor 
    update
}

proc FS_ModifyMenu {m togl i {spin {}}} {
    global fs
    
    $m add command -label "Surface Color" \
	-command [list FS_SetSurfColor $i $spin] -accelerator "Shift-s"

    $m add command -label "Cell Color" \
	-command [list FS_SetCellColor $i $spin] -accelerator "Shift-c"

    $m add separator
    $m add command -label "Lighting Parameters" -command [list glLight $togl] \
	-accelerator "Shift-l"

    $m add command -label "Depth-Cuing Parameters" \
	-command [list FS_ModDepthCuing $i $spin] -accelerator "Shift-d"

    $m add command -label "Anti-aliasing Parameters" \
	-command [list FS_ModAntiAlias $i $spin] -accelerator "Shift-a"
}

proc FS_ViewMenu {m togl i {spin {}}} {
    global fs

    #
    # Checkbuttons
    #
    $m add checkbutton -label "Show Toolbox" \
	-variable fs($spin,$i,show_toolbox_frame) \
	-command [list FS_ViewMenu:_show toolbox $i $spin]
    $m add checkbutton -label "Show Small Toolbox" \
	-variable fs($spin,$i,show_small_toolbox_frame) \
	-command [list FS_ViewMenu:_show small_toolbox $i $spin]
    $m add checkbutton -label "Show Status Frame" \
	-variable fs($spin,$i,show_status_frame) \
	-command [list FS_ViewMenu:_show status $i $spin]
}

proc FS_ViewMenu:_show {which i spin} {
    global fs

    set dh 0

    switch -exact -- $which {
	toolbox {
	    #
	    # TOOLBOX
	    #
	    if { $fs($spin,$i,show_toolbox_frame) } {
		eval pack $fs($spin,$i,toolbox_frame) $fs($spin,$i,toolbox_frame_pack) \
		    -before $fs($spin,$i,toolbox_frame_packbefore)
		set dh [expr -1 * $fs($spin,toolbox_frame_height)]
	    } else {
		pack forget $fs($spin,$i,toolbox_frame)
		set dh $fs($spin,toolbox_frame_height)
	    }
	}

	small_toolbox {
	    #
	    # SMALL-TOOLBOX
	    #
	    if { $fs($spin,$i,show_small_toolbox_frame) } {
		eval place $fs($spin,$i,small_toolbox_frame) $fs($spin,$i,toolbox_frame_place)
	    } else {
		place forget $fs($spin,$i,small_toolbox_frame)	
	    }
	}
	
	status {
	    #
	    # STATUS-FRAME
	    #
	    if { $fs($spin,$i,show_status_frame) } {
		eval pack $fs($spin,$i,status_frame) $fs($spin,$i,status_frame_pack)
		set dh [expr -1 * $fs($spin,status_frame_height)]
	    } else {
		pack forget $fs($spin,$i,status_frame)
		set dh $fs($spin,status_frame_height)
	    }
	}
    }

    set h [winfo height $fs($spin,$i,togl)]
    set fs($spin,$i,top_togl_dh) [expr $fs($spin,$i,top_togl_dh) - $dh]	
    set fs($spin,$i,togl_h)      [expr $h + $dh]

    $fs($spin,$i,togl) config -height $fs($spin,$i,togl_h)	    
    FS_ResizeWin . 0 0 [winfo toplevel $fs($spin,$i,togl)] $spin
}


proc FS_DisplayMenu {m togl i {spin {}}} {
    global fs

    if { $fs($spin,$i,celltype) == "para" } {
	set cell cell
    } else {
	set cell "first Brillouin zone"
    }

    #
    # Checkbuttons
    #
    $m add checkbutton -label "Transparent Fermi Surface" \
	    -variable fs($spin,$i,transparent) \
	    -command [list FS_Config $i $spin] -accelerator "t"

    $m add checkbutton -label "Display $cell" \
	-variable fs($spin,$i,displaycell) \
	-command [list FS_fsConfig $i $spin] -accelerator "c"

    $m add checkbutton -label "Crop Fermi Surface to first BZ" \
	-variable fs($spin,$i,cropbz) \
	-command [list FS_fsConfig $i $spin] -accelerator "p"

    if { $fs($spin,$i,celltype) == "para" } {
	$m entryconfig "Crop Fermi Surface to first BZ" -state disabled
    }
    $m add separator
    $m add checkbutton -label "Depth-Cuing" \
	-variable fs($spin,$i,depthcuing) -onvalue 1 -offvalue 0 \
	-command [list FS_DepthCuing $i $spin] -accelerator "d"

    $m add checkbutton -label "Anti-Aliasing" \
	-variable fs($spin,$i,antialias) -onvalue 1 -offvalue 0 \
	-command [list FS_AntiAlias $i $spin] -accelerator "a"

    #
    # CASCADES 
    #
    $m add separator
    $m add cascade -label "Cell type ..." -menu $m.celltype
    $m add cascade -label "Display $cell as ..." -menu $m.discell    
    $m add cascade -label "Surface Drawstyle ..."  -menu $m.draw
    $m add cascade -label "Surface Shademodel ..." -menu $m.shade

    $m add separator

    #$m add command -label "Surface Smoothing" \
    #	    -command [list FS_SurfSmooth $i $spin]


    # CELLTYPE CASCADE
    menu $m.celltype -tearoff 0
    $m.celltype add radiobutton -label "first Brillouin zone" \
	    -variable fs($spin,$i,text_celltype) \
	    -command [list celltype:FS_fsConfig $i $spin]
    $m.celltype add radiobutton -label "reciprocal primitive cell" \
	    -variable fs($spin,$i,text_celltype) \
	    -command [list celltype:FS_fsConfig $i $spin]
    
    # DISPLAYCELL CASCADE
    menu $m.discell -tearoff 0
    $m.discell add radiobutton -label "solid" \
	    -variable fs($spin,$i,celldisplaytype) \
	    -command [list FS_fsConfig $i $spin]
    $m.discell add radiobutton -label "wire" \
	    -variable fs($spin,$i,celldisplaytype) \
	    -command [list FS_fsConfig $i $spin]
    #$m.discell add radiobutton -label "rod" \
    #	    -variable fs($spin,$i,celldisplaytype) \
    #	    -command [list FS_fsConfig $i $spin]
    $m.discell add radiobutton -label "solidwire" \
	    -variable fs($spin,$i,celldisplaytype) \
	    -command [list FS_fsConfig $i $spin]
    #$m.discell add radiobutton -label "solidrod" \
    #	    -variable fs($spin,$i,celldisplaytype) \
    #	    -command [list FS_fsConfig $i $spin]

    #$m.discell entryconfig "rod"      -state disabled
    #$m.discell entryconfig "solidrod" -state disabled
    #########################################################################
    #/

    # DRAWSTYLE CASCADE
    menu $m.draw -tearoff 0
    $m.draw add radiobutton -label "solid" \
	    -variable fs($spin,$i,drawstyle) \
	    -command [list FS_Config $i $spin]
    $m.draw add radiobutton -label "wire" \
	    -variable fs($spin,$i,drawstyle) \
	    -command [list FS_Config $i $spin]
    $m.draw add radiobutton -label "dot" \
	    -variable fs($spin,$i,drawstyle) \
	    -command [list FS_Config $i $spin]

    # SHADEMODEL CASCADE
    menu $m.shade -tearoff 0
    $m.shade add radiobutton -label "smooth" \
	    -variable fs($spin,$i,shademodel) \
	    -command [list FS_Config $i $spin]
    $m.shade add radiobutton -label "flat" \
	    -variable fs($spin,$i,shademodel) \
	    -command [list FS_Config $i $spin]
}


proc FS_fsConfig {i {spin {}}} {
    global fs

    if { $fs($spin,$i,old_celltype) != $fs($spin,$i,celltype) } {
	FS:cry_surf $i $spin
    } else {
	FS:cry_surfconfig $i $spin
    }
    set fs($spin,$i,old_celltype) $fs($spin,$i,celltype)
}

proc FS_Config {i {spin {}}} {
    global fs

    # -level        $fs($spin,$i,isolevel)
    # -drawstyle    $fs($spin,$i,drawstyle)
    # -transparent  $fs($spin,$i,transparent)
    # -shademodel   $fs($spin,$i,shademodel)
    # -monocolor    $fs($spin,$i,monocolor)
    # -smoothsteps  $fs($spin,$i,smoothsteps)
    # -smoothweight $fs($spin,$i,smoothweight)

    FS:cry_surfconfig $i $spin    
}
    

proc celltype:FS_fsConfig {i {spin {}}} {
    global fs

    if { $fs($spin,$i,text_celltype) == "reciprocal primitive cell" } {
	set fs($spin,$i,celltype) para
    } else {
	set fs($spin,$i,celltype) bz
    }

    FS_fsConfig $i $spin
}

 
proc FS_SaveBXSF {i {spin {}} {multiband {}}} {
    global fs system
    
    set filetypes {
	{{BXSF}         {.bxsf} }
	{{All Files}     *      }
    }        
    set sfile [tk_getSaveFile \
		   -initialdir       $system(PWD) \
		   -title            "Save BXSF File" \
		   -defaultextension ".bxsf" \
		   -filetypes        $filetypes]
    if { $sfile == "" } {
	return
    }

    if { $multiband == "" } {
	# band index-identifiers strats from 0 not from 1 in ReadBandGrid !!!
	xc_writebandXSF $fs($spin,$i,ident) $fs(Efermi) $i $sfile
    } else {
	foreach band $fs($spin,bandlist) {
	    xc_writebandXSF $fs($spin,$band,ident) $fs(Efermi) $band ${sfile}.band-$band
	}
    }    
}


proc FS_SetSurfColor {i {spin {}}} {
    global fs

    if { ! [info exists fs($spin,$i,monocolor)] } {
	return
    }

    set fs($spin,$i,monocolor_R) [lindex $fs($spin,$i,monocolor) 0]
    set fs($spin,$i,monocolor_G) [lindex $fs($spin,$i,monocolor) 1]
    set fs($spin,$i,monocolor_B) [lindex $fs($spin,$i,monocolor) 2]
    set fs($spin,$i,monocolor_A) [lindex $fs($spin,$i,monocolor) 3]
        
    set fs($spin,$i,backmonocolor_R) [lindex $fs($spin,$i,backmonocolor) 0]
    set fs($spin,$i,backmonocolor_G) [lindex $fs($spin,$i,backmonocolor) 1]
    set fs($spin,$i,backmonocolor_B) [lindex $fs($spin,$i,backmonocolor) 2]
    set fs($spin,$i,backmonocolor_A) [lindex $fs($spin,$i,backmonocolor) 3]

    if { $spin == "" } {
	set t [xcToplevel [WidgetName] "Surface Colors for band #$i" "Surface Colors" . 0 0 1]
    } else {
	set t [xcToplevel [WidgetName] "Surface Colors for band #$i (spin: $spin)" "Surface Colors" . 0 0 1]
    }

    #
    # widgets
    # 
    set f1  [frame $t.f1]
    set f2  [frame $t.f2]
    set f21 [frame $f2.1 -relief groove -bd 2]
    set f22 [frame $f2.2 -relief groove -bd 2]
    set f23 [frame $f2.3]
    pack $f1 $f2 -side top -padx 5 -pady 5 -fill both -expand 1    
    pack $f21 $f22 $f23 -side left -padx 3 -pady 3 -fill both

    set fs($spin,$i,backcolor_frame) $f22

    RadioBut $f1 "Color model:" fs($spin,$i,colormodel) top top 1 0 \
	"set front-side color only" "set front- and back-side colors"

    setRGBAwidget $f21 "Front-side color:" \
	fs($spin,$i,monocolor_R) fs($spin,$i,monocolor_G) \
	fs($spin,$i,monocolor_B) fs($spin,$i,monocolor_A) \
	_UNKNOWN_
    
    setRGBAwidget $f22 "Back-side color:" \
	fs($spin,$i,backmonocolor_R) fs($spin,$i,backmonocolor_G) \
	fs($spin,$i,backmonocolor_B) fs($spin,$i,backmonocolor_A) \
	_UNKNOWN_

    trace variable fs($spin,$i,colormodel) w FS_SetSurfColor:_widget
    FS_SetSurfColor:_widget fs $spin,$i,colormodel w
    
    #
    # in bottom frame goes the "Close|Update" buttons
    #    
    set update  [button $f23.update -text "Update" -command [list FS_SetSurfColor:Update $i $spin]]
    set close   [button $f23.close  -text "Close"  -command [list CancelProc $t]]
    pack $update $close -side top -padx 5 -pady 5 -ipadx 3 -ipady 3 -fill x
}


proc FS_SetSurfColor:Update {i spin} {
    global fs
    
    set fs($spin,$i,monocolor)     [list \
					$fs($spin,$i,monocolor_R) \
					$fs($spin,$i,monocolor_G) \
					$fs($spin,$i,monocolor_B) \
					$fs($spin,$i,monocolor_A)]
    
    set fs($spin,$i,backmonocolor) [list \
					$fs($spin,$i,backmonocolor_R) \
					$fs($spin,$i,backmonocolor_G) \
					$fs($spin,$i,backmonocolor_B) \
					$fs($spin,$i,backmonocolor_A)]
    
    FS_Config $i $spin
}


proc FS_SetSurfColor:_widget {name1 name2 op} {
    global fs

    regsub -- {,colormodel$} $name2 {} spin_i
    
    if { $fs($name2) == "set front-side color only" } {
	xcDisableAll -disabledfg $fs($spin_i,backcolor_frame)
    } else {
	xcEnableAll -disabledfg $fs($spin_i,backcolor_frame)
    }
}


proc FS_SetCellColor {i {spin {}}} {
    global fs

    set t .fs_cellcolor
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Cell Color" "Cell Color"

    set fs($spin,$i,wirecellcolor_R) [lindex $fs($spin,$i,wirecellcolor) 0]
    set fs($spin,$i,wirecellcolor_G) [lindex $fs($spin,$i,wirecellcolor) 1]
    set fs($spin,$i,wirecellcolor_B) [lindex $fs($spin,$i,wirecellcolor) 2]
    set fs($spin,$i,wirecellcolor_A) [lindex $fs($spin,$i,wirecellcolor) 3]

    set fs($spin,$i,solidcellcolor_R) [lindex $fs($spin,$i,solidcellcolor) 0]
    set fs($spin,$i,solidcellcolor_G) [lindex $fs($spin,$i,solidcellcolor) 1]
    set fs($spin,$i,solidcellcolor_B) [lindex $fs($spin,$i,solidcellcolor) 2]
    set fs($spin,$i,solidcellcolor_A) [lindex $fs($spin,$i,solidcellcolor) 3]

    set f1 [frame $t.1]
    set f2 [frame $t.2]
    pack $f1 $f2 -side left -fill both -padx 5 -pady 5

    foreach type {wire solid} {
	set frame($type)  [frame $f1.$type -relief groove -bd 2]
	pack $frame($type) -side left -fill both -padx 5 -pady 0 -ipady 3 -expand 1	    
	
	setRGBAwidget $frame($type) "[string totitle $type]-cell color:" \
	    fs($spin,$i,${type}cellcolor_R) fs($spin,$i,${type}cellcolor_G) \
	    fs($spin,$i,${type}cellcolor_B) fs($spin,$i,${type}cellcolor_A) \
	    _UNKNOWN_
    }
    
    #
    # in bottom frame goes the "Close|Update" buttons
    #    
    set update  [button $f2.update -text "Update" -command [list FS_SetCellColor:Update $i $spin]]
    set close   [button $f2.close  -text "Close"  -command [list CancelProc $t]]
    pack $update $close -side top -padx 5 -pady 5 -ipadx 3 -ipady 3 -fill x
}
proc FS_SetCellColor:Update {i spin} {
    global fs
    
    set fs($spin,$i,wirecellcolor) [list \
					$fs($spin,$i,wirecellcolor_R) \
					$fs($spin,$i,wirecellcolor_G) \
					$fs($spin,$i,wirecellcolor_B) \
					$fs($spin,$i,wirecellcolor_A)]
    set fs($spin,$i,solidcellcolor) [list \
					$fs($spin,$i,solidcellcolor_R) \
					$fs($spin,$i,solidcellcolor_G) \
					$fs($spin,$i,solidcellcolor_B) \
					$fs($spin,$i,solidcellcolor_A)]
    FS_Config $i $spin
}

proc FS_SurfSmooth {i {spin {}}} {
    global fs fs_trial

    set fs_trial($spin,$i,smoothsteps)   $fs($spin,$i,smoothsteps) 
    set fs_trial($spin,$i,smoothweight)	 $fs($spin,$i,smoothweight)

    set t [xcToplevel [WidgetName] "Surface Smoothing" "SurfSmooth" . 20 20 1]

    message $t.m -aspect 800 \
	-relief groove -bd 2 \
	-text "Reasonable values for weight are between 0.1 and 1. Lighter weight will require more steps for smoothing, but will perturb the surface less !!!"
    pack $t.m -side top -padx 3m -pady 3m -ipadx 1m -ipady 1m 

    set f [frame $t.f]	        
    set e [FillEntries $t {"Smoothing steps:" "Smoothing weight:"} \
	    [list fs_trial($spin,$i,smoothsteps) \
	    fs_trial($spin,$i,smoothweight)] 17 7]
    set e1 [string trimright $e 1]
    set foclist "$e $e1"
    set varlist [list \
	[list fs_trial($spin,$i,smoothsteps) int] \
	[list fs_trial($spin,$i,smoothweight) real] ]

    button $t.b1 -text "Close"  -command [list CancelProc $t]
    button $t.b2 -text "Update" \
	    -command [list FS_SurfSmoothOK $t $foclist $varlist $i $spin]
    
    pack $f -side bottom -expand 1 -fill both  -padx 3m -pady 3m
    pack $t.b1 $t.b2 -side left -expand 1 -padx 2m -pady 2m
}
proc FS_SurfSmoothOK {t foclist varlist i {spin {}}} {
    global fs fs_trial
    
    if ![check_var $varlist $foclist] {
	return
    }
    set fs($spin,$i,smoothsteps)   $fs_trial($spin,$i,smoothsteps) 
    set fs($spin,$i,smoothweight)  $fs_trial($spin,$i,smoothweight)
    FS_Config $i $spin

    return
}


proc FS_AntiAlias {i {spin {}}} {
    global fs mody

    xc_newvalue $fs($spin,$i,togl) $mody(L_ANTIALIAS) $fs($spin,$i,antialias)

    # update display
    $fs($spin,$i,togl) render
}


proc FS_DepthCuing {i {spin {}}} {
    global fs mody

    xc_newvalue $fs($spin,$i,togl) $mody(L_FOG) $fs($spin,$i,depthcuing)

    # update display
    $fs($spin,$i,togl) render
}

proc FS_ModAntiAlias {i {spin {}}} {
    global fs
    set t .fs_antialias
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Anti-aliasing Parameters" "Antialias"
    glModParam:AntiAlias $t $t $fs($spin,$i,togl)
}

proc FS_ModDepthCuing {i {spin {}}} {
    global fs
    set t .fs_depthcuing
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Depth-Cuing Parameters" "Depth-Cuing"
    glModParam:DepthCuing $t $t $fs($spin,$i,togl)
}


proc FS_ToggleMenuCheckbutton {what spin cmd} {
    global fs

    set i [FS_getBandIndexFromNoteBookPage $spin]

    if { ! [info exists fs($spin,$i,$what)] } {
	return
    }
    
    if { $fs($spin,$i,$what) } {
	set fs($spin,$i,$what) 0
    } else {
    	set fs($spin,$i,$what) 1
    }
    eval $cmd $i $spin
}


proc FS_getBandIndexFromNoteBookPage {spin} {
    global fs
    set pageName [$fs($spin,nb) raise]
    if { $pageName == "multiband" } {
	set index [expr [lindex $fs($spin,bandlist) end] + 1]
    } else {
	set index [string trimleft $pageName band]
    }
    xcDebug -stderr "FS_getBandIndexFromNoteBookPage:: bandIndex = $index"
    return $index
}

proc FSbind_printTogl {spin} {
    global fs
    set i [FS_getBandIndexFromNoteBookPage $spin]
    printTogl $fs($spin,$i,togl)
}
proc FSbind_SetSurfColor {spin} {
    global fs
    if { [$fs($spin,nb) raise] == "multiband" } {
	return
    }
    set i [FS_getBandIndexFromNoteBookPage $spin]
    FS_SetSurfColor $i $spin
}
proc FSbind_SetCellColor {spin} {
    global fs
    if { [$fs($spin,nb) raise] == "multiband" } {
	return
    }
    set i [FS_getBandIndexFromNoteBookPage $spin]
    FS_SetCellColor $i $spin
}
proc FSbind_glLight {spin} {
    global fs
    set i [FS_getBandIndexFromNoteBookPage $spin]
    glLight $fs($spin,$i,togl)
}
proc FSbind_ModAntiAlias {spin} {
    global fs
    set i [FS_getBandIndexFromNoteBookPage $spin]
    FS_ModAntiAlias $i $spin
}
proc FSbind_ModDepthCuing {spin} {
    global fs
    set i [FS_getBandIndexFromNoteBookPage $spin]
    FS_ModDepthCuing $i $spin
}


proc FSbutton_SmallToolbox {button i spin} {
    global fs

    switch -exact -- $button {
	bz {
	    set fs($spin,$i,text_celltype) "first Brillouin zone"
	    celltype:FS_fsConfig $i $spin
	}
	para {
	    set fs($spin,$i,text_celltype) "reciprocal primitive cell"
	    celltype:FS_fsConfig $i $spin
	}
	nocrop {
	    if { $fs($spin,$i,nocropbz) } {
		set fs($spin,$i,cropbz) 0 
	    } else {
		set fs($spin,$i,cropbz) 1
	    }
	    FS_fsConfig $i $spin
	}
	nocell {
	    set fs($spin,$i,displaycell) 0
	    FS_fsConfig $i $spin
	}
	wirecell {
	    set fs($spin,$i,displaycell)     1
	    set fs($spin,$i,celldisplaytype) wire
	    FS_fsConfig $i $spin
	}
	solidcell {
	    set fs($spin,$i,displaycell)     1
	    set fs($spin,$i,celldisplaytype) solid
	    FS_fsConfig $i $spin
	}
	solidwirecell {
	    set fs($spin,$i,displaycell)     1
	    set fs($spin,$i,celldisplaytype) solidwire
	    FS_fsConfig $i $spin	    
	}
    }
}
