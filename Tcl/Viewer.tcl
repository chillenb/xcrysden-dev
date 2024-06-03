############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/Viewer.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
# ------                                                                    #
#  Modified by Eric Verfaillie ericverfaillie@yahoo.fr EV                   #
#  may 2004                                                                 #
#############################################################################

proc ViewMol {top {file {}} } {
    global light mode2D widsize maingeom fileselect style3D radio xcMisc system \
	dispmode xcFonts

    if { [winfo exists .mesa] } {
	return
    }

    set top .
    wm title $top "XCrySDen"
    wm iconname $top "XCrySDen"
    #wm geometry $top +100+100

    # ---------------
    # INITIALIZATIONS
    # ---------------

    set fac1 $xcMisc(resolution_ratio1)
    set fac2  $xcMisc(resolution_ratio2)
    set widsize(conW) [expr round(130 * $fac2)]
    set widsize(menH) [expr round(35 * $fac1)]
    set widsize(meaH) [expr round(35 * $fac1)]

    # make a global widget arrangement, that is menu, mesawin, control & mea

    PlaceGlobWin 0 [expr round(670 * $fac1)] [expr round(670 * $fac1)]
    #if { [winfo exists .title] } { 
    #wm deiconify .
    #}
    catch {wm deiconify .}
    focus $top
    #tkwait visibility .mesa

    #
    # create Menubar
    #

    ViewMolMenu .menu .mesa

    #
    # bindings
    #

    bind $top <Control-q> { exit_pr }
    bind $top <Control-l> { maximizeDisplay on }
    bind $top <Configure> { ResizeGlobWin }

    bind .mesa <Double-Button-1>  { maximizeDisplay on }
    bind .mesa <B1-Motion>        { bindCB {xc_B1motion .mesa %x %y }}
    bind .mesa <B2-Motion>        { bindCB {xc_B2motion .mesa %x %y }}
    bind .mesa <B1-ButtonRelease> { bindCB {xc_Brelease .mesa B1}; MouseZoomBrelease %W }
    bind .mesa <B2-ButtonRelease> { bindCB {xc_Brelease .mesa B2}; ResetCursor; }
    bind .mesa <Button-3>         { popupMenu %W %X %Y }
    bind .mesa <Shift-B1-Motion>         { MouseZoom %W %X %Y }
    bind .mesa <Shift-B1-ButtonRelease>  { MouseZoomBrelease %W }

    global tcl_platform
    if { $tcl_platform(platform) == "unix" } {
	bind .mesa <Button-4>         { MouseWheelZoom %W + }
	bind .mesa <Button-5>         { MouseWheelZoom %W - }    
    } else {
	bind . <MouseWheel>       { WindowsMouseWheel .mesa %D }    
    }
    

    # ------------------------------------------------------------------------
    #    CONTROL FRAME      
    # ------------------------------------------------------------------------
    # Control Frame will consist from more frames, such as:
    # special frame for DisplayModes, special frame for Translation, Rotation,
    # etc.
    # ------------------------------------------------------------------------

    set c [canvas .ctrl.c -yscrollcommand [list .ctrl.s set]]
    scrollbar .ctrl.s -width 8 -orient vertical -command [list .ctrl.c yview]
    pack .ctrl.s -side right -fill y
    pack .ctrl.c -side left -expand 1

    # here get the requested width of scrollbar
    set widsize(scrW) [winfo reqwidth .ctrl.s]
    set widsize(canW) [expr round(130 * $fac2) - $widsize(scrW) + 1]

    # create FRAME to hold every LABEL&ENTRY
    set fc [frame .ctrl.c.f -bd 0 -width $widsize(canW)]
    $c create window 0 0 -anchor nw -window $fc -width $widsize(canW) \
	    -tags frame

    # define various frames for above reason

    set f1 [frame $fc.fr1 -relief raised -bd 1]
    set f2 [frame $fc.fr2 -relief raised -bd 1 -class ctrl]
    set f3 [frame $fc.fr3 -relief raised -bd 1]

    # miscellaneous images

    image create photo ak -format gif -file $system(BMPDIR)/xcrysden_big.gif
    set fdummy [button $fc.fdummy -image ak -relief ridge -bd 4 \
    	    -anchor center -height 100 -command xcAbout]

    # ------------------------------------------------------------------------
    # FRAME #1 for translation
    # ------------------------------------------------------------------------

    ViewMolFrame1 $f1

    # ------------------------------------------------------------------------
    # FRAME #2 for rotation
    # ------------------------------------------------------------------------

    ViewMolFrame2 $f2

    pack $f1 $f2 -side top -fill x -padx 0 -pady 0 \
	    -ipadx 2 -ipady 2
    pack $f3 -side top -fill both -expand 1 -padx 0 -pady 0 \
	    -ipadx 2 -ipady 2
    pack $fdummy -side bottom -fill both -expand 1
    
    # ------------------------------------------------------------------------
    # FRAME #3 for DISPLAYMODES
    # do we want Lighting (3D or 2D DisplayModes)
    # ------------------------------------------------------------------------
    set f31 [frame $f3.1 -class Radios]
    set f32 [frame $f3.2 -class Radios]
    set light On
    #canvas $f3.canv -width $widsize(conW)
    canvas $f3.canv 
    pack $f31 $f32 -side top -fill x -expand 1
    pack $f3.canv -side top -expand 1

    set dispmode(mode3D_f2_packinfo) [pack info $f32]

    tkwait visibility $f3.canv    

    set rs [RadioButCmd $f31 "Lighting:" light Lighting top left 1 1 5 "On" "Off"]
    foreach r [lrange $rs 1 end] {
    	#$r config -indicatoron 0 -font $xcFonts(small)
	$r config -indicatoron 0
	pack configure $r -padx 1 -pady 3 -ipadx 5 -ipady 1 -fill x -expand 1
    }

    #Lighting On
    #$f3.canv configure -height [winfo reqwidth $f3.canv]

    # ------------------------------------------------------------------------
    #     CONTROL FRAME - END - END - END - END - END - END - END - END
    # ========================================================================


    # ------------------------------------------------------------------------
    #  CLOSE FRAME  
    # ------------------------------------------------------------------------
    set maingeom "Maxi"
    set b0 [button .close.0 -text "F" -width 1 -command maximizeDisplay]
    set b1 [button .close.1 -textvariable maingeom -width 3 -command MainGeom]
    set b2 [button .close.2 -text "Exit" -width 3 -command exit_pr]
    pack $b0 $b1 $b2 -side left -padx 0 -pady 0 -ipadx 0 -ipady 0 -expand 1 -fill x

    set mode2D(WF) On
    set mode2D(PL) Off
    set mode2D(BS1) Off
    set mode2D(BS2) Off

    # ------------------------------------------------------------------------
    # MEASURE FRAME      
    # ------------------------------------------------------------------------

    set meaF [frame .mea.f1 -relief ridge -bd 2 -class mea]
    set atmi [button $meaF.ainf \
	    -highlightthickness 1 -text "Atoms Info" -command \
	    [list PreSel .sel_dist .mesa "Atoms Info" \
	    "For \"Atom Info\" click on the atom" AtomInfo 20]]
    set flbd [button $meaF.dis \
	    -highlightthickness 1 -text "Distance" -command \
	    [list PreSel .sel_dist .mesa "Distance" \
	    "For \"Distance\" click on two atoms" Distance 2]]
    set flba [button $meaF.ang \
	    -highlightthickness 1 -text "Angle" -command \
	    [list PreSel .sel_dist .mesa "Angle" \
	    "For \"Angle\" click on three atoms" Angle 3]]
    set flbdi [button $meaF.dih \
	    -highlightthickness 1 -text "Dihedral" -command \
	    [list PreSel .sel_dist .mesa "Dihedral" \
	    "For \"Dihedral\" click on four atoms" Dihedral 4]]
    #set flbp [button $meaF.pln \
    #	    -highlightthickness 0 -text "Plane" -command \
    #	    [list PreSel .sel_plane .mesa "Plane" \
    #	    "For \"Plane\" click on three atoms" Plane 3]]
    #$flbp config -state disabled

    set style3D(draw)     "solid"
    set style3D(shade)    "smooth"
    set style3D(anaglyph) 0
    set style3D(stereo)   1
    set style3D(point)    "1"
    set modl [frame .mea.f2 -relief raised -bd 1]

    set r_solid  [radiobutton $modl.solid -image dm_solid -highlightthickness 1 \
		      -variable style3D(draw) -value solid -indicatoron 0 \
		      -selectcolor \#ff2222 -anchor center \
		      -command [list Style3D draw solid]]
    set r_wire   [radiobutton $modl.wire -image dm_wire -highlightthickness 1 \
		      -variable style3D(draw) -value wire -indicatoron 0 \
		      -selectcolor \#ff2222 -anchor center \
		      -command [list Style3D draw wire]]

    set r_anaglyph   [checkbutton $modl.anaglyph -image dm_anaglyph -highlightthickness 1 \
		      -variable style3D(anaglyph) -indicatoron 0 \
		      -selectcolor \#4444ff -anchor center \
		      -command [list Style3D draw anaglyph]]     
    set r_stereo  [checkbutton $modl.stereo -image dm_stereo -highlightthickness 1 \
		       -variable style3D(stereo) -indicatoron 0 \
		       -selectcolor \#4444ff -anchor center \
		       -command [list Style3D draw stereo]] 
    
    set separator_1 [frame $modl.s1 -width 2 -relief raised -bd 1]
    set separator_2 [frame $modl.s2 -width 2 -relief raised -bd 1]
    set separator_3 [frame $modl.s3 -width 2 -relief raised -bd 1]
    
    set r_smooth [radiobutton $modl.smooth -image dm_smooth -highlightthickness 1 \
		      -variable style3D(shade) -value smooth -indicatoron 0 \
		      -selectcolor \#44ff44 -anchor center \
		      -command [list Style3D shade smooth]]
    set r_flat   [radiobutton $modl.flat -image dm_flat -highlightthickness 1 \
		      -variable style3D(shade) -value flat -indicatoron 0 \
		      -selectcolor \#44ff44 -anchor center \
		      -command [list Style3D shade flat]]

    set rep_unit [radiobutton $modl.unit -image rep_unit -highlightthickness 1 \
		      -variable radio(unitrep) -value "cell" -indicatoron 0 \
		      -selectcolor \#4444ff -anchor center \
		      -command [list CellMode 1]]
    set rep_asym [radiobutton $modl.asym -image rep_asym -highlightthickness 1 \
		      -variable radio(unitrep) -value "asym" -indicatoron 0 \
		      -selectcolor \#4444ff -anchor center \
		      -command [list CellMode 1]]

    # reload button
    set separator_4 [frame $modl.s4 -width 2 -relief raised -bd 1]
    set b_reload [button $modl.reload -text "Reload" -command reloadCmd]
    
    set b_pack_option {-side left -fill y -padx 0 -pady 0 -ipadx 0 -ipady 0}
    set s_pack_option {-side left -fill y -padx 2 -pady 0 -ipadx 0 -ipady 0}
	    
    eval pack $separator_1 $s_pack_option
    eval pack $rep_unit $rep_asym $b_pack_option
    eval pack $separator_2 $s_pack_option
    eval pack $r_solid $r_wire $r_anaglyph $b_pack_option ; #$r_stereo 
    eval pack $separator_3 $s_pack_option
    eval pack $r_smooth $r_flat $b_pack_option
    
    eval pack $separator_4 $s_pack_option
    eval pack $b_reload $b_pack_option -pady 2
    reloadButton $b_reload
    
    global xcFonts
    foreach {wid text} {
	solid    "solid molecular draw-style"
	wire     "wire molecular draw-style"
	anaglyph "anaglyph molecular draw-style" 
	stereo   "stereo display mode"
	smooth   "smooth molecular shading"
	flat     "flat molecular shading"
	unit     "nicely-cut unit cell"
	asym     "translational asymmetric unit"
	xsf2     "update XSF2 structure"
    } {
	set path $modl.$wid
	if { [winfo exists $path] } {
	    DynamicHelp::register $path balloon $text
	}
    }

    pack $meaF -side left -fill y
    pack $modl -ipadx 5 -expand 1 -fill both -side left
    pack $atmi $flbd $flba $flbdi -side left -padx 0 -ipadx 0 -expand 1 -fill y

    #pack $atmi $flbd $flba $flbdi \
    #	    $styl $sty $shdl $shd -side left -padx 0 

    # before opening a structure, a xc_mesacontext must be made !!!!!    
    # xc_mesacontext .mesa

    #--------------------------------------#
    # configure control canvas & scrollbar #
    #--------------------------------------#
    update
    set width  [expr [winfo width .ctrl] - $widsize(scrW)]
    set height [winfo height .ctrl]
    set h      [winfo reqheight $fc]
    $c config -width $width -height $height \
	    -scrollregion "0 0 $width $h"     
    #----------------------------------------------------#
    # END -- configure control canvas & scrollbar -- END #
    #----------------------------------------------------#

    if { $file != {} } {
	# open file and display it !!!
	OpenXSFStruct .mesa $file
    }

    # disable all neccessary widgets
    xcUpdateState

    #
    # delete WELCOME window
    #
    destroyWelcome

    # load the attributes from the definition file ...
    defLoadAttributes

    # scrolling of the control toolbox
    scrollControlToolboxCmd {.ctrl.c .ctrl.s}

    global exit_viewer_win
    set exit_viewer_win .
    bind . <Destroy> { exit_viewer . }
}

proc exit_viewer {w} {        
    global exit_viewer_win

    if { ! [winfo exists $w] } {
	bind $w <Destroy> {}
	clean_exit
    }
}


proc scrollControlToolboxCmd {wlist} {
    
    set B4    [list .ctrl.c yview scroll -2 units]
    set B5    [list .ctrl.c yview scroll +2 units]
    # TODO: please tune the %D on windows
    set Wheel [list .ctrl.c yview scroll %D units]


    foreach w $wlist {
	mouseWheelScroll $w $B4 $B5 $Wheel
    }
}


proc ResizeGlobWin {} {
    global w h
    if { ! [info exists w] } { set w 0 }
    if { ! [info exists h] } { set h 0 }
    # update only if size of "toplevel ." has changed
    if { $w != [winfo width .] || $h != [winfo height  .] } {
	set w [winfo width .]
	set h [winfo height  .]
	PlaceGlobWin 1 $w $h
	update
	ConfigControlWid
    }
}


proc PlaceGlobWin {type w h} {
    global widsize unmapWin

    # type == 0 --> initial arrangement
    # type == 1 --> <Configure> rearrangement
    # type == display_resize --> display-window was requested to be resized
    # type == fullscreen-off --> fullscreen-mode was turned off

    if { $type == 0 } {
	# this is for toplevel width & height     
	set w [expr $w + $widsize(conW)]
	set h [expr $h + $widsize(menH) + $widsize(meaH)]
	. configure -width $w -height $h
	#tkwait visibility .    
    }
    
    # on TOP will be MENU; on left xc_mesawin -> that is a custom widget
    # made for XCrySDen;
    # on right will be frame, that will hold everything --> this frame's name
    # will be control !!!
    # On BOTTOM there will be a MEASURE frame -> this frame's name will
    # be "mea"
    # On RIGHT BOTTOM there will be another frame that will hold Close button 

    # some calculation for PLACEing the three widged mentioned above
    # frame 
    
    if { $type != "display_resize" } {
	set controlW $widsize(conW)
	set controlH [expr $h - $widsize(meaH)]
	set controlX [expr $w - $controlW]
	set controlY 0
	set widsize(conH) $controlH

	set menuX 0
	set menuY 0
	set menuW [expr $w - $controlW]
	#set menuW $w
	set menuH $widsize(menH)
	
	#new
	#set controlH [expr $controlH - $menuH]
	#set controlY $menuH
	#/
	
	set meaH $widsize(meaH)
	set meaW $menuW
	
	set mesaX 0
	set mesaY $menuH
	set mesaW [expr $w - $controlW]
	set mesaH [expr $h - $menuH - $meaH]
	
	set meaX 0
	set meaY [expr $menuH + $mesaH]
	
	set closeX [expr $w - $controlW]
	set closeY $meaY
	set closeW $controlW
	set closeH $meaH
    } else {
	set mesaX 0
	set mesaY $widsize(menH)
	set mesaW $w
	set mesaH $h
    
	set controlW $widsize(conW)
	set controlH $h
	set controlX $w
	set controlY 0
	set widsize(conH) $controlH

	set menuX 0
	set menuY 0
	set menuW $w
	set menuH $widsize(menH)
	
	set meaH $widsize(meaH)
	set meaW $menuW
		
	set meaX 0
	set meaY [expr $menuH + $mesaH]
	
	set closeX [expr $w - $controlW]
	set closeY $meaY
	set closeW $controlW
	set closeH $meaH

	set topW [expr $mesaW + $controlW]
	set topH [expr $menuH + $mesaH + $meaH]
	. config -width $topW -height $topH
    }
        
    # initial arrangement of GlobWid
    if { $type == 0 } {
	# definition of above widgets
	set menuWid  [frame .menu -relief raised -bd 1 -class MenuBar]
	#set mesaWid  [xc_mesawin .mesa -bg "#000"]
	
	global tcl_platform stereo_visual
	set stereo true
	if { ( $tcl_platform(platform) eq "windows" ) || ($tcl_platform(os) eq "Darwin" ) } {
            set stereo false
        }
        global toglOpt
        # N.B. -stereo is set to none at the moment        
	set mesaWid [togl .mesa \
			 -width          400 \
			 -height         400 \
			 -ident          .mesa \
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
                         -create  togl_create \
                         -display togl_display \
                         -reshape togl_reshape \
                         -destroy togl_destroy \
                         -timer   togl_timer ]

        puts stderr "Viewer: mesaWid = $mesaWid"
	
	set conWid   [frame .ctrl  -relief flat   -bd 0 -class Viewer]
	set meaWid   [frame .mea   -relief flat   -bd 0 -class Viewer]
	set closeWid [frame .close -relief raised -bd 1 -class Viewer]

	set unmapWin(frame,main) [frame .mesa.unmap -bg "#000"]

	# place the above widgets
	place $menuWid -x $menuX -y $menuY \
		-width $menuW -height $menuH -anchor nw
	place $mesaWid -x $mesaX -y $mesaY \
		-width $mesaW -height $mesaH -anchor nw
	place $conWid -x $controlX -y $controlY \
		-width $controlW -height $controlH -anchor nw
	place $meaWid -x $meaX -y $meaY \
		-width $meaW -height $meaH -anchor nw
	place $closeWid -x $closeX -y $closeY \
		-width $closeW -height $closeH -anchor nw

	pack $unmapWin(frame,main) -anchor sw -side bottom
	set unmapWin(packinfo,$unmapWin(frame,main)) \
		[pack info $unmapWin(frame,main)]

	after idle {
	    # check if stereo is supported
	    global stereo_visual
	    
	    #set stereo [xc_stereo]	
	    set stereo [lindex [.mesa configure -stereo] end]
	    if { $stereo != "none" } {
		#puts stderr "\n*** the hardware supports the stereo ***\n"
		set stereo_visual true
	    } else {
		#puts stderr "\n*** the hardware does not support the stereo ***\n"
		set stereo_visual false
	    }
	}
    } elseif { $type == "fullscreen-off" } {
	place .menu -x $menuX -y $menuY \
		-width $menuW -height $menuH -anchor nw
	place .mesa -x $mesaX -y $mesaY \
		-width $mesaW -height $mesaH -anchor nw
	place .ctrl -x $controlX -y $controlY \
		-width $controlW -height $controlH -anchor nw
	place .mea -x $meaX -y $meaY \
		-width $meaW -height $meaH -anchor nw
	place .close -x $closeX -y $closeY \
		-width $closeW -height $closeH -anchor nw
    } else {
	# configure above widgets
	place configure .menu -x $menuX -y $menuY \
		-width $menuW -height $menuH -anchor nw
	place configure .mesa -x $mesaX -y $mesaY \
		-width $mesaW -height $mesaH -anchor nw
	place configure .ctrl -x $controlX -y $controlY \
		-width $controlW -height $controlH -anchor nw
	place configure .mea -x $meaX -y $meaY \
		-width $meaW -height $meaH -anchor nw
	place configure .close -x $closeX -y $closeY \
		-width $closeW -height $closeH -anchor nw

	#place configure $unmapWin(frame,main) -width 100 -height 100 \
	#	-x 0 -y $mesaH -anchor sw
    }

    return
}


proc ViewMolFrame1 {f} {
    global translationStep system xcMisc widsize mody viewer

    set fac1 $xcMisc(resolution_ratio1)
    set fac2 $xcMisc(resolution_ratio2)
 
    #proc _translationStepChanged {name1 name2 op} {
    #	global translationStep
    #	# this is a way-around of uncorrected BUG
    #	puts stderr "----->translationStep: $translationStep"
    #	#set a
    #	if { [string match $translationStep "nan"] } {
    #	    set translationStep 0.05
    #	}
    #	#puts stderr "translationStep: $translationStep"
    #}
    ##initial value of translationStep
    #
    ##set translationStep 0.05
    #trace variable translationStep rw _translationStepChanged

    set translationStep 0.05

    # make a frame that will hold UP, DOWN, LEFT & RIGHT button
    set w [expr round(122 * $fac2) - $widsize(scrW)]
    set h [expr round(122 * $fac2) - $widsize(scrW)]
    # set w [expr round(122 * $fac2)]
    # set h [expr round(122 * $fac2)]
    set cpad [expr round(4 * $fac1)]
    set ff [frame $f.1 -relief groove -bd 2 -width $w -height $h]
    pack $ff -side top -expand 1 -padx $cpad -pady $cpad 
    
    # add something to $b bacause of button's borders 
    set b [expr round(26 * $fac2) + 8];
    set b1 [button $ff.b1 -image up -anchor center \
	    -anchor center -width $b -height $b]
    set b2 [button $ff.b2 -image down -anchor center \
	    -anchor center -width $b -height $b]
    set b3 [button $ff.b3 -image left -anchor center \
	    -anchor center -width $b -height $b]
    set b4 [button $ff.b4 -image right -anchor center \
	    -anchor center -width $b -height $b]
    set b5 [button $ff.b5 -image center -anchor center \
	    -anchor center -width $b -height $b]
    
    # now place the buttons
    set fw [expr $w / 2 - round(4 * $fac1)]
    set fh [expr $h / 2 - round(2 * $fac1)]
    set offset 0

    place $b1 -x $fw -y [expr $fh - ( $b + $offset )] \
	    -anchor center -width $b -height $b
    place $b2 -x $fw -y [expr $fh + ( $b + $offset )] \
	    -anchor center -width $b -height $b
    place $b3 -x [expr $fw - ( $b + $offset )] -y $fh \
	    -anchor center -width $b -height $b
    place $b4 -x [expr $fw + ( $b + $offset )] -y $fh \
	    -anchor center -width $b -height $b
    place $b5 -x $fw -y $fh -anchor center -width $b -height $b
        
    foreach {wid text} {
	b1    "translate up"
	b2    "translate down"
	b3    "translate left"
	b4    "translate right"
	b5    "translate back to center"
    } {
	set path $ff.$wid		
	DynamicHelp::register $path balloon $text
    }

        
    ###############
    # Zoom + Zoom -
    global tcl_platform
    if { $tcl_platform(platform) != "unix" } {
	set _w 6
    } else {
	set _w 4
    }
    set f2 [frame $f.2 -relief flat -bd 0]
    set zoom1 [button $f2.zoom1 -text "Zoom +" -width $_w ]
    set zoom2 [button $f2.zoom2 -text "Zoom -" -width $_w ]
    
    # define the bindings

    set viewer(bind) ""
    foreach {path name script} {
	b1 up     {Transl %W +y}
	b2 down   {Transl %W -y}
	b3 left   {Transl %W -x}
	b4 right  {Transl %W +x}
	b5 center {
	    bindCB [list xc_newvalue .mesa $mody(L_TR_XTRANSL) 0.0]
	    bindCB [list xc_newvalue .mesa $mody(L_TR_YTRANSL) 0.0]
	}
    } {	
	lappend viewer(bind) $name
	set viewer(bind,path,$name)   $ff.$path
	set viewer(bind,script,$name) $script
    }
    foreach {path name script} {
	zoom1 zoom1 {Transl %W +z}
	zoom2 zoom2 {Transl %W -z}
    } {
	lappend viewer(bind) $name
	set viewer(bind,path,$name)   $f2.$path
	set viewer(bind,script,$name) $script
    }

    ######################
    # SCALE for translationStep
    set f3 [frame $f.3 -bd 2 -relief groove]
    set sc [scale $f3.scale111 -from 0 -to 0.55 \
	    -length [expr round(115 * $fac2) - $widsize(scrW)] \
	    -variable translationStep \
	    -orient horizontal -tickinterval 0.25 \
	    -digits 2 -resolution 0.05 -showvalue true \
	    -width [expr round(10 * $fac2)]]
    set trlab   [label $f3.lab -text "Translation\nStep:" -justify right]
    set trentry [entry $f3.entry -relief sunken \
		     -textvariable translationStep \
		     -width [expr round(5 * $fac2)] \
                     -validate key -validatecommand {doubleOnly %P}]
    pack $sc -side bottom
    pack $trlab $trentry -side left -expand 1 -padx 0 -pady 0
    pack $zoom1 $zoom2 -side left -expand 1 
    pack $f2 -side top -fill x -expand 1 -pady 1

    set i5 [expr round(5 * $fac2)]
    set i15 [expr round(15 * $fac2)]
    pack $f3 -side bottom  -expand 1 -pady $i5 -padx $i5 -ipadx $i15 -ipady $i5

    # check if we need to adjust the width of .ctrl frame
    set w_zoom [expr [winfo reqwidth $zoom1] * 2]
    if { $w_zoom > $w } {
	set widsize(conW) [expr $w_zoom + $widsize(scrW) + 4]
	set w [winfo width .]
	set h [winfo height  .]
	PlaceGlobWin 1 $w $h
	.ctrl.c itemconfigure frame -width [expr $widsize(conW) - $widsize(scrW) + 1]
    }     
}


proc ViewMolFrame2 {f} {
    global rotstep downB1 xcMisc system xcFonts xcColors viewer

    #initial value of rotstep
    set rotstep 10

    set fac2 $xcMisc(resolution_ratio2)
    global tcl_platform
    if { $tcl_platform(platform) != "unix" } {
	set w [expr round(6 * $fac2)]
    } else {
	set w [expr round(4 * $fac2)]
    }
    set f1 [frame $f.1 -relief flat -bd 0]
    set rx1 [button $f1.rx1 -text "Rot +X" -width $w -highlightthickness 1]
    set rx2 [button $f1.rx2 -text "Rot -X" -width $w -highlightthickness 1]
    bind $rx1 <ButtonPress-1>   {Rotate %W +x}
    bind $rx2 <ButtonPress-1>   {Rotate %W -x}
    bind $rx1 <ButtonRelease-1> RelB1
    bind $rx2 <ButtonRelease-1> RelB1

    set f2 [frame $f.2 -relief flat -bd 0]
    set ry1 [button $f2.ry1 -text "Rot +Y" -width $w -highlightthickness 1] 
    set ry2 [button $f2.ry2 -text "Rot -Y" -width $w -highlightthickness 1] 
    bind $ry1 <ButtonPress-1>   {Rotate %W +y}
    bind $ry2 <ButtonPress-1>   {Rotate %W -y}
    bind $ry1 <ButtonRelease-1> RelB1
    bind $ry2 <ButtonRelease-1> RelB1

    set f3 [frame $f.3 -relief flat -bd 0]
    set rz1 [button $f3.rz1 -text "Rot +Z" -width $w -highlightthickness 1] 
    set rz2 [button $f3.rz2 -text "Rot -Z" -width $w -highlightthickness 1] 
    bind $rz1 <ButtonPress-1>   {Rotate %W +z}
    bind $rz2 <ButtonPress-1>   {Rotate %W -z}
    bind $rz1 <ButtonRelease-1> RelB1
    bind $rz2 <ButtonRelease-1> RelB1

    # t.k
    foreach {name script} {
	rx1 {Rotate %W +x}
	rx2 {Rotate %W -x}
    } {	
	lappend viewer(bind) $name
	set viewer(bind,path,$name)   $f1.$name
	set viewer(bind,script,$name) $script
    }
    foreach {name script} {
	ry1 {Rotate %W +y}
	ry2 {Rotate %W -y}
    } {	
	lappend viewer(bind) $name
	set viewer(bind,path,$name)   $f2.$name
	set viewer(bind,script,$name) $script
    }
    foreach {name script} {
	rz1 {Rotate %W +z}
	rz2 {Rotate %W -z}
    } {	
	lappend viewer(bind) $name
	set viewer(bind,path,$name)   $f3.$name
	set viewer(bind,script,$name) $script
    }
    #/

    set f4 [frame $f.4 -bd 2 -relief groove]
    set sc [scale $f.4.scale -from -180 -to 180 \
	    -length [expr round(115 * $fac2)] -variable rotstep \
	    -orient horizontal  -tickinterval 180 \
	    -digits 3 -resolution 1 -showvalue true \
	    -width [expr round(10 * $fac2)]]
    set rotlab [label $f4.lab -text "Rotation\nStep:" -justify right]
    set rotentry [entry $f4.entry -relief sunken \
                      -textvariable rotstep \
                      -width [expr round(5 * $fac2)] \
                      -validate key -validatecommand {doubleOnly %P}]
    pack $sc -side bottom 
    pack $rotlab $rotentry -side left -expand 1 -padx 0 -pady 0
    pack $f1 $f2 $f3 -side top -fill x -expand 1 -pady 1
    set i5 [expr round(5 * $fac2)]
    set i15 [expr round(15 * $fac2)]

    set f5 [frame $f.5 -bd 0 -relief flat]
    set b  [expr round(26 * $fac2)]
    foreach {row col bt im com text} [list \
	    0 0 $f5.xy rotXY {xc_rotate .mesa xy 0.0} "orient to XY plane" \
	    0 1 $f5.yz rotYZ {xc_rotate .mesa yz 0.0} "orient to YZ plane" \
	    0 2 $f5.xz rotXZ {xc_rotate .mesa xz 0.0} "orient to XZ plane" \
	    1 0 $f5.ac rotAC {xc_rotate .mesa ac 0.0} "orient to AC plane" \
	    1 1 $f5.ab rotAB {xc_rotate .mesa ab 0.0} "orient to AB plane" \
	    1 2 $f5.bc rotBC {xc_rotate .mesa bc 0.0} "orient to BC plane" ] {
	button $bt -image $im -command $com -width $b -height $b -anchor center
	grid $bt -row $row -column $col -padx 1 -pady 1
	DynamicHelp::register $bt balloon $text
    }
    pack $rx1 $rx2 $ry1 $ry2 $rz1 $rz2 -side left  \
	    -expand 1 -padx 0 -pady 0 -ipadx 0 -ipady 0
    pack $f5 -side bottom  -expand 1 -pady 2 -padx 1
    pack $f4 -side bottom  -expand 1 -pady $i5 -padx $i5 \
	    -ipadx $i15 -ipady $i5

    # discrete/continuous checkbutton
    set viewer(rot_zoom_button_mode) "Discrete"

    set rot_zoom_mode_w [RadioButCmd $f "Rotation+zoom\nbuttons mode:" \
			     viewer(rot_zoom_button_mode) Viewer:rotZoomButtonMode \
			     top top 1 1 5 \
			     "Discrete" "Click-and-hold" "Click-and-click"]
    set viewer(rot_zoom_button_mode_rbs) [lrange $rot_zoom_mode_w 1 end]
    foreach r $viewer(rot_zoom_button_mode_rbs) {
	$r config -indicatoron 0
    	#$r config -indicatoron 0 -font $xcFonts(small)
    	pack configure $r -padx 1 -pady 1 -ipadx 5 -ipady 1 -fill x -expand 1
    }
    Viewer:rotZoomButtonMode
}


proc Rotate {b dir} {
    global rotstep B1down err viewer

    # check if $rotstep is real number
    number rotstep real
    if { $err == 1 } { 
	return
    }
    
    if { $viewer(rot_zoom_button_mode) == "Click-and-hold" } {
	$b configure -relief sunken
	set B1down 1
	
	while { $B1down } {
	    xc_rotate .mesa $dir $rotstep
	    update 	
	}        

	$b configure -relief raised	
	return -code break

    } elseif { $viewer(rot_zoom_button_mode) == "Discrete" } {
	xc_rotate .mesa $dir $rotstep
	update 

    } elseif { $viewer(rot_zoom_button_mode) == "Click-and-click" } {
	if { ! [info exists viewer($b,$dir)] } {
	    set viewer($b,$dir) 0
	}
	if { ! $viewer($b,$dir) } {
	    Viewer:_clickAndClick start $b $dir
	    
	    set B1down 1	    
	    while { $B1down } {
		xc_rotate .mesa $dir $rotstep
		update 	
	    }        

	    Viewer:_clickAndClick end $b $dir
	    return -code break
	} else {
	    set viewer($b,$dir) 0
	    RelB1
	}
    }
}

proc Transl {b dir} {
    global translationStep B1down err viewer

    # check if $rotstep is real number
    number translationStep real
    if { $err == 1 } { 
	return
    }

    if { $viewer(rot_zoom_button_mode) == "Click-and-hold" } {
	$b configure -relief sunken
	set B1down 1
	
	while { $B1down } {
	    xc_translate .mesa $dir $translationStep
	    update 	
	}

	$b configure -relief raised
	return -code break

    } elseif { $viewer(rot_zoom_button_mode) == "Discrete" } {
	xc_translate .mesa $dir $translationStep
	update

    } elseif { $viewer(rot_zoom_button_mode) == "Click-and-click" } {
	if { ! [info exists viewer($b,$dir)] } {
	    set viewer($b,$dir) 0
	}
	if { ! $viewer($b,$dir) } {
	    Viewer:_clickAndClick start $b $dir
	    
	    set B1down 1	    
	    while { $B1down } {
		xc_translate .mesa $dir $translationStep
		update 	
	    }        
	    
	    Viewer:_clickAndClick end $b $dir
	    return -code break
	} else {
	    set viewer($b,$dir) 0
	    RelB1
	}
    }
}


proc Viewer:_clickAndClick {mode b dir} {
    global viewer
    if { $mode == "start" } {
	set viewer($b,$dir) 1
	foreach name $viewer(bind) {
	    $viewer(bind,path,$name) configure -state disabled
	}
	$b configure -relief sunken -state normal
	foreach rb $viewer(rot_zoom_button_mode_rbs) {
	    $rb config -state disabled
	}
    } elseif { $mode == "end" } {
	foreach name $viewer(bind) {
	    $viewer(bind,path,$name) configure -state normal -relief raised
	}	    
	foreach rb $viewer(rot_zoom_button_mode_rbs) {
	    $rb config -state normal
	}
    }
}


proc MouseZoom {W x y} {
    global mouseZoom

    # register Shift-B1-Motion event
    xc_ShiftB1motion $W 

    if { ! [info exists mouseZoom($W,oldX)] } {
	set mouseZoom($W,oldX) $x
	set mouseZoom($W,oldY) $y
    }

    #set w  [winfo width $W]
    #set fx [expr {double($mouseZoom($W,oldX) - $x)/double($w/2.0)}]
    #set f [expr {abs($fx) > abs($fy) ? $fx : $fy}]

    set h    [winfo height $W]
    set fy   [expr {double($mouseZoom($W,oldY) - $y)/double($h/2.0)}]    
    set sign [expr { $fy > 0 ? "+" : "-"}]

    cry_toglzoom $W $fy

    set mouseZoom($W,oldX) $x
    set mouseZoom($W,oldY) $y
}
proc MouseZoomBrelease {W} {
    global mouseZoom
    # unregister Shift-B1-Motion event
    xc_Brelease $W Shift-B1

    if { [info exists mouseZoom($W,oldX)] } { unset mouseZoom($W,oldX) }
    if { [info exists mouseZoom($W,oldY)] } { unset mouseZoom($W,oldY) }
}


proc MouseWheelZoom {W dir} {

    set fy 0.075
    cry_toglzoom $W ${dir}${fy}
}


proc WindowsMouseWheel {win step} {
    if { $step > 0 } {
	MouseWheelZoom $win +
    } elseif { $step < 0 } {
	MouseWheelZoom $win -
    }
}


proc RelB1 {} {
    global B1down
    
    set B1down 0
}


#
# arguments:
#      F  ... fullscreen
proc MainGeom { {F {}} } {
    global oldgeom maingeom global xcMisc

    if { ![info exists oldgeom] } { set oldgeom 0 }

    # go to Maxi
    if { $maingeom == "Maxi" } {
	set maingeom "Mini"
	set oldgeom [winfo geometry .]
	set w       [winfo screenwidth  .]
	set h       [winfo screenheight .]

	###############################
	# testing "FULLSCREEN" option #
	###############################

	if { $F != {} } {
	    wm withdraw .
	    wm overrideredirect . true
	    update
	    wm deiconify .
	}

	###############################
	#wm geometry . ${w}x${h}+$xcMisc(wm_rootXshift)+$xcMisc(wm_rootYshift)
	update
	wm attributes . -fullscreen 1
    } else {
	set maingeom "Maxi"
	###############################
	# testing "FULLSCREEN" option #
	###############################
	if { $F != {} } {
	    wm withdraw .
	    wm overrideredirect . false
	    update
	    wm deiconify .
	}
	###############################
	wm geometry . $oldgeom
	wm attributes . -fullscreen 0
	wm withdraw .
	wm overrideredirect . false
	update
	wm deiconify .
	wm geometry . $oldgeom
    }
}


proc maximizeDisplay {{action on}} {
    global widsize unmapWin maximizeDisplay

    switch -exact -- $action {
	on {
	    # nonFS_geometry == non full-screen geometry

	    set maximizeDisplay(nonFS_geometry) [wm geometry .]
	    bind . <Configure> {}
	    
	    set w [winfo screenwidth  .]
	    set h [winfo screenheight .]
	    wm geometry . ${w}x${h}+0+0
	    update 

	    wm withdraw .
	    wm overrideredirect . true	    
	    wm deiconify .
	    wm geometry . ${w}x${h}+0+0
	    	    
	    place configure .mesa -x 0 -y 0 -width $w -height $h -anchor nw	    
	    
	    foreach wid {.menu .ctrl .mea .close} {
		place forget $wid
	    }
	    place configure .mesa -x 0 -y 0 -width $w -height $h -anchor nw
	    
	    set maximizeDisplay(.mesa_bind_script) [bind .mesa <Button-3>]
            
	    bind .mesa <Button-3> { maximizeDisplay:popupMenu %W %X %Y }
	    bind .mesa <Double-Button-1> { 
		bindCB {xc_Brelease .mesa B1}
		maximizeDisplay off 
	    }
	    bind . <Control-l> { maximizeDisplay off }
	}

	off {

	    if { ![info exists maximizeDisplay(nonFS_geometry)] } {
		return
	    }

	    wm geometry . $maximizeDisplay(nonFS_geometry) 
	    wm withdraw .
	    wm overrideredirect . false
	    update
	    wm deiconify .
	    wm geometry . $maximizeDisplay(nonFS_geometry) 

	    PlaceGlobWin fullscreen-off [winfo width .] [winfo height .]

	    bind .mesa <Button-3> "bindCB [list $maximizeDisplay(.mesa_bind_script)]"
	    bind .mesa <Double-Button-1> { 
		bindCB {xc_Brelease .mesa B1}
		maximizeDisplay on 
	    }
	    bind . <Configure> {ResizeGlobWin}
	    bind . <Control-l> {maximizeDisplay on}
	}
	default {
	    ErrorDialog "wrong action $action for maximizeDisplay, should be on or off"
	}
    }
}


proc maximizeDisplay:popupMenu {W x y} {

    if { ! [winfo exists $W] } {
	return
    }

    if { [winfo exists $W.menu] } {
	destroy $W.menu
    }
    set m [menu $W.menu -tearoff 0]
    tk_popup $m $x $y

    $m add command -label "Exit from full-screen mode" \
	-command {maximizeDisplay off}
    $m add separator
    $m add command -label "Zoom" -command [list toglZoom "Zoom" .mesa]
    $m add cascade -label "Lighting Off Display Modes" -menu $m.light_off
    $m add cascade -label "Lighting On  Display Modes" -menu $m.light_on
    
    set moff [menu $m.light_off -tearoff 0]
    set mon  [menu $m.light_on  -tearoff 0]

    $moff add command -label "WireFrame"   -command [list Display2D WF]
    $moff add command -label "PointLine"   -command [list Display2D PL]
    $moff add command -label "Pipe&Ball"   -command [list Display2D PB]
    $moff add command -label "BallStick-1" -command [list Display2D BS1]
    $moff add command -label "BallStick-2" -command [list Display2D BS2]
    $moff add command -label "SpaceFill"   -command [list Display2D SF]

    $mon add command -label "Stick"     -command [list DisplayOver3D S]
    $mon add command -label "Pipe&Ball" -command [list DisplayOver3D PB]
    $mon add command -label "BallStick" -command [list DisplayOver3D BS]
    $mon add command -label "SpaceFill" -command [list DisplayOver3D SF]

    $m add separator
    
    popupMenu:popup $m .mesa
}



proc Style3D {what how} {
    global style3D

    if {$what == "draw"} {
	if { $how == "stereo" } {
	    global stereo_visual
	    if { $stereo_visual == "false" } {
		WarningDialog "Stereo mode is not supported by your computer"		
	    } 
	} elseif { $how != "anaglyph" } {
	    set style3D(draw) $how
	}
	xc_drawStyle3D .mesa $how 
    }    
    if {$what == "shade"} {
	set style3D(shade) $how
	xc_shadeModel3D .mesa $how
    }
    if {$what == "point"} {
	set style3D(point) $how
	puts stdout "style3D(point)"
	xc_pointSize .mesa $how
    }
}


proc bindCB {cmd} {
    eval $cmd
}


proc Viewer:rotZoomButtonMode {{value {}}} {
    global viewer

    #set viewer(bind,path,$name)   $ff.$path
    #set viewer(bind,script,$name) $script

    # "Discrete" "Click-and-hold" "Click-and-click"
    
    if { $viewer(rot_zoom_button_mode) == "Click-and-hold" } {
	foreach name $viewer(bind) {
	    bind $viewer(bind,path,$name) <ButtonPress-1>   $viewer(bind,script,$name)
	    bind $viewer(bind,path,$name) <ButtonRelease-1> RelB1
	}
    } elseif { $viewer(rot_zoom_button_mode) == "Discrete" } {
	foreach name $viewer(bind) {
	    bind $viewer(bind,path,$name) <ButtonPress-1>   {}
	    bind $viewer(bind,path,$name) <ButtonRelease-1> $viewer(bind,script,$name)
	}
    } elseif { $viewer(rot_zoom_button_mode) == "Click-and-click" } {
	foreach name $viewer(bind) {
	    bind $viewer(bind,path,$name) <ButtonPress-1>   {}
	    bind $viewer(bind,path,$name) <ButtonRelease-1> $viewer(bind,script,$name)
	}
    }
}

# used for entry validation; it does not allow empty string
proc doubleOnly {str} {
    if { $str == {} } {
        return false
    }
    return [string is double $str]
}
