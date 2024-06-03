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
# this proc is called withinh FS_Main. Its purpuse is the display of
# all merged FS bands.
#
proc FS_Multi {nb {spin {}}} {
    global fs xcMisc system
    
    set i [expr [lindex $fs($spin,bandlist) end] + 1]

    $nb insert $i multiband -text "Merged Bands" \
	-createcmd [list FS_RenderMultiSurface $fs($spin,bandlist) $spin]

    #
    # page container frame
    #
    set f    [$nb getframe multiband]
    set togl $f.togl$i

    #
    # toolbox frame
    #
    set ft [frame $f.container -relief raised -bd 1]
    pack $ft -side top -expand 0 -fill x -padx 0m -pady 0m 
    FS_Toolbox $ft $togl $spin $i multiband
    set fs($spin,$i,show_toolbox_frame)       1
    set fs($spin,$i,toolbox_frame)            $ft
    set fs($spin,$i,toolbox_frame_pack)       [pack info $ft]
    set fs($spin,$i,toolbox_frame_packbefore) $togl

    set fs($spin,$i,antialias)      0
    set fs($spin,$i,depthcuing)     0
    
    #
    # Togl
    #
    global toglOpt
    
    set fs($spin,$i,togl) [togl $togl \
			       -ident  $togl \
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

    # take care of togl's background
    FS_UserBackground $togl

    bind $fs($spin,$i,togl) <B1-Motion>        { xc_B1motion %W %x %y }
    bind $fs($spin,$i,togl) <B2-Motion>        { xc_B2motion %W %x %y }
    bind $fs($spin,$i,togl) <B1-ButtonRelease> { xc_Brelease %W B1; MouseZoomBrelease %W }
    bind $fs($spin,$i,togl) <B2-ButtonRelease> { xc_Brelease %W B2 }
    bind $fs($spin,$i,togl) <Button-3>         [list FS_PopupMenu %W %X %Y $i $spin multiband]
    bind $fs($spin,$i,togl) <Shift-B1-Motion>         {  MouseZoom %W %X %Y }
    bind $fs($spin,$i,togl) <Shift-B1-ButtonRelease>  {  MouseZoomBrelease %W }
    bind $fs($spin,$i,togl) <Button-4>  {  MouseWheelZoom %W +}
    bind $fs($spin,$i,togl) <Button-5>  {  MouseWheelZoom %W -}


    # ------------------------------------------------------------------------
    # DISPLAY-FUNCTION
    # ------------------------------------------------------------------------
    cry_dispfuncmultiFS $fs($spin,$i,togl) -togllist $fs($spin,togllist)
}


proc FS_RenderMultiSurface {bandlist {spin {}}} {
    global fs

    foreach i $bandlist {
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
}
