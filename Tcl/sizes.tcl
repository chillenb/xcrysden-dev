#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/sizes.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc ProbeResolution {} {
    global xcMisc xcFonts

    set w [winfo screenwidth  .]
    set h [winfo screenheight .]

    # Debugging: testing different screen resolutions:
    #-------------------------------------------------
    #set w 800
    #set h 600

    set xcMisc(resolution)  ${w}x${h}
    set xcMisc(resolutionX) $w
    set xcMisc(resolutionY) $h
    set xcMisc(resolution_ratio1) 1.0
    set xcMisc(resolution_ratio2) 1.0
    
    set xcFonts(normal_size)      [font actual TkDefaultFont -size]
    set xcFonts(small_size)       [expr $xcFonts(normal_size) - 2]
    set xcFonts(big_size)         [expr $xcFonts(normal_size) + 2]

    set xcMisc(resolution_ratio1) 0.65
    set xcMisc(resolution_ratio2) 0.75
    if { $xcMisc(resolutionX) > 700 } {
	set xcMisc(resolution_ratio1) 0.75
	set xcMisc(resolution_ratio2) 0.85
    } 
    if { $xcMisc(resolutionX) > 1000 } {
        set xcMisc(resolution_ratio1) 0.9
        set xcMisc(resolution_ratio2) 1.0
    }
    if { $xcMisc(resolutionX) > 1200 } {
	set xcMisc(resolution_ratio1) 1.0
	set xcMisc(resolution_ratio2) 1.0
    }

    #
    # now set the default fonts
    #
    foreach widclass {
	{button Button} {checkbutton Checkbutton} {entry Entry} 
	{label Label} {listbox Listbox} {menu Menu} 
	{menubutton Menubutton}	{message Message} 
	{radiobutton Radiobutton} {scale Scale} {text Text}} {
	set widget [lindex $widclass 0]
	set class  [lindex $widclass 1]
	$widget .test_widget
	set font \
	    [ModifyFont [.test_widget cget -font] .test_widget \
		 -size $xcFonts(normal_size) -default 1]

	###############################################
	if { $xcMisc(resolutionX) < 1000 } {
	    option add *$class.font $font interactive
	}
	###############################################
	
	if { $widget == "button" } {
	    set xcFonts(normal) $font
	    set xcFonts(small) [ModifyFont [.test_widget cget -font] .test_widget \
				    -size $xcFonts(small_size) -default 1]
	    set xcFonts(big) [ModifyFont [.test_widget cget -font] .test_widget \
				  -size $xcFonts(big_size) -default 1]
	} elseif { $widget == "entry" } {
	    set xcFonts(normal_entry) $font
	    set xcFonts(small_entry) \
		    [ModifyFont [.test_widget cget -font] .test_widget \
		    -size $xcFonts(small_size) -default 1]
	}
	destroy .test_widget
    }
    
    # configure the balloon help
    if { [info commands "DynamicHelp::configure"] != "" } {
	DynamicHelp::configure -delay 100 -font $xcFonts(normal)
    }

    xcDebug "Taking settings for $xcMisc(resolution) screen-resolution !!!"
    xcDebug "-------------------"
    xcDebug "Small  Font Size set to:   $xcFonts(small_size)"
    xcDebug "Normal Font Size set to:   $xcFonts(normal_size)"
    #xcDebug "Big    Font Size set to:"
}

proc SetImageSizes {} {
    global xcMisc

    foreach image $xcMisc(rescale_image_list) {
	set h [image height $image]
	set w [image width  $image]
	#xcDebug "h=$h"

	set f $xcMisc(resolution_ratio2)
	if { $f < 0.75 } { set f 0.75 }

	set hs [expr round( $h * $f )]
	set ws [expr round( $w * $f )]
	
	set x1 [expr {($w - $ws)/2}]
	set y1 [expr {($h - $hs)/2}]
	set x2 [expr {($w - $x1)}]
	set y2 [expr {($h - $y1)}]
	
	image      create photo swap_image -width $ws -height $hs
	swap_image copy   $image -from $x1 $y1 $x2 $y2
	$image     blank
	$image     configure -height $hs -width $ws
	$image     copy   swap_image	
	image      delete swap_image
	
	#$image configure -height $hs -width $ws

	# debugging only
	#xcDebug "Image dim:: ${w}x$h"
	#button .$image -image $image -anchor center
	#pack .$image -side left
    }
}
