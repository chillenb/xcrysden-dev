#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 177 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 177 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/TclisoControl.tcl                   
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc fontToplevelWidget {w labeltext {font {}} {allFonts {}}} {

    xcToplevel $w "Set Font Attributes" Font
    set wid [fontWidget $w.wid $labeltext $font $allFonts]
    pack $wid
    
    global $wid
    upvar \#0 $wid clientData
    set clientData(update) 0
    
    set upd [button $w.update -text "Close & Update" \
		 -command [list _fontWidgetDestroy $wid update]]
    set can [button $w.can    -text "Cancel" \
		 -command [list _fontWidgetDestroy $wid]]
    
    pack $can $upd -side left -expand 1 -padx 5 -pady 5
    
    tkwait window $wid

    set update $clientData(update)
    set font $clientData(font)
    #puts stderr "*** fontToplevelWidget :   font = $font"
    #puts stderr "*** fontToplevelWidget :   update = $update"
    unset clientData
    destroy $w
    
    if { $update } {
	return $font	
    } else {
	return {}
    }
}

proc _fontWidgetDestroy {w {update {}}} {
    global $w
    upvar \#0 $w clientData
    
    if { $update != {} } {
	set clientData(update) 1
    }
    # destrot the fontWidget
    destroy ${w}
}
proc _fontWidgetUnset {w} {
    # destroy the fontWidget data
    unset $w
}
    
proc fontWidget {w labeltext {font {}} {allFonts {}}} {
    global system
    global $w
    upvar \#0 $w data
    
    if { $font == {} } {
	#
	# create new font
	#
	set data(font) [font create]
	font configure $data(font) -size 14 -family $data(family)
    } elseif { $font eq "fixed" } {
	set data(font) [font create]
	eval {font configure $data(font)} [font actual $font]
    } else {
	set data(font) $font
    }

    set data(size)   [font actual $data(font) -displayof . -size]    
    #set data(family) [font actual $data(font) -displayof . -family]
    set data(family) [X11Font_MatchFamily [font actual $data(font) -displayof . -family]]

    if { $data(size) < 0 } {
	set scale [tk scaling]
	set data(size) [expr $data(size) * $scale]
    }

    if { [font actual $data(font) -displayof . -weight] == "normal" } {
	set data(bold) 0
    } else {
	set data(bold) 1
    }
    if { [font actual $data(font) -displayof . -slant] == "roman" } {
	set data(italic) 0
    } else {
	set data(italic) 1
    }	    
    set data(underline)  \
	[font actual $data(font) -displayof . -underline]  
    set data(overstrike) \
	[font actual $data(font) -displayof . -overstrike]

    set data(w) ${w}
    set cont [frame $data(w)]

    set f1   [frame $cont.f1 -relief groove -bd 2]
    pack $f1 -side top -expand 1 -fill both -padx 5 -pady 5
    
    # label for diplaying sample text
    set f1lf        [frame $f1.f1 \
			 -relief flat \
			 -class StressText \
			 -width  300 \
			 -height 100 \
			 -relief groove \
			 -bg \#ffffff \
			 -bd 2]
    set data(label) [label $f1lf.l \
			 -text $labeltext \
			 -font $data(font) \
			 -relief flat \
			 -bg \#ffffff \
			 -bd 2]
    pack $f1lf -side top -expand 0 -padx 5 -pady 5 -ipadx 0 -ipady 0
    pack propagate $f1lf false
    pack $data(label) -side top -expand 1 -ipadx 10 -ipady 2 -padx 10 -pady 5
    
    # container frames
    set f11  [frame $f1.1]
    set f12  [frame $f1.2]
    set f121 [frame $f1.2.1]
    set f122 [frame $f1.2.2]
    pack $f11  $f12  -side top  -fill both -padx 2 -pady 2
    pack $f121 $f122 -side left -fill both -padx 2

    #
    # Font Size (frame #1)
    #
    scale $f11.scale -from 4 -to 100 \
	-length 300 \
	-variable ${w}(size) \
	-orient horizontal -label "Font size:" \
	-digits 3 \
	-resolution 1 \
	-showvalue true \
	-highlightthickness 0 \
	-command [list fontUpdate $w]
    pack $f11.scale -fill y -expand 1

    #
    # Font Family (frame #2)
    #
    if { $allFonts != "" } {
	if { ! [info exists data(ComboBox.currentValue)] } {
	    set data(ComboBox.currentValue) $data(family)
	}
	set label       [label $f121.l -text "Font Family:"]
	set data(combo) [ComboBox $f121.combo \
			     -height   15 \
			     -values   [font families] \
			     -text     $data(ComboBox.currentValue) \
			     -editable false \
			     -modifycmd [list fontUpdate $w -1]]
	pack $label $data(combo) -side left -fill x -expand 1 -padx 1 -pady 3
    } else {
	#set fixed [font actual fixed -family]
	set familyMenu [xcMenuButton $f121 -labeltext "Font Family:" -labelwidth 12 \
			    -textvariable ${w}(family) \
			    -menu [list \
				       "Helvetica"  [list fontUpdate $w -1 Helvetica] \
				       "Times"      [list fontUpdate $w -1 Times] \
				       "Courier"    [list fontUpdate $w -1 Courier] \
				       "Fixed"      [list fontUpdate $w -1 Fixed]] 
			    ]
	pack $familyMenu -side left -fill x -expand 1
    }
    
    #
    # Miscellaneous attributes (frame #3)
    #
    xcCheckButtonRow $f122 4 \
	[list \
	     @$system(BMPDIR)/bold.xbm \
	     @$system(BMPDIR)/italic.xbm \
	     @$system(BMPDIR)/underline.xbm \
	     @$system(BMPDIR)/overstrike.xbm] \
	[list \
	     ${w}(bold) \
	     ${w}(italic) \
	     ${w}(underline) \
	     ${w}(overstrike)] \
	[list \
	     [list fontUpdate $w] \
	     [list fontUpdate $w] \
	     [list fontUpdate $w] \
	     [list fontUpdate $w]] \
	left
    return $w
}


proc fontUpdate {w {size -1} {family {}}} {
    global $w
    upvar \#0 $w data
 
    if { $data(bold) == 0 } {
	set weight normal
    } else {
	set weight bold
    }

    if { $data(italic) == 0 } {
	set slant roman
    } else {
	set slant italic
    }
    
    if { $size > -1 } {
	set data(size) $size
    }
    #if { [info exists data(combo)] } {	
    #	set data(family) [$data(combo) cget -text]
    #}
    if { $family != {} } {
	set data(family) $family
    }
        
    font configure  $data(font) \
	-size       $data(size) \
	-family     $data(family) \
	-weight     $weight \
	-slant      $slant \
	-underline  $data(underline) \
	-overstrike $data(overstrike)

    $data(label) config -font $data(font)
}
