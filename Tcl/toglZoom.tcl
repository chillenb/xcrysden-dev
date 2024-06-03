#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/toglZoom.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc toglZoom {title togl} {
    global toglZoom
    
    if { ! [info exists toglZoom(counter)] } {
	set toglZoom(counter) 1
    } else {
	incr toglZoom(counter)
    }

    if { ! [info exists toglZoom($togl)] } {
	set toglZoom($togl) 1.0
    }
    if { ! [info exists toglZoom($togl,zoomStep)] } {
	set toglZoom($togl,zoomStep) 0.05
    }
    if { ! [info exists toglZoom($togl,window)] } {
	set toglZoom($togl,window) .toglzoom$toglZoom(counter)
    }
    
    if { [winfo exists $toglZoom($togl,window)] } {
	return
    }

    set t [xcToplevel $toglZoom($togl,window) Zoom Zoom . 0 0 1]
    
    set f0 [frame $t.0 -class StressText]
    set l0 [label $f0.l -text $title]
    set f1 [frame $t.1]
    set f2 [frame $t.2]
    set f3 [frame $t.3]

    ###############
    # Zoom + Zoom -
    set zoom1 [button $f1.zoom1 -text "Zoom +"]
    set zoom2 [button $f1.zoom2 -text "Zoom -"]

    bind $zoom1 <ButtonPress-1>   [list toglZoom:zoom $togl $zoom1 +] 
    bind $zoom2 <ButtonPress-1>   [list toglZoom:zoom $togl $zoom2 -]
    bind $zoom1 <ButtonRelease-1> [list toglZoom:relB1] 
    bind $zoom2 <ButtonRelease-1> [list toglZoom:relB1]
    
    ######################
    # SCALE for translationStep
    set sc [scale $f2.scale -from 0 -to 0.95 \
		-length 200 \
		-variable toglZoom($togl,zoomStep) \
		-orient horizontal -label "Zoom Step:" -tickinterval 0.25 \
		-digits 3 -resolution 0.05 -showvalue true \
		-width 10]

    set close [button $f3.close -text Close \
		   -command [list CancelProc $t]]

    pack $f0 $l0 $f1 $f2 $f3 -side top -fill both -expand 1 -padx 10 -pady 10
    pack $zoom1 $zoom2 -side left -expand 1 -fill both
    pack $sc -expand 1 -pady 10
    pack $close -expand 1
}

proc toglZoom:zoom {togl button sign} {
    global toglZoom

    if { $sign == "+" } {
	set zoom $toglZoom($togl,zoomStep)
    } else {
	set zoom [expr -1.0 * $toglZoom($togl,zoomStep)]
    }

    #--Sunken
    $button configure -relief sunken
    #--
    set toglZoom(B1down) 1    
    while { $toglZoom(B1down) } {
	# a way-around of BUG
	#if { [string match $toglZoom($togl) "nan"] } {
	#    set toglZoom($togl) 1.0
	#}
	cry_toglzoom $togl $zoom
	update 	
    }
    #--Raised
    $button configure -relief raised
    
    return -code break
}

proc toglZoom:relB1 {} {
    global toglZoom
    set toglZoom(B1down) 0
}

proc toglZoom:discreteZoom {togl sign} {
    global toglZoom

    if { $sign == "+" } {
	set zoom $toglZoom($togl,zoomStep)
    } else {
	set zoom [expr -1.0 * $toglZoom($togl,zoomStep)]
    }
    
    cry_toglzoom $togl $zoom
    update 	
}
