#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/toolbox.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc FS_Toolbox {toolbox togl spin iband {multiband ""}} {
    global FS_Interpolation toglZoom fs
    
    # width == 700 * $xcMisc(resolution_ratio1)

    if { $multiband == "" } {
	#
	# Interpolation
	#
	if { ! [info exists FS_Interpolation($togl)] } {
	    if { [info exists fs($spin,$iband,interpolationdegree)] } {
		set FS_Interpolation($togl) [lindex $fs($spin,$iband,interpolationdegree) 0]
	    } else {
		set FS_Interpolation($togl) 1
	    }
	}

	set f_int [frame $toolbox.int -relief sunken -bd 1] 
	set l1  [label $f_int.l -text "Degree of\nInterpolation:" -font TkFixedFont -anchor w]
	set sc1 [scale $f_int.scale -from 1 -to 6 \
		     -font         TkFixedFont \
		     -length       90 \
		     -variable     FS_Interpolation($togl) \
		     -orient       horizontal \
		     -tickinterval 1 \
		     -digits       1 \
		     -resolution   1 \
		     -showvalue    false \
		     -width        10]
	set submit [button $f_int.submit -text Submit -font TkFixedFont \
			-command [list FS_Interpolation:submit $togl $spin $iband]]
	pack $f_int -side left -padx 1m -ipadx 0m -ipady 1m
	pack $l1 $sc1 $submit -side left -padx 1 -pady 0 -ipadx 0 -ipady 0
    }

    #
    # Zoom
    #
    if { ! [info exists toglZoom($togl)] } {
	set toglZoom($togl) 1.0
    }
    if { ! [info exists toglZoom($togl,zoomStep)] } {
	set toglZoom($togl,zoomStep) 0.05
    }
    
    set f_zoom [frame $toolbox.zoom -relief sunken -bd 1]

    # SCALE for translationStep
    
    set l2  [label $f_zoom.l -text "Zoom\nStep:" -font TkFixedFont -anchor w]
    set sc2 [scale $f_zoom.scale \
		 -from 0 -to 0.4 \
		 -font         TkFixedFont \
		 -length       90 \
		 -variable     toglZoom($togl,zoomStep) \
		 -orient       horizontal \
		 -tickinterval 0.2 \
		 -digits       1 \
		 -resolution   0.05 \
		 -showvalue    false \
		 -width        10]

    # zoom buttons

    set zoom1 [button $f_zoom.zoom1 -text "+" -font TkFixedFont]
    set zoom2 [button $f_zoom.zoom2 -text "-" -font TkFixedFont]

    #bind $zoom1 <ButtonPress-1>   [list toglZoom:zoom $togl $zoom1 +] 
    #bind $zoom2 <ButtonPress-1>   [list toglZoom:zoom $togl $zoom2 -]
    bind $zoom1 <ButtonRelease-1> [list toglZoom:discreteZoom $togl +] 
    bind $zoom2 <ButtonRelease-1> [list toglZoom:discreteZoom $togl -]
    
    pack $f_zoom -side left -padx 1m -ipadx 0m -ipady 1m
    pack $l2 $sc2 $zoom1 $zoom2 -side left -padx 1 -pady 0 -ipadx 0 -ipady 0

    if { $multiband == "" } {
	# Revert Sides + Revert normals buttons
	set f_revert [frame $toolbox.revert -relief sunken -bd 1]
	
	set side [button $f_revert.side -text "Revert Sides"   -font TkFixedFont \
		      -command [list FS_RevertIsoSides $spin $iband]]
	set norm [button $f_revert.norm -text "Revert Normals" -font TkFixedFont \
		      -command [list FS_RevertIsoNormals $spin $iband]]
	pack $f_revert   -side right -padx 1m -ipadx 1m
	pack $side $norm -side top -fill x -expand 1 -padx 1 -pady 0 -ipadx 0 -ipady 0
    }
}


proc FS_RevertIsoSides {spin iband} {
    global fs

    if { $fs($spin,$iband,frontface) == "CCW" } {
	set fs($spin,$iband,frontface) "CW"
    } else { 
	set fs($spin,$iband,frontface) "CCW"
    }
    FS_Config $iband $spin 
}


proc FS_RevertIsoNormals {spin iband} {
    global fs

    if { $fs($spin,$iband,revertnormals)  } {
	set fs($spin,$iband,revertnormals) 0
    } else { 
	set fs($spin,$iband,revertnormals) 1
    }    
    FS_Config $iband $spin 
}
