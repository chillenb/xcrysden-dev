#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/forces.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc forceVectors can {
    global check sInfo
    
    if { ! [info exists sInfo(lforce)] } {
	return
    }
    if { ! $sInfo(lforce) } {
	return
    }

    if { $check(forces) } {
	xc_forces $can on
	$can render
	xcSwapBuffers
    } else {
	xc_forces $can off
	$can render
	xcSwapBuffers
    }
}

proc forceVectorsSet can {
    global check forceVec sInfo

    if { ! [info exists sInfo(lforce)] } {
	return
    }
    
    set c [lindex [split $can .] end]
    if { [winfo exists .force$c] } {
	return
    }

    if { ![info exists forceVec(scalefunction)] } {
	global mody
	set forceVec(scalefunction) linear
	set forceVec(threshold)     0.0005
	set forceVec(lengthfactor)  200

	set forceVec(rod_thickf) [xc_getdefault $mody(D_FORCE_RODTHICKF)]
	set forceVec(arr_thickf) [xc_getdefault $mody(D_FORCE_ARRTHICKF)]
	set forceVec(arr_lenf)   [xc_getdefault $mody(D_FORCE_ARRLENF)]
	set forceVec(color)      [xc_getdefault $mody(D_FORCE_COLOR)]

	puts stderr "forceVec: color == $forceVec(color)"
    }

    # from here on: widget managing
    set t [xcToplevel .force$c "Forces: Settings" "Forces" . -0 0 1]
    set f1 [frame $t.f1 -relief ridge -bd 2]
    set f2 [frame $t.f2 -relief ridge -bd 2]
    set f3 [frame $t.f3]
    pack $f1 $f2 $f3 -side top -expand 1 -fill x -padx 3m -pady 3m -ipadx 3m -ipady 3m
    pack $f3 -side top -expand 1

    #
    # FRAME-1: Scale-Function + Threshold + Length-factor
    #
    set mb [xcMenuButton $f1 \
	    -labeltext  "Scale Function:" \
	    -labelwidth 15 \
	    -textvariable forceVec(scalefunction) \
	    -menu { \
	    {Linear}             {set forceVec(scalefunction) linear} \
	    {Natural Logarithm}  {set forceVec(scalefunction) log} \
	    {Decadic Logaritm}   {set forceVec(scalefunction) log10} \
	    {Square Root}        {set forceVec(scalefunction) sqrt} \
	    {Cubic  Root}        {set forceVec(scalefunction) root3} \
	    {Exponential}        {set forceVec(scalefunction) exp} \
	    {exp(x*x)}           {set forceVec(scalefunction) exp2} \
	}]

    set m1 [message $f1.m1 -justify left -aspect 800 \
	    -text "Threshold::\ndo not show the forces whose magnitude is smaller then the specified threshold.\n\n\nLength Factor::\nassumed unit for force is Hartree/Angstrom\nforceVector_length == forceVectorSize * LengthFactor"]
    pack $mb $m1 -side top -fill x -anchor w -expand 1
    FillEntries $f1 {Treshold {Length Factor}} \
	    {forceVec(threshold) forceVec(lengthfactor)} 13 10 top left

    #
    # FRAME-2: attributes + Threshold + Length-factor
    #
    FillEntries $f2 {
	"Vectors thickness factor:"
	"Thickness factor for arrow-cap:"
	"Length factor for arrow-cap:"
    } {forceVec(rod_thickf) forceVec(arr_thickf) forceVec(arr_lenf)} 31 10 top left
    button $f2.button -text "Set vector's color" -command [list forceVectors_Color $t]
    button $f2.reset  -text "Reset vector's attributes" -command forceVectors_ResetAttributes

    pack $f2.button -side left  -padx 5 -pady 5
    pack $f2.reset  -side right -padx 5 -pady 5

    #
    # FRAME-3: Close + Update buttons
    #

    set b1 [button $f3.close -text "Close" -command [list DestroyWid $t]]
    set b2 [button $f3.update -text "Update" \
	    -command [list forceUpdate $can]]
    
    pack $b1 $b2 -side left -padx 3m -pady 3m -expand 1
}

proc forceUpdate can {
    global forceVec mody
    
    xc_forces $can scalefunction  $forceVec(scalefunction)
    xc_forces $can threshold      $forceVec(threshold)
    xc_forces $can lengthfactor   $forceVec(lengthfactor)

    xc_newvalue .mesa $mody(R_FORCE_RODTHICKF) $forceVec(rod_thickf)
    xc_newvalue .mesa $mody(R_FORCE_ARRTHICKF) $forceVec(arr_thickf)
    xc_newvalue .mesa $mody(R_FORCE_ARRLENF)   $forceVec(arr_lenf)
    eval xc_newvalue .mesa $mody(R_FORCE_COLOR)     $forceVec(color)

    $can render
    xcSwapBuffers
}


proc forceVectors_Color {parent} {
    global forceVec
    
    set t [xcToplevel [WidgetName] \
	       "Set Color of Vectors" "Vector's color" $parent 0 0 1]

    set init_color [rgb_f2h $forceVec(color)]
    xcModifyColor $t "Set Color of Vectors" $init_color \
	groove left left 100 100 70 5 20
    set forceVec(colorID) [xcModifyColorGetID]
    
    set ok  [DefaultButton [WidgetName $t] -text "OK" \
		 -command [list forceVectors_ColorOK $t]]
    set can [button [WidgetName $t] -text "Cancel" \
		 -command [list destroy $t]]

    pack $ok $can -padx 10 -pady 10 -expand 1
}


proc forceVectors_ColorOK {t} {
    global forceVec mody_col mody

    set id $forceVec(colorID)

    set alpha [lindex $forceVec(color) 3]
    set forceVec(color) [list $mody_col($id,red) $mody_col($id,green) $mody_col($id,blue) $alpha]

    # update vector's color
    eval xc_newvalue .mesa $mody(L_FORCE_COLOR) $forceVec(color)
    
    destroy $t
}


proc forceVectors_ResetAttributes {} {
    global forceVec mody

    # reset attributes
    xc_resetvar .mesa $mody(R_FORCE_RODTHICKF)
    xc_resetvar .mesa $mody(R_FORCE_ARRTHICKF)
    xc_resetvar .mesa $mody(R_FORCE_ARRLENF)  
    xc_resetvar .mesa $mody(R_FORCE_COLOR)    

    # now reload default values to forceVec
    set forceVec(rod_thickf) [xc_getdefault $mody(D_FORCE_RODTHICKF)]
    set forceVec(arr_thickf) [xc_getdefault $mody(D_FORCE_ARRTHICKF)]
    set forceVec(arr_lenf)   [xc_getdefault $mody(D_FORCE_ARRLENF)]
    set forceVec(color)      [xc_getdefault $mody(D_FORCE_COLOR)]     
}
