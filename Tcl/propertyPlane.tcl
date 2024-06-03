#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/propertyPlane.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc IsoControl2D {} {
    global isosurf nxdir nydir nzdir periodic xcFonts XCTrace prop \
	    isoControl unmapWin
    
    if { [winfo exists .iso2D] } { return }

    #
    # initializations
    #
    set isoControl(plane) {}
    set isoControl(datagridDim) 2
    set isoControl(colorplane)  1
    set isoControl(isoline)     1

    IsoControl_InitVar
    reloadButtonDisable
    
    set t [xcToplevel .iso2D \
	    "Property-plane Controls" "Property-PlaneControls" . -0 0 1]
    
    xcRegisterUnmapWindow $t $unmapWin(frame,main) \
	    isosurf_control -image unmap-isosurf
    bind $t <Unmap> [list xcUnmapWindow unmap %W $t \
	    $unmapWin(frame,main) isosurf_control]
    bind $t <Map>   [list xcUnmapWindow map %W $t \
	    $unmapWin(frame,main) isosurf_control]
    
    ########################
    ### ColorPLANE frame ###
    ########################
    set f0 [frame $t.f0 -relief raised -bd 2]
    frame $f0.1 -relief groove -bd 2
    pack $f0 -side top -fill both -expand 1
    pack $f0.1 -padx 2 -pady 5 -ipady 3 -fill x
 
    scale $f0.1.sc -from 1 -to 5 -length 100 \
	    -variable isosurf(3Dinterpl_degree) -orient horizontal \
	    -label "Degree of biCubic Spline:" \
	    -tickinterval 1 -resolution 1 \
	    -width 7 -sliderlength 20 \
	    -showvalue true \
	    -font $xcFonts(small)
    pack $f0.1.sc -side top -padx 2 -pady 5 -fill both -expand 1
    

    set f1 [frame $t.f1 -relief raised -bd 2]
    frame $f1.1 -relief groove -bd 2
    set f11 [frame $f1.1.f1]
    set f12 [frame $f1.1.f2]
    xcMenuEntry $f11 "Select color basis:" 30  \
	    isoControl(cpl_basis) {MONOCHROME RAINBOW RGB GEOGRAPHIC BLUE-WHITE-RED BLACK-BROWN-WHITE} \
	    -labelwidth 17 \
	    -labelfont $xcFonts(small) \
	    -entryfont $xcFonts(small_entry) \
	    -labelanchor w
    xcMenuEntry $f12 "Select scale function:" 30  \
	    isoControl(cpl_function) \
	    {LINEAR LOG10 SQRT 3th-ROOT EXP(x) EXP(x^2)} \
	    -labelwidth 17 \
	    -labelfont $xcFonts(small) \
	    -entryfont $xcFonts(small_entry) \
	    -labelanchor w
    pack $f1 -side top -fill both -expand 1
    pack $f1.1 -padx 2 -pady 5 -ipady 3 -fill x
    pack $f11 $f12 -side top -in $f1.1 -fill x
    
    #
    # DISPLAY/RANGES/EXPAND/ISOLINE
    #    
    set f2  [frame $t.f2 -relief raised -bd 2]
    pack $f2 -side top -fill both -expand 1 
    
    IsoControl_DispRanExpIso $f2
    
    ########################################
    # BOTTOM FRAME
    set f3  [frame $t.f3 -relief raised -bd 2]
    pack $f3 -side top -fill x -expand 1
    
    set hid [button $f3.hid -text "Hide" \
	    -command [list IsoControl_Hide $t]]
    set can [button $f3.can -text "Close" \
	    -command [list IsoControlCan $t 2]]
    set sav [button $f3.sav -text "Save Grid" \
	    -command IsoControlSave]
    set sub [button $f3.sub -text "Submit" \
	    -command UpdatePropertyPlane]
    pack $hid $can $sav $sub -side left -pady 5 -expand 1
}


proc IsoControl_DispRanExpIso f {
    global isosurf nxdir nydir nzdir periodic xcFonts XCTrace prop \
	    isoControl xcColors
    
    set f11 [frame $f.1]
    pack $f11 -side top -fill both -expand 1 -padx 4 -pady 8
    set mfb [frame $f11.mfb]
    set mf  [frame $f11.mf -relief raised -bd 2]
    pack $mfb $mf -side top -padx 2 -fill both -expand 1
    
    set df  [frame $mf.d -relief groove -bd 2]
    set rf  [frame $mf.r -relief groove -bd 2]
    set ef  [frame $mf.e -relief groove -bd 2]
    set ifr [frame $mf.i]
    
    set d [button $mfb.d -text "Display" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow display \
	    $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr]]
    set r [button $mfb.r -text "Ranges" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow ranges \
	    $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr]]
    set e [button $mfb.e -text "Expand" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow expand \
	    $mfb.d $mfb.r $mfb.e $mfb.i\
	    $df $rf $ef $ifr]]
    set i [button $mfb.i -text "Isoline" \
	    -highlightthickness 0 -bd 1 \
	    -command [list IsoControl_IsoLineShow isoline \
	    $mfb.d $mfb.r $mfb.e $mfb.i\
	    $df $rf $ef $ifr]]

    if { $periodic(dim) == 0 } {
	$e config -state disabled
    }
    pack $d $r $e $i -side left -padx 0 -pady 0 -fill both -expand 1
    
    IsoControl_IsoLineShow display $mfb.d $mfb.r $mfb.e $mfb.i \
	    $df $rf $ef $ifr
    
    
    #
    # DISPLAY
    #
    set l     [label $df.l -text "Property-plane display option:"]
    set lfont [ModifyFont [$df.l cget -font] $df.l -underline 1]
    $l configure -font $lfont
    
    set ck1 [checkbutton $df.c1 \
	    -text "display color-plane" \
	    -variable isoControl(colorplane) \
	    -width 21 \
	    -anchor w]
    set ck2 [checkbutton $df.c2 \
	    -text "display isolines" \
	    -variable isoControl(isoline) \
	    -width 21 \
	    -anchor w]
    set ck3 [checkbutton $df.c3 \
		 -text "lighting of color-plane" \
		 -variable isoControl(colorplane_lighting) \
		 -width 21 \
		 -onvalue 1 \
		 -offvalue 0 \
		 -anchor w]
    set ck4 [checkbutton $df.c4 \
	    -text "transparent color-plane" \
	    -variable isoControl(cpl_transparency) \
	    -width 21 \
	    -onvalue 1 -offvalue 0 \
	    -anchor w]
    set ck5 [checkbutton $df.c5 \
		 -text "display thermometer" \
		 -variable isoControl(cpl_thermometer) \
		 -width 21 \
		 -onvalue 1 -offvalue 0 \
		 -anchor w]
    set ck6 [checkbutton $df.c6 \
		 -text "thermometer in toplevel" \
		 -variable isoControl(cpl_thermoTplw) \
		 -width 21 \
		 -onvalue 1 -offvalue 0 \
		 -anchor w]
    set thermo [frame $df.f]
    grid configure $l   -column 0 -row 0 -columnspan 2 
    grid configure $ck1 -column 0 -row 1
    grid configure $ck4 -column 0 -row 2
    grid configure $ck2 -column 1 -row 1
    grid configure $ck3 -column 1 -row 2
    grid configure $ck5    -column 0 -row 3
    grid configure $ck6    -column 1 -row 3
    grid configure $thermo -column 0 -row 4 -columnspan 2

    # thermometer-widgets
    set tf1  [frame $thermo.1 -relief groove -bd 2]
    set tf11 [frame $tf1.1]
    set tf12 [frame $tf1.2]
    label $tf11.__l -text "Thermometer settings:" -anchor w
    pack $tf11.__l -side top -fill x
    FillEntries $tf11 {
	"Format string:" 
	"Label:" 
	"No. of tics:"
    } {
	isoControl(cpl_thermoFmt)
	isoControl(cpl_thermoLabel)
	isoControl(cpl_thermoNTics)
    } 14 20
    button $tf12.font -text "Set Font" -command isoControl_thermoFont
    pack $tf1 -side top -fill both -expand 1 -padx 5 -pady 5
    pack $tf11 $tf12 -side left -fill x -expand 1 -padx 5 -pady 5
    pack $tf12.font -side top -fill x -expand 1

    #
    # RANGES
    #
    set frame_color [lindex \
	    [GetWidgetConfig frame -background] end]
    FillEntries $rf {"Minimum 3D grid value:"} \
	    isosurf(minvalue) \
	    25 9 top left \
	    -e_relief flat -e_state disabled -e_bg $frame_color \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)
    FillEntries $rf {"Maximum 3D grid value:"} \
	    isosurf(maxvalue) \
	    25 9 top left \
	    -e_relief flat -e_state disabled -e_bg $frame_color \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)

    set isoControl(2Dlowvalue_entry) [FillEntries $rf \
	    {"Lowest rendered value:"} isoControl(2Dlowvalue) \
	    25 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    set isoControl(2Dhighvalue_entry) [FillEntries $rf \
	    {"Highest rendered value:"} isoControl(2Dhighvalue) \
	    25 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    set isoControl(2Dnisoline_entry) [FillEntries $rf \
	    {"Number of isolines:"} \
	    isoControl(2Dnisoline) \
	    25 9 top left \
	    -l_font $xcFonts(small) -e_font $xcFonts(small_entry)]
    
    
    #
    # EXPAND
    #
    set f21  [frame $ef.top]
    set f22  [frame $ef.bot]
    pack $f21 $f22 -side top -fill both -expand 1

    set enable "\
	    $f22.scX configure -foreground $xcColors(enabled_fg); \
            $f22.scY configure -foreground $xcColors(enabled_fg); \
            $f22.scZ configure -foreground $xcColors(enabled_fg)"
    set disable " \
	    $f22.scX configure -foreground $xcColors(disabled_fg); \
            $f22.scY configure -foreground $xcColors(disabled_fg); \
            $f22.scZ configure -foreground $xcColors(disabled_fg)"
    if { $periodic(dim) > 0 } {
	label $f21.l -text "Expand Property-plane:" -relief flat
	$f21.l configure -font $lfont
	set r1 [radiobutton $f21.r1 \
		-text "do not expand" \
		-variable isosurf(2Dexpand) \
		-value "none" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f22; catch {eval $disable}"]
	set r2 [radiobutton $f21.r2 \
		-text "to whole structure" \
		-variable isosurf(2Dexpand) \
		-value "whole" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcDisableAll $f22; catch {eval $disable}"]
	set r3 [radiobutton $f21.r3 \
		-text "separately in each direction" \
		-variable isosurf(2Dexpand) \
		-value "specify" \
		-anchor w \
		-font $xcFonts(small) \
		-command "xcEnableAll $f22; catch {eval $enable}"]
	
	set XCTrace(2DscX) [scale $f22.scX -from 1 -to $nxdir -length 100 \
		-variable isosurf(2Dexpand_X) -orient horizontal \
		-label "repeat in X-dir:" -tickinterval 1 -resolution 1 \
		-showvalue true \
		-font $xcFonts(small) \
		-width 7 -sliderlength 20]
	# TRACE -- nxdir
	trace variable nxdir w xcTrace
	pack $XCTrace(2DscX) -side left -pady 5 -expand 1 -fill x
	if { $periodic(dim) > 1} {
	    set XCTrace(2DscY) [scale $f22.scY -from 1 -to $nydir \
		    -length 100 \
		    -variable isosurf(2Dexpand_Y) -orient horizontal \
		    -label "repeat in Y-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nydir
	    trace variable nydir w xcTrace
	    pack $XCTrace(2DscY) -side left -pady 5 -expand 1 -fill x
	}
	if { $periodic(dim) > 2 } {
	    set XCTrace(2DscZ) [scale $f22.scZ -from 1 -to $nzdir \
		    -length 100 \
		    -variable isosurf(2Dexpand_Z) -orient horizontal \
		    -label "repeat in Z-dir:" -tickinterval 1 -resolution 1 \
		    -showvalue true \
		    -font $xcFonts(small) \
		    -width 7 -sliderlength 20]
	    # TRACE -- nzdir
	    trace variable nzdir w xcTrace
	    pack $XCTrace(2DscZ) -side left -pady 5 -expand 1 -fill x 
	}
	# get correct state for EXPAND options
	if { $isosurf(2Dexpand) == "none" || $isosurf(2Dexpand) == "whole" } {
	    xcDisableAll $f22; catch {eval $disable}
	} else {
	    xcEnableAll $f22; catch {eval $enable}
	}
	pack $f21.l -side top -expand 1
	pack $r1 $r2 $r3 -side top -fill x

    }

    #
    # Isoline
    #
    set mf1 [frame $ifr.1 -relief groove -bd 2]
    set mf2 [frame $ifr.2 -relief groove -bd 2]
    set mf3 [frame $ifr.3 -relief groove -bd 2]
    pack $mf1 $mf2 $mf3 -side top -expand 1 -fill both -padx 3 -pady 3
    
    set isoControl(bmc) [button $mf1.b -text "set color" \
	    -command IsoControl_SetIsolineColor]
    
    RadioButVarCmd $mf1 "Isoline Color:" isoControl(isoline_color) \
	    IsoControl_IsolineColor left top 0 0 \
	    {monocolor} {property color}
    pack $isoControl(bmc) -side left -padx 2

    RadioBut $mf2 "Isoline Stipple:" isoControl(isoline_stipple) left top 0 1 \
	    {no stipple} {stipple negative} {full stipple}    
    # this is temporal, since "Isoline Stipple" is not yet working
    xcDisableAll $mf2

    # here goes isoline_width entry
    set isoControl(2Disolinewidth_entry) [FillEntries $mf3 \
					      {"Isoline width:"} \
					      isoControl(isoline_width) \
					      14 10 left left]
}



proc UpdatePropertyPlane {{var {}}} {
    global prop isosurf periodic isosign isodata isosurf_struct \
	    isoControl isoControl_struct vec dif_isosurf

    SetWatchCursor

    if { ![check_var {
	{isoControl(2Dlowvalue)    real} 
	{isoControl(2Dhighvalue)   real}
	{isoControl(2Dnisoline)    posint}
	{isoControl(isoline_width) posint}
    } [list \
	   $isoControl(2Dlowvalue_entry) \
	   $isoControl(2Dhighvalue_entry) \
	   $isoControl(2Dnisoline_entry) \
	   $isoControl(2Disolinewidth_entry)]] } {
	ResetCursor
	return
    }
    if { $isoControl(2Dnisoline) > $isoControl(max_allowed_2Dnisoline) } {
	tk_dialog [WidgetName] WARNING "WARNING: more then $isoControl(max_allowed_2Dnisoline) isolines was requested; Maximum number is $isoControl(max_allowed_2Dnisoline) !!! Setting number of isolines to $isoControl(max_allowed_2Dnisoline)" \
		warning 0 OK
	set isoControl(2Dnisoline) $isoControl(max_allowed_2Dnisoline)
    }

    if { $isoControl(cpl_basis) == {} } {
	tk_dialog [WidgetName] ERROR "ERROR: You forgot to specify color basis. Please Do it !" error 0 OK
	ResetCursor
	return
    }
    if { $isoControl(cpl_function) == {} } {
	tk_dialog [WidgetName] ERROR "ERROR: You forgot to specify scale function. Please Do it !" error 0 OK
	ResetCursor	
	return
    }

    ###################################
    # this is needed if spin is changed
    set upd 0; 
    if { $isosurf_struct(isosign)    != $isosign || \
	    $isosurf_struct(isodata) != $isodata || \
	    $isosurf_struct(spin)    != $isosurf(spin) } {	    	    
	SetXC_Iso 2D
	eval $isosign
	xc_iso end isosign
	eval $isodata
	xc_isodata make
	set upd 1
    }

    if { $isosurf_struct(3Dinterpl_degree) != $isosurf(3Dinterpl_degree) } {
	xc_iso interpolate $isosurf(3Dinterpl_degree)
	set upd 1
    }
    
    if { \
	     $isoControl_struct(cpl_basis)        != $isoControl(cpl_basis) || \
	     $isoControl_struct(cpl_function)     != $isoControl(cpl_function) || \
	     $isoControl_struct(colorplane)       != $isoControl(colorplane) || \
	     $isoControl_struct(isoline)          != $isoControl(isoline)	     || \
	     $isoControl_struct(colorplane_lighting) != $isoControl(colorplane_lighting)	     || \
	     $isoControl_struct(cpl_transparency) != $isoControl(cpl_transparency)  || \
	     $isoControl_struct(cpl_thermometer) != $isoControl(cpl_thermometer) || \
	     $isoControl_struct(cpl_thermoTplw)  != $isoControl(cpl_thermoTplw) || \
	     $isoControl_struct(cpl_thermoFmt)   != $isoControl(cpl_thermoFmt) || \
	     $isoControl_struct(cpl_thermoLabel) != $isoControl(cpl_thermoLabel) || \
	     $isoControl_struct(cpl_thermoNTics) != $isoControl(cpl_thermoNTics) || \
	     $isoControl_struct(cpl_thermoFont)  != $isoControl(cpl_thermoFont) || \
	     $isosurf_struct(2Dexpand)            != $isosurf(2Dexpand)  	     || \
	     $isosurf_struct(2Dexpand_X)          != $isosurf(2Dexpand_X)	     || \
	     $isosurf_struct(2Dexpand_Y)          != $isosurf(2Dexpand_Y)	     || \
	     $isosurf_struct(2Dexpand_Z)          != $isosurf(2Dexpand_Z)	     || \
	     $isoControl_struct(anim_step)        != $isoControl(anim_step)	     || \
	     $isoControl_struct(current_slide)    != $isoControl(current_slide)     || \
	     $isoControl_struct(2Dlowvalue)       != $isoControl(2Dlowvalue)        || \
	     $isoControl_struct(2Dhighvalue)      != $isoControl(2Dhighvalue)	     || \
	     $isoControl_struct(2Dnisoline)       != $isoControl(2Dnisoline)	     || \
	     $isoControl_struct(isoline_color)    != $isoControl(isoline_color)     || \
	     $isoControl_struct(isoline_width)    != $isoControl(isoline_width)     || \
	     $isoControl_struct(isoline_monocolor)!= $isoControl(isoline_monocolor) || \
	     $isoControl_struct(isoline_stipple)  != $isoControl(isoline_stipple) } {
	set upd 1
    }

    if { $upd } {
	if { $isoControl(colorplane) && !$isoControl(isoline) } {
	    set what "colorplane"
	} elseif { !$isoControl(colorplane) && $isoControl(isoline) } {
	    set what "isoline"
	} elseif { $isoControl(colorplane) && $isoControl(isoline) } {
	    set what "both"
	} else {
	    set what "none"
	}
	
	#
	# this must be the same as in "isosurf.h"
	#
	set COLORBASE_MONO           0
	set COLORBASE_RAINBOW        1
	set COLORBASE_RGB            2
	set COLORBASE_GEOGRAPHIC     3
	set COLORBASE_BLUE_WHITE_RED 4
	set COLORBASE_BLACK_BROWN_WHITE 5
	
	set SCALE_FUNC_LIN      0
	set SCALE_FUNC_LOG      1
	set SCALE_FUNC_LOG10    2
	set SCALE_FUNC_SQRT     3
	set SCALE_FUNC_ROOT3    4
	set SCALE_FUNC_GAUSS    5
	set SCALE_FUNC_SLATER   6
	
	set basis COLORBASE_MONO
	set func  SCALE_FUNC_LIN

	switch -exact -- $isoControl(cpl_basis) {
	    MONOCHROME     { set basis $COLORBASE_MONO }
	    RAINBOW        { set basis $COLORBASE_RAINBOW }
	    RGB	           { set basis $COLORBASE_RGB }
	    GEOGRAPHIC     { set basis $COLORBASE_GEOGRAPHIC }
	    BLUE-WHITE-RED { set basis $COLORBASE_BLUE_WHITE_RED }
	    BLACK-BROWN-WHITE { set basis $COLORBASE_BLACK_BROWN_WHITE }
	}
	switch -exact -- $isoControl(cpl_function) {
	    LINEAR     { set func $SCALE_FUNC_LIN    } 
	    LOG        { set func $SCALE_FUNC_LOG    } 
	    LOG10      { set func $SCALE_FUNC_LOG10  } 
	    SQRT       { set func $SCALE_FUNC_SQRT   } 
	    3th-ROOT   { set func $SCALE_FUNC_ROOT3  } 
	    EXP(x)     { set func $SCALE_FUNC_SLATER } 
	    EXP(x^2)   { set func $SCALE_FUNC_GAUSS  } 
	}
	
	if { $isoControl(isoline_color) == "monocolor" } {
	    set linecolor [concat monocolor \
		    [rgb_h2f $isoControl(isoline_monocolor)]]
	} else {
	    set linecolor polycolor
	}
	switch -exact -- $isoControl(isoline_stipple)  {
	    {no stipple}        { set linedash nodash   }	   
	    {stipple negative}  { set linedash negdash  }
	    {full stipple}      { set linedash fulldash }
	}
	
	xc_iso isoplaneconfig 0 \
	    -isoplanemin      $isoControl(2Dlowvalue) \
	    -isoplanemax      $isoControl(2Dhighvalue) \
	    -isolinecolor     $linecolor \
	    -isolinedash      $linedash \
	    -isolinewidth     $isoControl(isoline_width) \
	    -isolinenlevels   $isoControl(2Dnisoline) \
	    -isoplanelighting $isoControl(colorplane_lighting)
	
	xc_iso isoplane 0 $basis $func $what 
	if { $what != "none" } {
	    xc_isoplane .mesa 0 -planetype $what \
		    -transparency $isoControl(cpl_transparency) \
		    -render after
	    
	    ########################################
	    # now take care of xc_isoexpand
	    # isosurf(expand) ......... none/whole/specify
	    #                                      |-> isosurf(expand_[X,Y,Z])
	    if { $isosurf(2Dexpand) == "none" } {
		if { $periodic(dim) == 3 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} elseif { $periodic(dim) == 2 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} elseif { $periodic(dim) == 1 } {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		} else {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		}
	    } elseif { $isosurf(2Dexpand) == "whole" } {
		set expand2D "whole"
	    } else {
		# expand == specify
		if { $periodic(dim) == 3 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf(2Dexpand_X) \
			    $isosurf(2Dexpand_Y) \
			    $isosurf(2Dexpand_Z)]
		} elseif { $periodic(dim) == 2 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf(2Dexpand_X) $isosurf(2Dexpand_Y) 0]
		} elseif { $periodic(dim) == 1 } {
		    set expand2D [list 1 0 0  0 1 0  0 0 1 \
			    $isosurf(2Dexpand_X) 0 0]
		} else {
		    set expand2D {1 0 0  0 1 0  0 0 1    0 0 0}
		}
	    }
	
	    xcDebug "Expand2D::::periodic == $periodic(dim)"
	    xcDebug "Expand2D::::expand2D == $expand2D"
	    xc_isoexpand .mesa 0 -repeattype default \
		    -shape parapipedal -expand $expand2D -render after
	}

	############################
	# take care of thermometer #
	############################
	if { $isoControl(cpl_thermometer) } {
	    # render thermometer !!!	   
	    foreach wid [list .mesa.thermo .thermo] {
		if { [winfo exists $wid] } {
		    destroy $wid
		}
	    }
	    if { $isoControl(cpl_thermoTplw) == 1 } {
		set prefix {}
	    } else {
		set prefix .mesa
	    }
	    set tw [thermometerWidget $prefix.thermo {} \
			$isoControl(cpl_thermoLabel) \
			$isoControl(cpl_thermoFont) \
			$isoControl(cpl_thermoFmt) \
			$isoControl(cpl_basis) \
			$isoControl(cpl_function) \
			$isoControl(2Dlowvalue) \
			$isoControl(2Dhighvalue) \
			$isoControl(cpl_thermoNTics) \
			$isoControl(cpl_thermoTplw)]
	    if { $isoControl(cpl_thermoTplw) == 0 } {
		pack $tw -anchor nw -side top
	    }
	} else {
	    if { [winfo exists .mesa.thermo] } {
		destroy .mesa.thermo
	    } elseif { [winfo exists .thermo] } {
		destroy .thermo
	    }
	}
    }
    
    # update the display
    xc_display .mesa

    # now updata isosurf_struct global variable
    Set_UpdateIsosurf_Struct
    
    ResetCursor
}
