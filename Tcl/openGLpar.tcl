#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/openGLpar.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

########################################
# openGL(what).......this is "-what" flag for xc_getGLparam/xc_setGLparam 
#                    commands
# openGL(type).......material/light/lightmodel
# openGL(ambient_R, ambient_G, ambient_B, ambient_A
#        diffuse_R, ...
#        specular_R, ...
#        shininess
#        emission_R, ...)
# openGL(src_blend)     SOURCE BLEND FUNCTION
# openGL(dst_blend)     DESTINATION BLEND FUNCTION
# openGL(isoside)       CCW/CW

proc SetOpenGLPar {item type f {test {}}} {
    global openGL xcFonts multi_widget_list xcMisc

    xcDebug "00000000>"
    if { $type == "ISOSURF" } {
	# TRANSPARENCY
	# here QUERY the current BLEND function values 
	switch -glob -- $item {
	    *POS_FRONT_COLOR*  {
		set openGL(type) material
		set openGL(what) {isosurf pos front}
	    }
	    *POS_BACK_COLOR*   {
		set openGL(type) material
		set openGL(what) {isosurf pos back}
	    }
	    *NEG_FRONT_COLOR*  {
		set openGL(type) material
		set openGL(what) {isosurf neg front}
	    }
	    *NEG_BACK_COLOR*   {
		set openGL(type) material
		set openGL(what) {isosurf neg back}
	    }
	    *ONE_FRONT_COLOR*  {
		set openGL(type) material
		set openGL(what) {isosurf one front}
	    }
	    *ONE_BACK_COLOR*   {
		set openGL(type) material
		set openGL(what) {isosurf one back}
	    }
	    *_BLEND_*          {
		set openGL(type) blendfunc
		set openGL(what) isosurf
	    }
	}

	if [string match "*_COLOR_*" $item] {
	    xcDebug "1111111>"
	    set com [concat xc_getGLparam material -what $openGL(what)]
	    set ambient           [eval $com -get ambient]
	    set diffuse           [eval $com -get diffuse]
	    set specular          [eval $com -get specular]
	    set emission          [eval $com -get emission]
	    set openGL(shininess) [eval $com -get shininess]
	    
	    SetRGBA {ambient diffuse specular emission} \
		    [list $ambient $diffuse $specular $emission]
	}
    }

    xcDebug "2222222>"

    if { [string match "*_BLEND_*" $item] } {
	# now trim $item $type for "_" character
	set type [string trim $type _]
	set item [string trim $item _]
	set multi_widget_list(post) $f.f1

	set bldf              [xc_getGLparam blendfunc -what isosurf]
	set openGL(src_blend) [lindex $bldf 0]
	set openGL(dst_blend) [lindex $bldf 1]

	set f1  [frame $f.f1 -relief flat]
	set f11 [frame $f1.1 -relief groove -bd 2]
	set f12 [frame $f1.2 -relief groove -bd 2]
	xcMenuEntry $f11 "Specify source blending type:" 30 \
		openGL(src_blend) $openGL(src_blend_list) \
		-labelwidth 30 -labelanchor w
	button $f11.b1 -text "Default" \
		-command [list GetOpenGLPar _SRC_${item}_${type}_ \
		openGL(src_blend) dummy2 dummy3 dummy4 default]
	xcMenuEntry $f12 "Specify destination blending type:" 30 \
		openGL(dst_blend) $openGL(dst_blend_list) \
		-labelwidth 30 -labelanchor w
	button $f12.b1 -text "Default" \
		-command [list GetOpenGLPar _DST_${item}_${type}_ \
		openGL(dst_blend) dummy2 dummy3 dummy4 default]
	pack $f1 -side top -fill both -expand 1
	pack $f11 $f12 -side top -fill x -padx 2 -pady 10 -ipady 3
	pack $f11.b1 $f12.b1 -side left -padx 5 -pady 5
    }

    set side FRONT
    if { [string match "*_BACK_*" $item] } { set side BACK }

    # here COLOR means Material-Properties; FRONT_SIDE
    if { [string match "*_COLOR_*" $item] || $test == "test" } {
	set type [string trim $type _]
	set item [string trim $item _]
	xcDebug "3333333>"
	set multi_widget_list(post) $f.f2
	set f2  [frame $f.f2 -relief flat]
	set f2t [frame $f.f2.t -relief flat]
	set f21 [frame $f.f2.1 -relief groove -bd 2]
	set f22 [frame $f.f2.2 -relief groove -bd 2]
	set f23 [frame $f.f2.3 -relief groove -bd 2]
	set f24 [frame $f.f2.4 -relief groove -bd 2]	    
	set f25 [frame $f.f2.5 -relief groove -bd 2]
	
	set ltop [label $f2t.l -text "Color parameters for $side SIDE"]
	# GL_AMBIENT
	setRGBAwidget $f21 "Ambient color:" \
		openGL(ambient_R) openGL(ambient_G) \
		openGL(ambient_B) openGL(ambient_A) \
		_${item}_${type}_AMBIENT_
	# GL_DIFFUSE
	setRGBAwidget $f22 "Diffuse color:" \
		openGL(diffuse_R) openGL(diffuse_G) \
		openGL(diffuse_B) openGL(diffuse_A) \
		_${item}_${type}_DIFFUSE_    
	# GL_SPECULAR
	setRGBAwidget $f23 "Specular color:" \
		openGL(specular_R) openGL(specular_G) \
		openGL(specular_B) openGL(specular_A) \
		_${item}_${type}_SPECULAR_
	# GL_EMISSION
	setRGBAwidget $f24 "Emission color:" \
		openGL(emission_R) openGL(emission_G) \
		openGL(emission_B) openGL(emission_A) \
		_${item}_${type}_EMISSION_
	# GL_SHININESS
	set l   [label $f25.l -text "Shininess:"]
	set sc  [scale $f25.sc -from 0 -to 128 \
		     -length [expr 100 * $xcMisc(resolution_ratio2)] \
		     -variable openGL(shininess) -orient horizontal \
		     -label "specular exponent:" -resolution 1 \
		     -showvalue true -width 7 -sliderlength 20 \
		     -font $xcFonts(small)]
	set def [button $f25.b -text "Default" -font $xcFonts(small) \
		-command \
		[list GetOpenGLPar _${item}_${type}_SHININESS_ \
		openGL(shininess) dummy2 dummy3 dummy4 default]]
	
	pack $f2 $f2t -side top -fill both -expand 1
	pack $ltop -side top -expand 1
	pack $f21 $f22 $f23 $f24 $f25 -side left -fill both \
		-padx 2 -pady 10 -ipady 3 -expand 1
	pack $l -side top -fill x 
	pack $sc -side top -fill y -expand 1 
	pack $def -side top -pady 3
	xcDebug "44444444>"
 	if { $test == "test" } {
	    # findout the geometry of $f2 window
	    tkwait visibility $f2
	    set wid [winfo width  $f2]
	    set hig [winfo height $f2]
	    set hlt [$f2 cget -highlightthickness]
	    $f config -width  [expr $wid + 2 * $hlt + 6]
	    $f config -height [expr $hig + 2 * $hlt + 6]
	    pack propagate $f false
	}
    }
}



proc SetRGBA {vars values} {
    global openGL

    xcDebug "SetRGBA_values> $values"
    set n 0
    foreach item $vars {
	set val [lindex $values $n]
	xcDebug "SetRGBA_val> $val"
	set openGL(${item}_R) [lindex $val 0]
	set openGL(${item}_G) [lindex $val 1]
	set openGL(${item}_B) [lindex $val 2]
	set openGL(${item}_A) [lindex $val 3]
	xcDebug "openGL:: openGL(${item}_R) == $openGL(${item}_R)"
	incr n
    }
}



proc setRGBAwidget {parent l_text varR varG varB varA what} {
    global xcFonts xcMisc
    upvar $varR R
    upvar $varG G
    upvar $varB B
    upvar $varA A

    # set current color:
    set color [SetRGBColor $R $G $B] 
    set color_frame  [frame $parent.cf -width 20 -height 20 \
	    -relief sunken -bd 2 \
	    -background $color]
    
    set l   [label $parent.l -text $l_text]

    set scR [scale $parent.scR -from 0 -to 1 \
		 -length [expr 100 * $xcMisc(resolution_ratio2)] \
		 -variable $varR -orient horizontal \
		 -label "red component:" -resolution 0.01 \
		 -showvalue true -width 7 -sliderlength 20 \
		 -font $xcFonts(small) \
		 -command [list SetColorFrame $color_frame $varR $varG $varB]]    
    
    set scG [scale $parent.scG -from 0 -to 1 \
		 -length [expr 100 * $xcMisc(resolution_ratio2)] \
		 -variable $varG -orient horizontal \
		 -label "green component:" -resolution 0.01 \
		 -showvalue true -width 7 -sliderlength 20 \
		 -font $xcFonts(small) \
		 -command [list SetColorFrame $color_frame $varR $varG $varB]]
    set scB [scale $parent.scB -from 0 -to 1 \
		 -length [expr 100 * $xcMisc(resolution_ratio2)] \
		 -variable $varB -orient horizontal \
		 -label "blue component:" -resolution 0.01 \
		 -showvalue true -width 7 -sliderlength 20 \
		 -font $xcFonts(small) \
		 -command [list SetColorFrame $color_frame $varR $varG $varB]]
    set scA [scale $parent.scA -from 0 -to 1 \
		 -length [expr 100 * $xcMisc(resolution_ratio2)] \
		 -variable $varA -orient horizontal \
		 -label "alpha component:" -resolution 0.01 \
		 -showvalue true -width 7 -sliderlength 20 \
		 -font $xcFonts(small)]

    set def [button $parent.b -text "Default" -font $xcFonts(small) \
	    -command [list GetOpenGLPar $what $varR $varG $varB $varA \
	    default $color_frame]]

    pack $l -side top -fill x
    pack $scR $scG $scB $scA -side top -fill y -expand 1
    pack $color_frame $def -side left -pady 3 -expand 1
    
    if { $what == "_UNKNOWN_" } {
	pack forget $def
	pack configure $color_frame -fill x -padx 3
    }
    
    ### t.k: temporarily
    #if { $l_text == "Emission color:" } {
    #	xcDisableAll $parent
    #}
    ###
}



proc SetColorFrame {f R G B {val {}}} {
    # this proc is also used by scale widget and this require aditional 
    # parameter, thatwhy $val is needed
    upvar $R r
    upvar $G g
    upvar $B b
        
    if { ! [info exists r] } { return }; # BUG work-around

    set color [SetRGBColor $r $g $b] 
    $f configure -background $color
}



proc SetRGBColor {r g b} {
    set r [d2h [expr int($r * 255)]]
    set g [d2h [expr int($g * 255)]]
    set b [d2h [expr int($b * 255)]]
    xcDebug "RGB: #$r$g$b"
    return "#$r$g$b"
}



proc UpdateOpenGLPar {{t {}}} {
    global openGL

    puts stderr "UpdateOpenGLPar: $openGL(type), [info exists openGL(specular_R)]"

    if { $openGL(type) == "material" && [info exists openGL(specular_R)] } {
	########################################
	# first update COLOR parameters
	set com [concat xc_setGLparam material -what $openGL(what)]
	lappend specular [list $openGL(specular_R) $openGL(specular_G) \
		$openGL(specular_B) $openGL(specular_A)]
	lappend ambient  [list $openGL(ambient_R) $openGL(ambient_G)   \
		$openGL(ambient_B) $openGL(ambient_A)]
	lappend diffuse  [list $openGL(diffuse_R) $openGL(diffuse_G)   \
		$openGL(diffuse_B) $openGL(diffuse_A)] 
	lappend emission [list $openGL(emission_R) $openGL(emission_G)   \
		$openGL(emission_B) $openGL(emission_A)]
	eval $com \
	    -shininess $openGL(shininess) \
	    -specular  $specular \
	    -ambient   $ambient  \
	    -diffuse   $diffuse  \
	    -emission  $emission
	    	
	# now render the changes
	.mesa render
    } 

    if { $openGL(type) == "blendfunc" && [info exists openGL(src_blend)] } {
	########################################
	# update blendfunc
	xcDebug "\nxc_setGLparam blendfunc \
		-what  $openGL(what) \
		-sfunc $openGL(src_blend) \
		-dfunc $openGL(dst_blend)\n"
	xc_setGLparam blendfunc \
		-what  $openGL(what) \
		-sfunc $openGL(src_blend) \
		-dfunc $openGL(dst_blend)
	# now render the changes
	.mesa render
    }
}
   


proc GetOpenGLPar {item var1 var2 var3 var4 {default {}} {color_frame {}}} {
    upvar #0 $var1 v1
    upvar #0 $var2 v2
    upvar #0 $var3 v3
    upvar #0 $var4 v4

    if { $default == "default" } { 
	set get "-get def_" 
    } else {
	set get "-get "
    }

    if [string match "*_COLOR_*" $item] {
	set type "material"
	if [string match "*_STRUCTURE_*" $item] {
	    append what "-what structure "
	} elseif [string match "*_ISOSURF_*" $item] {
	    append what "-what isosurf "
	    if [string match "*_ONE_*" $item] {
		append what "one "
	    }
	    if [string match "*_POS_*" $item] {
		append what "pos "
	    }
	    if [string match "*_NEG_*" $item] {
		append what "neg "
	    }
	    if [string match "*_FRONT_*" $item] {
		append what "front "
	    }
	    if [string match "*_BACK_*" $item] {
		append what "back "
	    }
	}
	switch -glob -- $item {
	    *_AMBIENT_*      { append get "ambient" }
	    *_DIFFUSE_*      { append get "diffuse" }
	    *_EMISSION_*     { append get "emission" }
	    *_SPECULAR_*     { append get "specular" }
	    *_SHININESS_*    { append get "shininess" }
	}	
    }

    if [string match {*_LIGHT[0-9]_*} $item] {
	set type "light"
	regexp -- {_LIGHT[0-9]_} $item light
	regexp -- {[0-9]} $light n
	set what "-light $n "
	switch -glob -- $item {	    
	    *_AMBIENT_*        { append get "ambient" }
	    *_DIFFUSE_*        { append get "diffuse" }
	    *_SPECULAR_*       { append get "specular" }
	    *_FRACT-POSITION_* { append get "fract_position" }
	    *_SPOT-DIR_*       { append get "spot_dir" }
	    *_SPOT-EXP_*       { append get "spot_exp" }
	    *_SPOT-CUTOFF_*    { append get "spot_cutoff" }
	    *_CONST-ATTEN_*    { append get "const_atten" }
	    *_LIN-ATTEN_*      { append get "lin_atten" }
	    *_QUAD-ATTEN_*     { append get "quad_atten" }
	}
    }

    if [string match "*_LIGHTMODEL_*" $item] {
	set type "lightmodel"
	set what {}
	switch -glob -- $item {
	    *_TWO-SIDE_*         { append get "two_side" }
	    *_TWO-SIDE-ISO_*     { append get "two_side_iso" }
	    *_AMBIENT_*          { append get "ambient" }
	    *_LOCAL-VIEWER_*     { append get "local_viewer" }
	    *_IS-LIGHT-ENABLED*  { 
		regexp -- {_IS-LIGHT-ENABLED[0-9]_} $item light
		regexp -- {[0-9]} $light n
		append get [list [list is_light_enabled $n]]
	    }
	}	
    }

    if [string match "*_BLEND_*" $item] {
	# this has to be comleted
	if [string match "*_ISOSURF_*" $item] {
	    set bldf [xc_getGLparam blendfunc -what isosurf]
	    xcDebug "bldf:: $bldf"
	}
	if { [string match "*_SRC_*" $item] == 1 && $default == "default" } {
	    set v1 [lindex $bldf 2]
	    return $v1
	} 
	if [string match "*_SRC_*" $item] {
	    set v1 [lindex $bldf 0]
	    return $v1
	}
	if {  [string match "*_DST_*" $item] == 1 && $default == "default" } {
	    set v1 [lindex $bldf 3]
	    return $v1
	}
	if [string match "*_DST_*" $item] {
	    set v1 [lindex $bldf 1]
	    return $v1
	}
    }
    
    if { [string match "*_FOG_*" $item] } {
	set type "fog"
	set what ""
	switch -glob -- $item {	    
	    *_FOGCOLOR_*   { append get "fogcolor" }
	}
    }

    #puts stderr "\nopenGLpar>  xc_getGLparam $type $what $get\n"

    set com    [concat xc_getGLparam $type $what $get]
    set result [eval $com]
    #xcDebug "GetOpenGLPar#2:: result = $result"
    set n_res  [llength $result]
    for {set i 1} {$i <= $n_res} {incr i} {
	set ii [expr $i - 1]	
	if { $i == 1 } { set v1 [lindex $result $ii] }
	if { $i == 2 } { set v2 [lindex $result $ii] }
	if { $i == 3 } { set v3 [lindex $result $ii] }
	if { $i == 4 } { set v4 [lindex $result $ii] }
    }

    if { $color_frame != {} } { 
	SetColorFrame $color_frame v1 v2 v3
    }
    return $result
}

