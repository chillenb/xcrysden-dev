#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/glModParam.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc glModParam {} {
    global glModParam

    #
    # this proc is used for setting: 
    #
    # 1.) material-properties for atoms
    # 2.) fog parameters
    # 3.) Lighting-On AntiAlias parameters
    #

    set t .glpar    
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Parameters Setting" "Setting"

    set nb [NoteBook $t.nb]
    pack $nb -expand 1 -fill both


    #
    # Material-properties for atoms (specular color and shininess)
    #
    $nb insert 0 materials -text "Material Properties"
    set page [$nb getframe materials]
    glModParam:Materials $t $page

    $nb insert 1 depthcuing -text "Depth-Cuing Properties"
    set page [$nb getframe depthcuing]
    glModParam:DepthCuing $t $page

    $nb insert 2 antialias -text "AntiAlias Properties"
    set page [$nb getframe antialias]
    glModParam:AntiAlias $t $page

    $nb raise materials
}


# ------------------------------------------------------------------------
#
# Atomic MATERIALS properties
#
# ------------------------------------------------------------------------

proc glModParam:Materials {t page} {
    global glModParam xcFonts mody
    
    set sf [frame $page.stress -class StressText -relief groove -bd 2]
    set m  [message $sf.m -justify left -anchor center -width 300 \
		-text {Here you can set the atomic specular, emission and shininess material properties.}]
    pack $sf -padx 3 -pady 3 -fill x
    pack $m -side top -expand 1 -fill x -padx 0 -pady 0 -ipadx 3 -ipady 3
    
    set specular  [xc_getGLparam material -what structure -get specular]
    set emission  [xc_getGLparam material -what structure -get emission]
    set shininess [xc_getGLparam material -what structure -get shininess]

    set glModParam(material,emission_R) [lindex $emission 0]
    set glModParam(material,emission_G) [lindex $emission 1]
    set glModParam(material,emission_B) [lindex $emission 2]
    set glModParam(material,emission_A) [lindex $specular 3]

    set glModParam(material,specular_R) [lindex $specular 0]
    set glModParam(material,specular_G) [lindex $specular 1]
    set glModParam(material,specular_B) [lindex $specular 2]
    set glModParam(material,specular_A) [lindex $specular 3]

    set glModParam(material,shininess)  $shininess

    set f_color [frame $page.f]
    pack $f_color -side top -expand 1 -fill both -padx 3 -pady 3

    # EMISSION & SPECULAR colors
    foreach type {emission specular} {
	set frame($type)  [frame $f_color.$type -relief groove -bd 2]
	pack $frame($type) -side left -fill both -padx 2 -pady 0 -ipady 3 -expand 1	    

	setRGBAwidget $frame($type) "[string totitle $type] color:" \
	    glModParam(material,${type}_R) glModParam(material,${type}_G) \
	    glModParam(material,${type}_B) glModParam(material,${type}_A) \
	    _COLOR_STRUCTURE_[string toupper $type]_
    }

    # GL_SHININESS
    set frame(shininess)  [frame $f_color.shininess -relief groove -bd 2]
    pack $frame(shininess) -side left -fill both -padx 2 -pady 0 -ipady 3 -expand 1	    

    set l   [label $frame(shininess).l -text "Shininess:"]
    set sc  [scale $frame(shininess).sc -from 0 -to 128 -length 100 \
		 -variable glModParam(material,shininess) -orient horizontal \
		 -label "specular exponent:" -resolution 1 \
		 -showvalue true -width 7 -sliderlength 20 \
		 -font $xcFonts(small)]    
    set def [button $frame(shininess).b -text "Default" -font $xcFonts(small) \
		 -command \
		 [list GetOpenGLPar _COLOR_STRUCTURE_SHININESS_ \
		      glModParam(material,shininess) dummy2 dummy3 dummy4 default]]    
    pack $l  -side top -fill x 
    pack $sc  -side top -fill y -expand 1 
    pack $def -side top -pady 3

    # AMBIENT_BY_DIFFUSE factor:
    set glModParam(material,ambient_by_diffuse) [xc_getvalue $mody(D_AMBIENT_BY_DIFFUSE)]
    
    set frame(ambient_by_diffuse) [frame $page.abd -relief groove -bd 2]
    pack $frame(ambient_by_diffuse) -side top -fill x -padx 2 -pady 0 -ipady 3 -expand 1	    
    
    set sc [scale $frame(ambient_by_diffuse).sc -from 0 -to 2.0 -length 300 \
		-variable glModParam(material,ambient_by_diffuse) -orient horizontal \
		-label "Ambient/Diffuse factor:" -resolution 0.01 \
		-tickinterval 0.25 -digits 3 \
		-showvalue true -width 7 -sliderlength 20]
    pack $sc -side top -fill x -padx 5 -pady 5        

    #
    # in bottom frame goes the "Close|Update" buttons
    #
    set bottom [frame $page.bottom]
    pack $bottom -side top -fill x -padx 3 -pady 3
    
    set close   [button $bottom.close  -text "Close"  -command [list CancelProc $t]]
    set update  [button $bottom.update -text "Update" -command [list glModParam:Materials:Update]]
    pack $close $update -side left -padx 3 -pady 3 -ipadx 3 -ipady 3 -expand 1	
}

proc glModParam:Materials:Update {{t {}}} {
    global glModParam mody
        
    set emission [list $glModParam(material,emission_R) $glModParam(material,emission_G) \
		      $glModParam(material,emission_B) $glModParam(material,emission_A)]
    set specular [list $glModParam(material,specular_R) $glModParam(material,specular_G) \
		      $glModParam(material,specular_B) $glModParam(material,specular_A)]

    xc_setGLparam material -what structure \
	-shininess $glModParam(material,shininess) \
	-specular  $specular \
	-emission  $emission

    xc_newvalue .mesa $mody(L_AMBIENT_BY_DIFFUSE) $glModParam(material,ambient_by_diffuse)

    # now render the changes
    .mesa render     
}




# ------------------------------------------------------------------------
#
# DEPTH-CUING properties
#
# ------------------------------------------------------------------------

proc glModParam:DepthCuing {t page {togl .mesa}} {
    global glModParam xcFonts mody

    # this must be the same as in struct.h
    set XC_FOG_BGCOLOR     0
    set XC_FOG_CUSTOMCOLOR 1
    set XC_FOG_LINEAR      2
    set XC_FOG_EXP         3
    set XC_FOG_EXP2        4

    set glModParam(fog,mode)	    [xc_getvalue $mody(GET_FOG_MODE)]        
    set glModParam(fog,colormode)   [xc_getvalue $mody(GET_FOG_COLORMODE)]   
    set glModParam(fog,density)	    [xc_getvalue $mody(GET_FOG_DENSITY)]     
    set glModParam(fog,ort_start_f) [xc_getvalue $mody(GET_FOG_ORT_START_F)] 
    set glModParam(fog,ort_end_f)   [xc_getvalue $mody(GET_FOG_ORT_END_F)]   
    set glModParam(fog,persp_f1)    [xc_getvalue $mody(GET_FOG_PERSP_F1)]    
    set glModParam(fog,persp_f2)    [xc_getvalue $mody(GET_FOG_PERSP_F2)]    

    #puts stderr "DEBUG> glModParam(fog,colormode) == $glModParam(fog,colormode)"
    #puts stderr "DEBUG> glModParam(fog,mode) == $glModParam(fog,mode)"

    if { $glModParam(fog,colormode) == $XC_FOG_BGCOLOR } { 
	set glModParam(fog,colormode_txt) "use background color" 
    } else {
	set glModParam(fog,colormode_txt) "use custom color"
    }
    if { $glModParam(fog,mode) == $XC_FOG_LINEAR } { 
	set glModParam(fog,mode_txt) "Linear" 
    } elseif { $glModParam(fog,mode) == $XC_FOG_EXP } {
	set glModParam(fog,mode_txt) "Exponential"
    } else {
	set glModParam(fog,mode_txt) "Gaussian"
    }

    foreach f {1 2} {
	set frame($f) [frame $page.$f -relief groove -bd 2]
	pack $frame($f) -side top -expand 1 -fill both -padx 3 -pady 3
    }

    set f1l  [frame $frame(1).left]
    set f1r  [frame $frame(1).right -relief groove -bd 2]
    pack $f1l -side left  -padx 3 -pady 3 -fill both
    pack $f1r -side right -padx 6 -pady 8 -fill y

    set sf [frame $f1l.stress -class StressText -relief groove -bd 2]
    set m  [message $sf.m -justify left -anchor center -width 200 \
		-text {Depth-Cuing, also known as fog, make objects fade into the distance.}]
    pack $sf -padx 5 -pady 5 -fill x
    pack $m -side top -expand 1 -fill x -padx 0 -pady 0 -ipadx 3 -ipady 3

    set f1l1 [frame $f1l.1]
    set f1l2 [frame $f1l.2]
    pack $f1l1 $f1l2 -side top -fill x -pady 3 -padx 3
    
    RadioBut $f1l1 "Fog color-mode:" glModParam(fog,colormode_txt) top top 1 0 \
	"use background color" "use custom color"
    
    RadioBut $f1l2 "Fog mode:" glModParam(fog,mode_txt) top top 1 0 \
	"Linear" "Exponential" "Gaussian"     

    trace variable glModParam(fog,colormode_txt) w glModParam:DepthCuing:_widget
    trace variable glModParam(fog,mode_txt)      w glModParam:DepthCuing:_widget

    #set m1 [xcMenuButton $f1l1 \
    #		-labeltext   "Fog color mode:" \
    #		-labelwidth   15 \
    #		-side         top \
    #		-textvariable glModParam(fog,colormode) \
    #		-menu {
    #		    "use background color" {glModParam:DepthCuing:Set colormode XC_FOG_BGCOLOR}
    #		    "use custom color"     {glModParam:DepthCuing:Set colormode XC_FOG_CUSTOMCOLOR}
    #		}]
    #set m2 [xcMenuButton $f1l2 \
    #		-labeltext   "Fog mode:" \
    #		-labelwidth   15 \
    #		-textvariable glModParam(fog,mode) \
    #		-side         top \
    #		-menu {
    #		    "Linear"      {glModParam:DepthCuing:Set mode XC_FOG_LINEAR}
    #		    "Exponential" {glModParam:DepthCuing:Set mode XC_FOG_EXP}
    #		    "Gaussian"    {glModParam:DepthCuing:Set mode XC_FOG_EXP2}
    #		}]
    #pack $m1 $m2 -side top -padx 3 -expand 1 -fill x
    
    set fogcolor [xc_getGLparam fog -get fogcolor]
    set glModParam(fog,color_R) [lindex $fogcolor 0]
    set glModParam(fog,color_G) [lindex $fogcolor 1]
    set glModParam(fog,color_B) [lindex $fogcolor 2]
    set glModParam(fog,color_A) [lindex $fogcolor 3]
    
    set glModParam(fog,fogcolor_framewidget) $f1r
    setRGBAwidget $f1r "Fog color:" \
	glModParam(fog,color_R) glModParam(fog,color_G) \
	glModParam(fog,color_B) glModParam(fog,color_A) \
	_FOG_FOGCOLOR_
    
    set glModParam(fog,density_entrywidget) [FillEntries $frame(2) {
	"Fog DENSITY:" 
	"Fog START factor (for Orthographic projection):"
	"Fog END   factor (for Orthographic projection):"
	"Fog NEAR  factor (for Perspective projection):"
	"Fog FAR   factor (for Perspective projection):"
    } {
	glModParam(fog,density) 
	glModParam(fog,ort_start_f) glModParam(fog,ort_end_f) 
	glModParam(fog,persp_f1)    glModParam(fog,persp_f2) 
    } 40 10]
    #puts stderr "DEBUG> fog: [winfo exists $glModParam(fog,density_entrywidget)]"

    #
    # in bottom frame goes the "Close|Update" buttons
    #
    set bottom [frame $page.bottom]
    pack $bottom -side top -fill x -padx 3 -pady 3
    
    set close   [button $bottom.close  -text "Close"  -command [list CancelProc $t]]
    set update  [button $bottom.update -text "Update" -command [list glModParam:DepthCuing:Update $togl]]
    pack $close $update -side left -padx 3 -pady 3 -ipadx 5 -ipady 5 -expand 1	

    glModParam:DepthCuing:_widget glModParam fog,colormode_txt init
    glModParam:DepthCuing:_widget glModParam fog,mode_txt      init
}

proc glModParam:DepthCuing:_widget {name1 name2 op} {
    global glModParam

    regsub -- {1.entry1$} $glModParam(fog,density_entrywidget) {} core
    set start_entrywidget ${core}2.entry2
    set end_entrywidget   ${core}3.entry3
    set near_entrywidget  ${core}4.entry4
    set far_entrywidget   ${core}5.entry5

    if { $name2 == "fog,colormode_txt" } {
	if { [winfo exists $glModParam(fog,fogcolor_framewidget)] } {
	    if { $glModParam(fog,colormode_txt) == "use background color" } {
		xcDisableAll -disabledfg $glModParam(fog,fogcolor_framewidget)
	    } else {
		# "use custom color"
		xcEnableAll -disabledfg $glModParam(fog,fogcolor_framewidget)
	    }
	}
    } else {
	# *_entrywidget 	
	if { [winfo exists $glModParam(fog,density_entrywidget)] } {
	    if { $glModParam(fog,mode_txt) == "Linear" } {
		xcDisableEntry $glModParam(fog,density_entrywidget)
		xcEnableEntry  $start_entrywidget $end_entrywidget \
		    $near_entrywidget $far_entrywidget
	    } else {
		xcEnableEntry  $glModParam(fog,density_entrywidget)
		xcDisableEntry $start_entrywidget $end_entrywidget \
		    $near_entrywidget $far_entrywidget
	    }
	}
    }
}

proc glModParam:DepthCuing:Update {togl} {
    global glModParam mody

    # this must be the same as in struct.h
    set XC_FOG_BGCOLOR     0
    set XC_FOG_CUSTOMCOLOR 1
    set XC_FOG_LINEAR      2
    set XC_FOG_EXP         3
    set XC_FOG_EXP2        4

    #
    # FOG-COLOR mode
    #
    if { $glModParam(fog,colormode_txt) == "use background color" } {
	set glModParam(fog,colormode) $XC_FOG_BGCOLOR
	xcDisableAll $glModParam(fog,fogcolor_framewidget)
    } else {
	# "use custom color"
	set glModParam(fog,colormode) $XC_FOG_CUSTOMCOLOR
	set fogcolor [list $glModParam(fog,color_R) $glModParam(fog,color_G) \
			  $glModParam(fog,color_B) $glModParam(fog,color_A)]
	xcEnableAll $glModParam(fog,fogcolor_framewidget)
	
	xc_setGLparam fog -color $fogcolor
    }
    
    #
    # FOG-MODE
    #
    if { $glModParam(fog,mode_txt) == "Linear" } {
	set glModParam(fog,mode) $XC_FOG_LINEAR
    } elseif { $glModParam(fog,mode_txt) ==  "Exponential" } {
	set glModParam(fog,mode) $XC_FOG_EXP
    } else {
	# Gaussian
	set glModParam(fog,mode) $XC_FOG_EXP2
    }   
 
    if { ! [string is double $glModParam(fog,density) ] } {
	WarningDialog "wanted double, but got $glModParam(fog,density), for \"Fog DENSITY\""
	return
    }
    if { ! [string is double $glModParam(fog,ort_start_f) ] } {
	WarningDialog "wanted double, but got $glModParam(fog,ort_start_f), for \"Fog START factor\""
	return
    }
    if { ! [string is double $glModParam(fog,ort_end_f) ] } {
	WarningDialog "wanted double, but got $glModParam(fog,ort_end_f), for \"Fog END factor\""
	return
    }
    if { ! [string is double $glModParam(fog,persp_f1) ] } {
	WarningDialog "wanted double, but got $glModParam(fog,persp_f1), for \"Fog NEAR factor\""
	return
    }
    if { ! [string is double $glModParam(fog,persp_f2) ] } {
	WarningDialog "wanted double, but got $glModParam(fog,persp_f2), for \"Fog FAR factor\""
	return
    }
    
    xc_newvalue $togl $mody(SET_FOG_MODE)        $glModParam(fog,mode)
    xc_newvalue $togl $mody(SET_FOG_COLORMODE)   $glModParam(fog,colormode)
    xc_newvalue $togl $mody(SET_FOG_DENSITY)     $glModParam(fog,density)
    xc_newvalue $togl $mody(SET_FOG_ORT_START_F) $glModParam(fog,ort_start_f)
    xc_newvalue $togl $mody(SET_FOG_ORT_END_F)   $glModParam(fog,ort_end_f)
    xc_newvalue $togl $mody(SET_FOG_PERSP_F1)    $glModParam(fog,persp_f1)
    xc_newvalue $togl $mody(SET_FOG_PERSP_F2)    $glModParam(fog,persp_f2)

    # now render the changes
    $togl render     
}

# ------------------------------------------------------------------------
#
# ANTI-ALIAS properties
#
# ------------------------------------------------------------------------

proc glModParam:AntiAlias {t page {togl .mesa}} {
    global glModParam mody
    
    set glModParam(antialias,degree) [xc_getvalue $mody(GET_ANTIALIAS_DEGREE)]        
    set glModParam(antialias,offset) [xc_getvalue $mody(GET_ANTIALIAS_OFFSET)]

    set sf [frame $page.stress -class StressText -relief groove -bd 2]
    set m  [message $sf.m -justify left -anchor center -width 300 \
		-text {Here you can set the two parameters (degree & offset) for multi-sampling anti-aliasing, which is used for Lighting-On display mode. The degree of anti-aliasing determine the number of samplings (degree=1 --> 9 samplings; degree=2 --> 25 samplings, degree=3 --> 36 samplings, etc). The larger the degree the more time consuming the anti-aliasing. The multi-sampling is very time consuming, so its typical use is for printing.}]
    pack $sf -padx 3 -pady 3 -fill x
    pack $m -side top -expand 1 -fill x -padx 0 -pady 0 -ipadx 3 -ipady 3

    foreach i {1 2} {
	set f($i) [frame $page.f$i -relief groove -bd 2]
	pack $f($i) -side top -fill x -padx 3 -pady 5
    }
    
    set sc1 [scale $f(1).sc -from 1 -to 5 -length 300 \
		 -variable glModParam(antialias,degree) -orient horizontal \
		 -label "Degree of multi-sample Antialiasing:" -resolution 1 \
		 -tickinterval 1 -digits 1 \
		 -showvalue true -width 7 -sliderlength 20]

    set sc2 [scale $f(2).sc -from 0.01 -to 1.99 -length 300 \
		 -variable glModParam(antialias,offset) -orient horizontal \
		 -label "Offset for multi-sample Antialiasing:" -resolution 0.01 \
		 -tickinterval 0.25 -digits 3 \
		 -showvalue true -width 7 -sliderlength 20]
    pack $sc1 $sc2 -side top -fill y -expand 1 -padx 5 -pady 5
    
    #
    # in bottom frame goes the "Close|Update" buttons
    #
    set bottom [frame $page.bottom]
    pack $bottom -side bottom -fill x -padx 3 -pady 3
       
    set close   [button $bottom.close  -text "Close"  -command [list CancelProc $t]]
    set update  [button $bottom.update -text "Update" -command [list glModParam:AntiAlias:Update $togl]]
    pack $close $update -side left -padx 3 -pady 3 -ipadx 5 -ipady 5 -expand 1	
}
proc glModParam:AntiAlias:Update {togl} {
    global glModParam mody
    
    xc_newvalue $togl $mody(SET_ANTIALIAS_DEGREE) $glModParam(antialias,degree)
    xc_newvalue $togl $mody(SET_ANTIALIAS_OFFSET) $glModParam(antialias,offset)

    # now render the changes
    $togl render     
}
