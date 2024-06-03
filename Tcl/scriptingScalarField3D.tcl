#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scriptingScalarField3D.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


# ------------------------------------------------------------------------
#****c* Scripting/scripting::scalarField3D
#
# NAME
# scripting::scalarField3D
#
# PURPOSE

# Encapsulate the scripting interface for manipulating the 3D scalar
# field and display of isosurfaces and colorplane+contour plots. A
# typical usage of the interface is the following. First the scalar
# field is loaded with the "load" command. Then with the "configure"
# command the isosurface and/or contours display parameters are set.
# Finally the isosurface and/or contours are displayed using the
# "render" command. It is possible to make several snapshots of
# isosurface and/or contours by changing configuration options, this
# would look as:
#
# load
# configure ...
# render
# configure ..
# render
# ...

#
# COMMANDS

#
# -- scripting::scalarField3D::load

# Loads the scalar field.

#
# -- scripting::scalarField3D::configure

# Configure the display parameters for isosurface and/or contours
# plot.

#
# -- scripting::scalarField3D::render

# Renders the sosurface and/or contours.

#****
# ------------------------------------------------------------------------

namespace eval scripting::scalarField3D {
    variable scalarField3D
    
    set scalarField3D(loaded)     0
    set scalarField3D(configured) 0
    set scalarField3D(basalplane) -1
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField3D::load
#
# NAME
# scripting::scalarField3D::load
#
# USAGE
# scripting::scalarField3D::load
#
# PURPOSE

# This proc loads the 3D scalar field. Prior to calling this routine
# the XSF file with the description of 3D scalar filed should be
# loaded by "scripting::open --xsf file.xsf" call.

#
# RETURN VALUE
# Undefined.
#
# WARNINGS
# Prior to calling this routine the XSF file with the description of
# 3D scalar filed should be loaded by "scripting::open --xsf file.xsf"
# call.
#
# EXAMPLE
# scripting::open --xsf file.xsf
# scripting::scalarField3D::load
#
#****
# ------------------------------------------------------------------------

proc scripting::scalarField3D::load {} {
    variable scalarField3D
    
    #
    # load datagrid 
    #    
    set dim [DataGridOK [DataGrid]]
    if { $dim != "3D" } {
	error "datagrid is not 3D, but $dim"
    }
    set scalarField3D(loaded) 1

    # load the scalar-field 3D defaults

    IsoControl_InitVar
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField3D::configure
#
# NAME
# scripting::scalarField3D::configure
#
# USAGE
# scripting::scalarField3D::configure -option value ?-option value? ...
#
# PURPOSE

# This proc configures the display parameters of isosurface and/or
# contours plot.

# ARGUMENTS
# args -- various configuration "-option value" pairs (see Options).
#
# OPTIONS
# ------------------------------------------------------------------------
#   OPTION::                 ALLOWED-VALUES + Description        
# ------------------------------------------------------------------------
#   -isosurface              0|1                                 
#                            0 = do not render isosurface
#                            1 = render isosurface
#
#   -interpolation_degree    integer                             
#                            degree of scalar-field interpolation
#
#   -isolevel                real                                
#                            the isovalue of isosurface
#
#   -plusminus               0|1
#                            0 = display only the isosurface of isovalue 
#                                specified by -isolevel
#                            1 = display the two isosurfaces of +-isovalue
#
#   -revertsides             pos|neg|{pos neg}
#                            pos = revert the front- and back-side of isovalue 
#                                  isosurface
#                            neg = revert the front- and back-side of -isovalue
#                                  isosurface
#   -revertnormals           pos|neg|{pos neg}
#                            pos = revert the normals of isovalue isosurface
#                            neg = revert the normals of -isovalue isosurface
#
#   -expand                  none|whole|specify
#                            none = do not expand the isosurface along the 
#                                   periodic directions
#                            whole = expand the isosurface over the
#                                    whole structure along the periodic 
#                                    directions
#                            specify = expand the isosurface along the
#                                      periodic directions as specified by 
#                                      -expand_X, -expand_Y, and -expand_Z 
#                                      factors
#
#   -expand_X                positive-integer
#                            expand  isosurface n-times along the 1st
#                            periodic dimension
#
#
#   -expand_Y                positive-integer
#                            expand  isosurface n-times along the 2nd
#                            periodic dimension
#
#   -expand_Z                positive-integer
#                            expand  isosurface n-times along the 3rd
#                            periodic dimension
#
#   -basalplane              1|2|3
#                            show the ith basal plane (1=xy, 2=xz, 3=yz)
#                            as colorplane and/or isolines (as specified by
#                            -colorplane and -isoline options)
#
#   -colorbasis              MONOCHROME|RAINBOW|RGB|GEOGRAPHIC|BLUE-WHITE-RED|BLACK-BROWN-WHITE
#                            the color basis for the colorplane
#
#   -scalefunction           LINEAR|LOG10|SQRT|3th-ROOT|EXP(x)|EXP(x^2)
#                            the scalefunctions for contour/colorplane plots
#
#   -expand2D                none|whole|specify
#                            none = do not expand the contour/colorplane plots
#                                   along the periodic directions
#                            whole = expand the contour/colorplane plots
#                                    over the whole structure along the 
#                                    periodic directions
#                            specify = expand the contour/colorplane plots
#                                    along the periodic directions as specified
#                                    by -expand2D_X, -expand2D_Y, and 
#                                    -expand2D_Z factors
#
#   -expand2D_X              positive-integer
#                            expand  contour/colorplane n-times along the 1st
#                            periodic dimension
#
#
#   -expand2D_Y              positive-integer
#                            expand  contour/colorplane n-times along the 2nd
#                            periodic dimension
#
#   -expand2D_Z              positive-integer
#                            expand  contour/colorplane n-times along the 3rd
#                            periodic dimension
#
#   -colorplane              0|1
#                            do not display|display the colorplane
#
#   -isoline                 0|1
#                            do not display|display the isolines
#             
#   -colorplane_lighting     0|1
#                            0 = do not perform lighting for colorplane
#                            1 = perform lighting for colorplane
#
#   -cpl_transparency        0|1
#                            0 = render colorplane as non-transparent
#                            1 = render colorplane as transparent
#
#   -cpl_thermometer         0|1
#                            0 = do not make a legend (i.e. thermometer) for 
#                                colorplane colors
#                            1 = make a legend (i.e. thermometer) for 
#                                colorplane colors
#
#   -2Dlowvalue              real
#                            minimum rendered value of colorplane/isolines
#
#   -2Dhighvalue             real
#                            maximum rendered value of colorplane/isolines
#              
#   -2Dnisoline              positive-integer
#                            number of isoline
#
#   -anim_step               positive-integer
#                            animation step
#
#   -current_slide           positive-integer
#                            iD (i.e. sequential number) of colorplane
#
#   -isoline_color           monocolor|{property color}
#                            monocolor = all isolines have the same color
#                            {property color} = isolines are colorer according
#                            to color-basis
#
#   -isoline_width           positive-integer
#                            width (i.e. thickness) of isolines
#
#   -isoline_monocolor       #rgb
#                            color of the monolor-type isolines
#
# ------------------------------------------------------------------------
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::scalarField3D::configure \
#     -isosurface           1 \
#     -interpolation_degree 2 \
#     -isolevel             0.1 \
#     -plusminus            1 \
#     -revertsides          {pos neg} \
#     -revertnormal         {pos neg} \
#     -expand               specify \
#     -expand_X 	    1 \
#     -expand_Y             1 \
#     -expand_Z             1 \
#     -basalplane           1 \
#     -colorbasis           BLUE-WHITE-RED \
#     -scalefunction        LINEAR \
#     -expand2D             specify \
#     -expand2D_X 	    1 \
#     -expand2D_Y           1 \
#     -expand2D_Z           1 \
#     -colorplane           1 \
#     -isoline              1 \
#     -colorplane_lighting  0 \
#     -cpl_transparency     0 \
#     -cpl_thermometer      1 \
#     -2Dlowvalue           -0.1 \
#     -2Dhighvalue          +0.1 \
#     -2Dnisoline           11 \
#     -anim_step            1 \
#     -current_slide        30 \
#     -isoline_color        monocolor \
#     -isoline_width        3 \
#     -isoline_monocolor    \#000000
# 
#****
# ------------------------------------------------------------------------

proc scripting::scalarField3D::configure {args} {
    variable scalarField3D
    global isosurf isoControl prop

    #
    # assuming the 3D grid !!!
    #

    set scalarField3D(configured) 1

    # if the colorplane/isolines are to be rendered, then "-basalplane
    # whichplane" should be specigied, and the iplane will be assigned
    # accordingly

    set iplane $scalarField3D(basalplane)

    # parse the options

    set i 0
    foreach option $args {
	incr i
	
	if { $i%2 } {
            set tag $option
	} else {	    
            switch -glob -- $tag {
		"-isosurface" {
		    switch -- $option {
			on - 1  { set isoControl(isosurf)  1 }
			off - 0 { set isoControl(isosurf)  0 }
			default {
			    ErrorDialog "wrong -isosurface value $option, must be 0 or 1"
			}
		    }
		}
		
		"-interpolation_degree" {
		    if { ! [string is integer $option] } {
			ErrorDialog "expected integer for -interpolationdegree option, but got $option"
		    } elseif { $option < 1 } {
			ErrorDialog "-interpolationdegree value should be greater then 0, but got $option"
		    } else {
			set isosurf(3Dinterpl_degree) $option
		    }
		}
		
		"-isolevel" {
		    if { ! [string is double $option] } {
			ErrorDialog "expected real-number for -isolevel option, but got $option"
		    } else {
			set prop(isolevel) $option
		    }
		}

		"-plusminus" {
		    switch -- $option {
			yes - on  - 1 { set prop(pm_isolevel) 1 }
			no  - off - 0 { set prop(pm_isolevel) 0 }
			default {
			    ErrorDialog "wrong -isosurface value $option, must be 0 or 1"
			}
		    }
		    IsoControlCommand
		}
		
		"-revertside*" {
		    foreach side $option {
			switch -- [string tolower $side] {
			    pos { RevertIsoSides pos }
			    neg { RevertIsoSides neg }
			    default {
				ErrorDialog "the value of -revertsides must be pos or neg, but got $option"	
			    }
			}
		    }
		}
		
		"-revertnormal*" {
		    foreach side $option {
			switch -- [string tolower $side] {
			    pos { RevertIsoNormals pos }
			    neg { RevertIsoNormals neg }
			    default {
				ErrorDialog "the value of -revertnormals must be pos or neg, but got $option"	
			    }
			}
		    }
		}

		"-basalplane" {
		    switch $option {
			1 { set iplane 1 }
			2 { set iplane 2 }
			3 { set iplane 3 }
			default {
			    ErrorDialog "the value of -basalplane must be 1, 2, or 3, but got $option" 
			    continue
			}
		    }
		    set scalarField3D(basalplane) $iplane
		    set t .iso
		    set ft  $t.ft
		    set fb1 $t.fb1
		    set fb2 $t.fb2
		    IsoControl_Show colorplane$iplane $fb1 $fb2 $ft.b1 $ft.b2 $ft.b3 $ft.b4
		}

		"-colorbasis" {
		    switch -exact -- $option {
			MONOCHROME - RAINBOW - RGB - GEOGRAPHIC - BLUE-WHITE-RED - BLACK-BROWN-WHITE { set cpl_basis $option }			    
			default {
			    ErrorDialog "wrong colorbasis $option, must be ..." 
			    continue
			}  
		    }
		}

		"-scalefunction" {
		    switch -exact -- $option {
			LINEAR - LOG10 - SQRT - 3th-ROOT - EXP(x) - EXP(x^2) {set cpl_function $option }
			default {			    
			    ErrorDialog "wrong scale-function $option, must be ..." 
			    continue
			}  
		    }
		}

		"-expand" -
		"-expand2D" {
		    set var [string trimleft $tag -]
		    switch -- $option {
			none - whole - specify { set $var $option }
			default {
			    ErrorDialog "wrong $tag value $option, must be none, whole, specify" 
			    continue
			}
		    }
		}
		"-expand_X" -	
		"-expand_Y" -
		"-expand_Z" -
		"-expand2D_X" -	
		"-expand2D_Y" -
		"-expand2D_Z" {
		    if { ! [string is integer $option] } {
			ErrorDialog "expected integer for $tag option, but got $option"
		    } elseif { $option < 1 } {
			ErrorDialog "$tag value should be greater then 0, but got $option"
		    } else {
			set var [string trimleft $tag -]
			set $var $option
		    }
		}
		
		"-colorplane" -
		"-isoline" -
		"-colorplane_lighting" -
		"-cpl_transparency" -
		"-cpl_thermometer" -
		"-2Dlowvalue" -
		"-2Dhighvalue" -
		"-2Dnisoline" -
		"-anim_step" -    
		"-current_slide" -
		"-isoline_color" -
		"-isoline_width" -
		"-isoline_monocolor" {
		    set var [string trimleft $tag -]
		    set $var $option
		}
	    }
	}
    }
    if { $i%2 } {
	ErrorDialog "scripting::scalarField3D::configure called with an odd number of arguments !!!"
    }

    if { [info exists expand]   } { set isosurf(expand)   $expand }
    if { [info exists expand_X] } { set isosurf(expand_X) $expand_X }
    if { [info exists expand_Y] } { set isosurf(expand_Y) $expand_Y }
    if { [info exists expand_Z] } { set isosurf(expand_Z) $expand_Z }
	

    if { $iplane > 0 } {

	if { [info exists cpl_basis] } {           set isoControl($iplane,cpl_basis)           $cpl_basis }
	if { [info exists cpl_function] } {        set isoControl($iplane,cpl_function)        $cpl_function }
	if { [info exists colorplane] } {          set isoControl($iplane,colorplane)          $colorplane }
	if { [info exists isoline] } {             set isoControl($iplane,isoline)             $isoline }
	if { [info exists colorplane_lighting] } { set isoControl($iplane,colorplane_lighting) $colorplane_lighting }
	if { [info exists cpl_transparency] } {    set isoControl($iplane,cpl_transparency)    $cpl_transparency }
	if { [info exists cpl_thermometer] } {     set isoControl($iplane,cpl_thermometer)     $cpl_thermometer }
	if { [info exists cpl_thermoTplw] } {      set isoControl($iplane,cpl_thermoTplw)      $cpl_thermoTplw }
	
	if { [info exists 2Dlowvalue]  } { set isoControl($iplane,2Dlowvalue)  $2Dlowvalue }
	if { [info exists 2Dhighvalue] } { set isoControl($iplane,2Dhighvalue) $2Dhighvalue }
	if { [info exists 2Dnisoline]  } { set isoControl($iplane,2Dnisoline)  $2Dnisoline }
	
	if { [info exists expand2D]   } { set isosurf($iplane,2Dexpand)   $expand2D }
	if { [info exists expand2D_X] } { set isosurf($iplane,2Dexpand_X) $expand2D_X }
	if { [info exists expand2D_Y] } { set isosurf($iplane,2Dexpand_Y) $expand2D_Y }
	if { [info exists expand2D_Z] } { set isosurf($iplane,2Dexpand_Z) $expand2D_Z }
	
	if { [info exists anim_step] } {        
	    set isoControl($iplane,anim_step) $anim_step 
	}
	if { [info exists current_slide] } {    
	    set isoControl($iplane,current_slide) $current_slide 
	    set isoControl(current_text_slide) "Current slide:  $isoControl($iplane,current_slide) / $isoControl($iplane,nslide)"	
	}
	
	if { [info exists isoline_color] } {    set isoControl($iplane,isoline_color)      $isoline_color }
	if { [info exists isoline_width] } {    set isoControl($iplane,isoline_width)      $isoline_width }
	if { [info exists isoline_monocolor] } {set isoControl($iplane,isoline_monocolor)  $isoline_monocolor }
	if { [info exists isoline_stipple] } {  set isoControl($iplane,isoline_stipple)    $isoline_stipple }

	# register the colorplane values
	IsoControl_SetColorPlaneVar 2nd $iplane
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField3D::render
#
# NAME
# scripting::scalarField3D::render
#
# USAGE
# scripting::scalarField3D::render
#
# PURPOSE

# This proc displays the either isosurface or contours/colorplanes or
# both (depending on the configuration).

#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::scalarField3D::render
#
#****
# ------------------------------------------------------------------------

proc scripting::scalarField3D::render {} {
    variable scalarField3D

    if { !$scalarField3D(loaded) } {
	ErrorDialog "can't render isosurface/contours, the scalar field was not loaded. Call scripting::scalarField3D::load and scripting::scalarField3D::configure before calling scripting::scalarField3D::render"
	return
    }

    if { !$scalarField3D(configured) } {
	ErrorDialog "can't render isosurface/contours, since they were not configured. Call scripting::scalarField3D::configure before calling scripting::scalarField3D::render"
	return
    }

    #############
    UpdateIsosurf
}
