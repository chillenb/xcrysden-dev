#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scriptingScalarField2D.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ------------------------------------------------------------------------
#****c* Scripting/scripting::scalarField2D
#
# NAME
# scripting::scalarField2D
#
# PURPOSE

# Encapsulate the scripting interface for manipulating the 2D scalar
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
# -- scripting::scalarField2D::load

# Loads the scalar field.

#
# -- scripting::scalarField2D::configure

# Configure the display parameters for isosurface and/or contours
# plot.

#
# -- scripting::scalarField2D::render

# Renders the sosurface and/or contours.

#****
# ------------------------------------------------------------------------

namespace eval scripting::scalarField2D {
    variable scalarField2D
    
    set scalarField2D(loaded)     0
    set scalarField2D(configured) 0
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField2D::load
#
# NAME
# scripting::scalarField2D::load
#
# USAGE
# scripting::scalarField2D::load
#
# PURPOSE
# This proc load the 2D scalar field.
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::scalarField2D::load
#
#****
# ------------------------------------------------------------------------

proc scripting::scalarField2D::load {} {
    variable scalarField2D
    
    #
    # load datagrid 
    #    
    set dim [DataGridOK [DataGrid]]
    if { $dim != "2D" } {
	error "datagrid is not 2D, but $dim"
    }
    set scalarField2D(loaded) 1

    # 
    # check that the grid is really 2D
    #

    # load the scalar-field 2D defaults
    IsoControl_InitVar
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField2D::configure
#
# NAME
# scripting::scalarField2D::configure
#
# USAGE
# scripting::scalarField2D::configure -option value ?-option value? ...
#
# PURPOSE

# This proc configures the display parameters of colorplane/contours
# plot.

#
# ARGUMENTS
# args -- various configuration "-option value" pairs (see Options).
#
# OPTIONS
# ------------------------------------------------------------------------
#   OPTION::                 ALLOWED-VALUES + Description        
# ------------------------------------------------------------------------
#   -interpolation_degree    integer                            
#                            degree of scalar-field interpolation
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
# scripting::scalarField2D::configure \
#     -interpolation_degree 2 \
#     -colorbasis           RAINBOW \
#     -scalefunction        LOG10 \
#     -expand2D             specify \
#     -expand2D_X 	    1 \
#     -expand2D_Y           1 \
#     -expand2D_Z           1 \
#     -colorplane           1 \
#     -isoline              1 \
#     -colorplane_lighting  0 \
#     -cpl_transparency     0 \
#     -cpl_thermometer      1 \
#     -2Dlowvalue           +0.00001 \
#     -2Dhighvalue          +0.1 \
#     -2Dnisoline           5 \
#     -isoline_color        monocolor \
#     -isoline_width        3 \
#     -isoline_monocolor    \#ffffff
#
#****
# ------------------------------------------------------------------------

proc scripting::scalarField2D::configure {args} {
    variable scalarField2D
    global isosurf isoControl prop

    #
    # assuming the 2D grid !!!
    #

    set scalarField2D(configured) 1

    # parse the options
    
    set i 0
    foreach option $args {
	incr i
	
	if { $i%2 } {
            set tag $option
	} else {	    
            switch -glob -- $tag {
		"-interpolation_degree" {
		    if { ! [string is integer $option] } {
			ErrorDialog "expected integer for -interpolationdegree option, but got $option"
		    } elseif { $option < 1 } {
			ErrorDialog "-interpolationdegree value should be greater then 0, but got $option"
		    } else {
			set isosurf(3Dinterpl_degree) $option
		    }
		}
		
		"-colorbasis" {
		    switch -exact -- $option {
			MONOCHROME - RAINBOW - RGB - GEOGRAPHIC - BLUE-WHITE-RED { set cpl_basis $option }			    
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

		"-expand2D" {
		    switch -- $option {
			none - whole - specify { set expand2D $option }
			default {
			    ErrorDialog "wrong -expand2D value $option, must be none, whole, specify" 
			    continue
			}
		    }
		}
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
	ErrorDialog "scripting::scalarField2D::configure called with an odd number of arguments !!!"
    }
    
    if { [info exists 2Dexpand]   } { set isosurf(2Dexpand)       $2Dexpand }
    if { [info exists 2Dexpand_X] } { set isosurf(2Dexpand_X)     $2Dexpand_X }
    if { [info exists 2Dexpand_Y] } { set isosurf(2Dexpand_Y)     $2Dexpand_Y }
    if { [info exists 2Dexpand_Z] } { set isosurf(2Dexpand_Z)     $2Dexpand_Z }   	

    if { [info exists cpl_basis] } {           set isoControl(cpl_basis)           $cpl_basis }
    if { [info exists cpl_function] } {        set isoControl(cpl_function)        $cpl_function }
    if { [info exists colorplane] } {          set isoControl(colorplane)          $colorplane }
    if { [info exists isoline] } {             set isoControl(isoline)             $isoline }
    if { [info exists colorplane_lighting] } { set isoControl(colorplane_lighting) $colorplane_lighting }
    if { [info exists cpl_transparency] } {    set isoControl(cpl_transparency)    $cpl_transparency }
    if { [info exists cpl_thermometer] } {     set isoControl(cpl_thermometer)     $cpl_thermometer }
    if { [info exists cpl_thermoTplw] } {      set isoControl(cpl_thermoTplw)      $cpl_thermoTplw }    
    if { [info exists 2Dlowvalue] } {          set isoControl(2Dlowvalue)          $2Dlowvalue }
    if { [info exists 2Dhighvalue] } {         set isoControl(2Dhighvalue)         $2Dhighvalue }
    if { [info exists 2Dnisoline] } {          set isoControl(2Dnisoline)          $2Dnisoline }    
    if { [info exists isoline_color] } {       set isoControl(isoline_color)       $isoline_color }
    if { [info exists isoline_width] } {       set isoControl(isoline_width)       $isoline_width }
    if { [info exists isoline_monocolor] } {   set isoControl(isoline_monocolor)   $isoline_monocolor }
    if { [info exists isoline_stipple] } {     set isoControl(isoline_stipple)     $isoline_stipple }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::scalarField2D::render
#
# NAME
# scripting::scalarField2D::render
#
# USAGE
# scripting::scalarField2D::render
#
# PURPOSE
#
# This proc displays the contours or colorplane or both (depending on
# the configuration).
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::scalarField2D::render
#
#****
# ------------------------------------------------------------------------

proc scripting::scalarField2D::render {} {
    variable scalarField2D

    if { !$scalarField2D(loaded) } {
	ErrorDialog "can't render contours, the scalar field was not loaded. Call scripting::scalarField2D::load and scripting::scalarField2D::configure before calling scripting::scalarField2D::render"
	return
    }

    if { !$scalarField2D(configured) } {
	ErrorDialog "can't render contours, since they were not configured. Call scripting::scalarField2D::configure before calling scripting::scalarField2D::render"
	return
    }

    #############
    UpdatePropertyPlane
}
