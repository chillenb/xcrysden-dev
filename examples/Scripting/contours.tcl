# ------------------------------------------------------------------------
#****** ScriptingExamples/contours.tcl ***
#
# NAME
# contours.tcl -- plots contours+colorplane from XSF
#
# USAGE
# xcrysden --script contours.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to produce colorplane and
# isoline plot by loading the XSF file which contains the description
# of the 2D scalar field.

#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Sometime in February 2003
# 
# SOURCE


# ------------------------------------------------------------------------
# load the structure (the argument to scription:open is the command line
# ------------------------------------------------------------------------

scripting::open --xsf $env(XCRYSDEN_TOPDIR)/examples/XSF_Files/mol-urea2D.xsf



# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
# uncomment this is for Lighting-On modes:
#-----
scripting::lighting On

# choose a 3D-display mode

#scripting::displayMode3D Stick
scripting::displayMode3D Pipe&Ball
#scripting::displayMode3D BallStick
#scripting::displayMode3D SpaceFill
#----



# ------------------------------------------------------------------------
# zoom and rotate the structure 
# ------------------------------------------------------------------------

scripting::zoom +0.1 
scripting::rotate x -90 
scripting::rotate y -90 


# ------------------------------------------------------------------------
# load the 3D scalar field
# ------------------------------------------------------------------------

scripting::scalarField2D::load


#
# description of scripting::scalarField2D::configure options:
#
# ------------------------------------------------------------------------
#   OPTION::                 ALLOWED-VALUES + Description          STATUS
# ------------------------------------------------------------------------
#   -interpolation_degree    integer                             OPTIONAL
#                            degree of scalar-field interpolation
#
#   -colorbasis              MONOCHROME|RAINBOW|RGB|GEOGRAPHIC|BLUE-WHITE-RED
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
#   -expand2D_X	             positive-integer
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
#
# Example:
#
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


# ------------------------------------------------------------------------
# configure, i.e., specify how to render the scalar field
# ------------------------------------------------------------------------
scripting::scalarField2D::configure \
    -interpolation_degree 3 \
    -colorbasis           BLUE-WHITE-RED \
    -scalefunction        LINEAR \
    -colorplane           1 \
    -isoline              1 \
    -colorplane_lighting  0 \
    -cpl_transparency     0 \
    -cpl_thermometer      1 \
    -2Dlowvalue           -0.1 \
    -2Dhighvalue          +0.1 \
    -2Dnisoline           11 \
    -isoline_color        monocolor \
    -isoline_width        3 \
    -isoline_monocolor    \#000000


# ------------------------------------------------------------------------
# hide the propertyplane control window
# ------------------------------------------------------------------------

IsoControl_Hide .iso2D
# alternatively you could do
#wm withdraw .iso2D


# ------------------------------------------------------------------------
# render the 2D scalar field as requested by 
# scripting::scalarField2D::configure
# ------------------------------------------------------------------------

scripting::scalarField2D::render


# ------------------------------------------------------------------------
# now lets print to file what we have on the display window
# ------------------------------------------------------------------------

# this will query the filename:
scripting::printToFile; # here the colorplane legend is NOT printed 

# # while this will print directly to print.png:
# scripting::printToFile print.png; # here the colorplane legend is printed 

# # while this will print directly to print.png:
# scripting::printToFile print.png windowdump; # here the colorplane legend is printed 

#****
# ------------------------------------------------------------------------
