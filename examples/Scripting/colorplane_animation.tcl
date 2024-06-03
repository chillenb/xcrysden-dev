# ------------------------------------------------------------------------
#****** ScriptingExamples/colorplane_animation.tcl ***
#
# NAME
# colorplane_animation.tcl -- a simple example for animating a contour plots
#
# USAGE
# xcrysden --script colorplane_animation.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to produce an animation
# of colorplane/contour plots by loading the XSF file which contains
# the description of the 3D scalar field. In this example the contours
# of the LUMO molecular-orbital are animated trough the CO molecule
# for all three three basal planes (i.e. XY, XZ, and YZ).

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

scripting::open --xsf $env(XCRYSDEN_TOPDIR)/examples/XSF_Files/CO_lumo.xsf.gz



# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
# this is Lighting-On modes:
#-----
scripting::lighting On

# choose a given 3D-display mode

#scripting::displayMode3D Stick
scripting::displayMode3D Pipe&Ball
#scripting::displayMode3D BallStick
#scripting::displayMode3D SpaceFill
#----


# ------------------------------------------------------------------------
# maximize the XCRYSDEN main window
# ------------------------------------------------------------------------

#scripting::mainWindow maximize
scripting::displayWindow fullscreen

# ------------------------------------------------------------------------
# load the 3D scalar field
# ------------------------------------------------------------------------

scripting::scalarField3D::load


# ------------------------------------------------------------------------
# hide the isosourface control window
# ------------------------------------------------------------------------

IsoControl_Hide
# alternatively you could do
#wm withdraw .iso


# ------------------------------------------------------------------------
# this is preliminary, used for the sake of better visualizing
# effect. See below the desciption of the two functions used here.
# ------------------------------------------------------------------------

scripting::scalarField3D::configure \
    -isosurface           0 \
    -interpolation_degree 2 \
    -basalplane           2 \
    -colorbasis           BLUE-WHITE-RED \
    -scalefunction        LINEAR \
    -colorplane           1 \
    -isoline              1 \
    -colorplane_lighting  0 \
    -cpl_transparency     1 \
    -cpl_thermometer      1 \
    -2Dlowvalue           -0.2 \
    -2Dhighvalue          +0.2 \
    -2Dnisoline           11 \
    -anim_step            1 \
    -current_slide        1 \
    -isoline_color        {property color} \
    -isoline_width        3 \

#    -isoline_monocolor    \#000000            
scripting::scalarField3D::render



# ------------------------------------------------------------------------
# zoom and rotate the structure 
# ------------------------------------------------------------------------

scripting::zoom -0.08 8 
scripting::rotate x -5 18
scripting::rotate y +5 15


# do the animation show forever

while {1} {

#
# description of scripting::scalarField3D::configure options:
#
# ------------------------------------------------------------------------
#   OPTION::                 ALLOWED-VALUES + Description          STATUS
# ------------------------------------------------------------------------
#   -isosurface              0|1                                 REQUIRED
#                            0 = do not render isosurface
#                            1 = render isosurface
#
#   -interpolation_degree    integer                             OPTIONAL
#                            degree of scalar-field interpolation
#
#   -isolevel                real                                REQUIRED
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
#   -basalplane              0|1|2
#                            show the ith basal plane (0=xy, 1=xz, 2=yz)
#                            as colorplane and/or isolines (as specified by
#                            -colorplane and -isoline options)
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
#
# Example:
#
# scripting::scalarField3D::configure \
#     -isosurface           1 \
#     -interpolation_degree 2 \
#     -isolevel             0.1 \
#     -plusminus            1 \
#     -revertsides          {pos neg} \
#     -revertnormal         {pos neg} \
#     -basalplane           0 \
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


    # ------------------------------------------------------------------------
    # configure, i.e., specify how to render the scalar field
    #
    # render 2nd basal plane
    # ------------------------------------------------------------------------
    scripting::scalarField3D::configure \
	-isosurface           0 \
	-interpolation_degree 2 \
	-basalplane           2 \
	-colorbasis           BLUE-WHITE-RED \
	-scalefunction        LINEAR \
	-colorplane           1 \
	-isoline              1 \
	-colorplane_lighting  0 \
	-cpl_transparency     1 \
	-cpl_thermometer      1 \
	-2Dlowvalue           -0.2 \
	-2Dhighvalue          +0.2 \
	-2Dnisoline           11 \
	-anim_step            1 \
	-current_slide        1 \
	-isoline_color        {property color} \
	-isoline_width        3 \

#	-isoline_monocolor    \#000000            
    
    # ------------------------------------------------------------------------
    # render the 3D scalar field as requested by 
    # scripting::scalarField3D::configure
    # ------------------------------------------------------------------------
    
    scripting::scalarField3D::render
    

    
    # ------------------------------------------------------------------------
    # now let us do the animation
    # ------------------------------------------------------------------------
    
    for {set i 1} {$i <= 49} {incr i} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    for {set i 49} {$i > 0} {incr i -1} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    scripting::scalarField3D::configure -colorplane 0 -isoline 0 -cpl_thermometer 0
    
    
    # now render the 1st basal plane

    scripting::scalarField3D::configure \
	-basalplane           1 \
	-colorbasis           BLUE-WHITE-RED \
	-scalefunction        LINEAR \
	-colorplane           1 \
	-isoline              1 \
	-colorplane_lighting  0 \
	-cpl_transparency     1 \
	-cpl_thermometer      1 \
	-2Dlowvalue           -0.2 \
	-2Dhighvalue          +0.2 \
	-2Dnisoline           11 \
	-anim_step            49 \
	-current_slide        1 \
	-isoline_color        {property color} \
	-isoline_width        3 \

#	-isoline_monocolor    \#000000
    scripting::scalarField3D::render
    scripting::rotate y 5 18    

    for {set i 1} {$i <= 49} {incr i} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    for {set i 49} {$i > 0} {incr i -1} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    scripting::scalarField3D::configure -colorplane 0 -isoline 0 \
	-cpl_thermometer 0    


    # now render the 3rd basal plane
    scripting::scalarField3D::configure \
	-basalplane           3 \
	-colorbasis           BLUE-WHITE-RED \
	-scalefunction        LINEAR \
	-colorplane           1 \
	-isoline              1 \
	-colorplane_lighting  0 \
	-cpl_transparency     1 \
	-cpl_thermometer      1 \
	-2Dlowvalue           -0.2 \
	-2Dhighvalue          +0.2 \
	-2Dnisoline           11 \
	-anim_step            1 \
	-current_slide        1 \
	-isoline_color        {property color} \
	-isoline_width        3 \

#	-isoline_monocolor    \#000000
    scripting::scalarField3D::render
    scripting::rotate y -5 18    
    scripting::rotate x +5 26 

    for {set i 1} {$i <= 49} {incr i} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    for {set i 49} {$i > 0} {incr i -1} {
	scripting::scalarField3D::configure -current_slide $i
	scripting::scalarField3D::render
    }
    scripting::rotate x -5 26
    scripting::scalarField3D::configure -colorplane 0 -isoline 0 \
	-cpl_thermometer 0      
}

#****
# ------------------------------------------------------------------------
