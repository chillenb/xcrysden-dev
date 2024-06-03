# ------------------------------------------------------------------------
#****** ScriptingExamples/movie.tcl ***
#
# NAME
# movie.tcl -- a simple example for making a MPEG movie
#
# USAGE
# xcrysden --script movie.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to produce an animation
# of a molecular structure and concomitantly creating an MPEG movie.

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

scripting::open --xyz $env(XCRYSDEN_TOPDIR)/examples/XYZ/mol2.xyz


# ------------------------------------------------------------------------
# resize the window
# ------------------------------------------------------------------------

scripting::displayWindow resize 400 400


# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
# this is for Lighting-On modes:
#-----
scripting::lighting On

# choose a 3D-display mode

#scripting::displayMode3D Stick
#scripting::displayMode3D Pipe&Ball
scripting::displayMode3D BallStick
#scripting::displayMode3D SpaceFill
#----

# #---
# # this is for Lighting-Off modes:
# #---
# scripting::lighting Off
# 
# # choose a 2D-display mode
# 
# #scripting::displayMode2D WireFrame  
# #scripting::displayMode2D PointLine  
# #scripting::displayMode2D Pipe&Ball  
# #scripting::displayMode2D BallStick-1 
# scripting::displayMode2D BallStick-2  
# #scripting::displayMode2D SpaceFill  
# #---



# ------------------------------------------------------------------------
# perform initial preparation for movie
# ------------------------------------------------------------------------


# 
# ::scripting::makeMovie::init --
#
#
# Description of options:
#
# ------------------------------------------------------------------------
#  OPTION               ALLOWED-VALUES + DESCRIPTION
# ------------------------------------------------------------------------
#  -gif_transp          0|1 --> make oblique|transparent animated-GIF
#
#  -gif_minimize        0|1 --> don't-minimize|minimize animateg-GIF
#
#  -gif_global_colormap 0|1 --> don't-use|use global colormap for animated-GIF
#
#  -movieformat         gif|mpeg --> create Animated-GIF|MPEG
#
#  -dir                 tmp|pwd --> put temporary (i.e. frame) files to 
#                                   scratch(tmp) or current working
#                                   directory(pwd) 
#
#  -frameformat         PPM|PNG|JPEG --> format of the frame-files 
#
#  -firstframe          positive-integer --> repeat first frame n-times
#
#  -lastframe           positive-integer --> repeat first frame n-times
#
#  -delay               positive-integer --> time dalay between frames 
#                                            in 1/100 sec
#  -save_to_file        file --> if specified the movie will be saved to file
#                                otherwise the filename will be queried
# 

scripting::makeMovie::init \
    -movieformat mpeg \
    -dir         tmp \
    -frameformat PPM \
    -firstframe  10 \
    -lastframe   10 \
    -delay       0

scripting::makeMovie::begin


# ------------------------------------------------------------------------
# first zoom the structure slowly
# ------------------------------------------------------------------------

repeat 10 {
    scripting::makeMovie::makeFrame
    scripting::zoom +0.03
}
    


# ------------------------------------------------------------------------
# rotate a few times
#
# usage: 
#    scripting::rotate x|y|z step_size number_of_times
# or
#    scripting::rotate xy|xz|yz step_size1 step_size2 number_of_times
#
# ------------------------------------------------------------------------

repeat 20 {
    scripting::makeMovie::makeFrame    
    scripting::rotate x +3
}

repeat 20 {
    scripting::makeMovie::makeFrame    
    scripting::rotate y +3
}

repeat 20 {
    scripting::makeMovie::makeFrame
    scripting::rotate z -3     
}

# ------------------------------------------------------------------------
# now finish the movie and save
# ------------------------------------------------------------------------

scripting::makeMovie::end

#****
# ------------------------------------------------------------------------
