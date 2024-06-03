# ------------------------------------------------------------------------
#****** ScriptingExamples/animation.tcl ***
#
# NAME
# animation.tcl -- a simple example for animating a molecular structure
#
# USAGE
# xcrysden --script animation.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to produce an animation
# of a molecular structure. In this example the structure is rotated
# infinitely in different directions (see the scripting:rotate below).

#
# WARNINGS

# In this example XCRYSDEN switches to fullscreen mode. To exit from
# fullscreen mode double-click the first mouse button. Note also that
# in fullscreen mode the right-mouse button triggers the pop-up menu

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

scripting::open --xyz $env(XCRYSDEN_TOPDIR)/examples/XYZ/mol1.xyz


# ------------------------------------------------------------------------
# switch to fullscreen mode
# ------------------------------------------------------------------------

scripting::displayWindow fullscreen


# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
# uncomment this is for Lighting-On modes:
#-----
# scripting::lighting On
# 
# # choose a 3D-display mode
# 
# #scripting::displayMode3D Stick
# #scripting::displayMode3D Pipe&Ball
# scripting::displayMode3D BallStick
# #scripting::displayMode3D SpaceFill
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
# first zoom the structure slowly
# ------------------------------------------------------------------------

scripting::zoom +0.02 5


# ------------------------------------------------------------------------
# rotate a few times
#
# usage: 
#    scripting::rotate x|y|z step_size number_of_times
# or
#    scripting::rotate xy|xz|yz step_size1 step_size2 number_of_times
#
# ------------------------------------------------------------------------

scripting::rotate x +3 20
scripting::rotate y +3 20
scripting::rotate z -3 10    

# ------------------------------------------------------------------------
# rotate infinitely
# ------------------------------------------------------------------------

while {1} {
    scripting::rotate xy +5 +1 3
    scripting::rotate xy +5 +3 5
    scripting::rotate xy +5 +5 20
    scripting::rotate xy +5 +3 5
    scripting::rotate xy +5 +1 3
    scripting::rotate xy +5 +0 1
    scripting::rotate xz +5 +1 3
    scripting::rotate xz +5 +3 5
    scripting::rotate xz +5 +5 20
    scripting::rotate xz +5 +3 5
    scripting::rotate xz +5 +1 3
    scripting::rotate xz +5 +0 1
}

#****
# ------------------------------------------------------------------------
