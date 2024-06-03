# ------------------------------------------------------------------------
#****** ScriptingExamples/build_crystal.tcl ***
#
# NAME
# build_crystal.tcl -- animates and builds a crystal structure
#
# USAGE
# xcrysden --script build_crystal.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a simple script. It renders a crystal structure, and builds
# step by step bigger portions of the crystal.

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

scripting::open --xsf $env(XCRYSDEN_TOPDIR)/examples/XSF_Files/ZnS.xsf


# ------------------------------------------------------------------------
# resize the display window to 400x400
# ------------------------------------------------------------------------

scripting::displayWindow resize 400 400


# ------------------------------------------------------------------------
# let us make a white background
# ------------------------------------------------------------------------

set myParam(BACKGROUND)      {1.00 1.00 1.00}
scripting::load_myParam


# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
#  this is for Lighting-On modes:
#-----
scripting::lighting On

# choose a 3D-display mode

#scripting::displayMode3D Stick
#scripting::displayMode3D Pipe&Ball
#scripting::displayMode3D BallStick
scripting::displayMode3D SpaceFill
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
# show th crystal-cells and switch to perspective projection
# ------------------------------------------------------------------------

scripting::display on crystal-cells
scripting::display on perspective


# ------------------------------------------------------------------------
# first zoom the structure slowly
# ------------------------------------------------------------------------

scripting::zoom +0.03 15


# ------------------------------------------------------------------------
# rotate a few times
#
# usage: 
#    scripting::rotate x|y|z step_size number_of_times
# or
#    scripting::rotate xy|xz|yz step_size1 step_size2 number_of_times
#
# ------------------------------------------------------------------------

scripting::rotate x +3 10
scripting::rotate y +3 10
scripting::rotate z -3 10    


# ------------------------------------------------------------------------
# now build the crystal slowly
# ------------------------------------------------------------------------

for {set i 1} {$i < 3} {incr i} {
    for {set j $i} {$j < 3} {incr j} {
	for {set k $j} {$k < 3} {incr k} {
	    if { $i*$j*$k != 8 } {
		scripting::buildCrystal $i $j $k
		wait 100
	    }
	}
    }
}

for {set i 2} {$i < 4} {incr i} {
    for {set j $i} {$j < 4} {incr j} {
	for {set k $j} {$k < 4} {incr k} {
	    if { $i*$j*$k != 27 } {
		scripting::buildCrystal $i $j $k
		update
	    }
	}
    }
}

for {set i 3} {$i < 5} {incr i} {
    for {set j $i} {$j < 5} {incr j} {
	for {set k $j} {$k < 5} {incr k} {
	    scripting::buildCrystal $i $j $k
	    update
	}
    }
}
#****
# ------------------------------------------------------------------------
