# ------------------------------------------------------------------------
#****** ScriptingExamples/atomic_labels.tcl ***
#
# NAME
# atomic_labels.tcl -- shows how to edit atomic labels
#
# USAGE
# xcrysden --script atomic_labels.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2003
#
# PURPOSE

# This is a scripting example that shows how to edit atomic labels
# using the scripting::atomicLabels namespace interface. There are two
# levels of atomic-labels, so-called, global and custom (see
# scripting::atomicLabels for explanation of the two items), and in
# this example we will edit both of them.

#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Fri Mar  7 17:12:59 CET 2003
# 
# SOURCE


# ------------------------------------------------------------------------
# load the structure (the argument to scription:open is the command line
# ------------------------------------------------------------------------

scripting::open --wien_struct $env(XCRYSDEN_TOPDIR)/examples/WIEN_struct_files/fe2p.struct


# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

#-----
# uncomment this is for Lighting-On modes:
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
# display atomic-labels, crystal-cells, unicolor-bonds, perspective
# ------------------------------------------------------------------------

foreach item { atomic-labels crystal-cells unicolor-bonds perspective } {
    scripting::display on $item
}


# ------------------------------------------------------------------------
# rotate and zoom the structure
# ------------------------------------------------------------------------

scripting::zoom     +0.60
scripting::rotate x +60
scripting::rotate y +20
scripting::rotate z +10


# ------------------------------------------------------------------------
# scripting::atomicLabels::global -- 
#
# let us change the global font. Create a new font by "font create ..."
# Tk mechanism and then load it by scripting::atomLabels::global function
# ------------------------------------------------------------------------

set global_font [font create -family Times -size 30 -weight bold]

scripting::atomicLabels::global -tkfont $global_font -brightcolor \#ffff55
wait 500

# ------------------------------------------------------------------------
# scripting::atomicLabels::atomID --
#
# lets make few CUSTOM-ATOMIC-LABELS. Custom atomic label are those one
# that were edited explicitly by scripting::atomicLabels::atomID or
# scripting::atomicLabels::atomType commands (unless cleared by 
# scripting::atomicLabels::clear).
# ------------------------------------------------------------------------

scripting::atomicLabels::atomID 21 -label "Fe(0th)" -brightcolor \#aaffaa
wait 500
scripting::atomicLabels::atomID 4  -label "Fe(1st)" -brightcolor \#aaffaa
wait 500
scripting::atomicLabels::atomID 5  -label "Fe(2nd)" -brightcolor \#aaffaa
wait 500
scripting::atomicLabels::atomID 17 -label "Fe(3rd)" -brightcolor \#aaffaa
wait 500
scripting::atomicLabels::atomID 11 -label "Fe(4th)" -brightcolor \#aaffaa
wait 500

# ------------------------------------------------------------------------
# scripting::atomicLabels::atomType --
#
# lets change the labels for all P-atoms
# ------------------------------------------------------------------------

set type_font [font create -family Helvetica -size 17 -weight bold]

scripting::atomicLabels::atomType P -label "Phosphorus" \
    -tkfont $type_font -brightcolor \#ffffff
wait 500


# ------------------------------------------------------------------------
# scripting::atomicLabels::clear --
#
# lets clear some custom labels, hence they will become global
# ------------------------------------------------------------------------

scripting::atomicLabels::clear atomID   21
wait 500
scripting::atomicLabels::clear atomType P
wait 500
scripting::atomicLabels::clear all
wait 500


# ------------------------------------------------------------------------
# lets make phosphorus labels again custom ones; change label to
# P(yellow)
# ------------------------------------------------------------------------

scripting::atomicLabels::atomType P -label "P(yellow)" \
    -tkfont $type_font -brightcolor \#ffffff


# ------------------------------------------------------------------------
# scripting::atomicLabels::render --
#
# lets toggle the display of some labels
# ------------------------------------------------------------------------

repeat 10 {
    scripting::atomicLabels::render global off
    wait 500
    scripting::atomicLabels::render global on
    wait 500
    scripting::atomicLabels::render custom off
    wait 500
    scripting::atomicLabels::render custom on    
    wait 500
    scripting::atomicLabels::render atomType Fe off
    wait 500
    scripting::atomicLabels::render atomType Fe on
    wait 500
    scripting::atomicLabels::render atomID 4 off
    wait 500
    scripting::atomicLabels::render atomID 4 on
    wait 500
    scripting::atomicLabels::render atomID 5 off
    wait 500
    scripting::atomicLabels::render atomID 5 on
    wait 500
    scripting::atomicLabels::render atomID 17 off
    wait 500
    scripting::atomicLabels::render atomID 17 on
    wait 500
    scripting::atomicLabels::render atomID 11 off
    wait 500
    scripting::atomicLabels::render atomID 11 on
    wait 500    
}

# ------------------------------------------------------------------------
#****
# ------------------------------------------------------------------------
