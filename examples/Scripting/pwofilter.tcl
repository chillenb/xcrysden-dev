# ------------------------------------------------------------------------
#****** ScriptingExamples/pwfilter.tcl ***
#
# NAME
# pwfilter.tcl -- check of PWscf Output filter
#
# USAGE
# xcrysden --script pwofilter.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2004
#
# CREATION DATE
# Sometime in February 2004
# 
# SOURCE

global system

set dir [file join $env(XCRYSDEN_TOPDIR) examples PWSCF_files]
cd $dir
#scripting::filter::pwscfOutput -lc EthAl001-2x2.out 2 {
#    1 13   2 6   3 1 
#}
scripting::filter::pwscfOutput -lc CH3Rh111.out 3

# ------------------------------------------------------------------------
# display the structure in appropriate display-mode
# ------------------------------------------------------------------------

scripting::lighting On
scripting::displayMode3D spacefill

# ------------------------------------------------------------------------
# zoom and rotate the structure 
# ------------------------------------------------------------------------

scripting::zoom   0.4
scripting::rotate x +60 
scripting::rotate y +20 


# ------------------------------------------------------------------------
# display 3x3x1 cell
# ------------------------------------------------------------------------
scripting::buildCrystal 3 3 1

#******