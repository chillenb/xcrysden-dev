# ------------------------------------------------------------------------
#****** ScriptingExamples/pw_filters.tcl ***
#
# NAME
# pw_filters.tcl -- check PWscf I/O filters
#
# USAGE
# xcrysden --script pw_filters.tcl
#
# COPYRIGHT
# Anton Kokalj (C) 2004
#
# PURPOSE

# This is a multiScript that tests the PWscf I/O filters. It uses a multiScript
# utility (see example: multiScript.tcl for more detail).

#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Sometime in February 2004
# 
# SOURCE


scripting::multiScript {

    # ------------------------------------------------------------------------
    # This is the MAJOR script 
    #
    # It opens files, one at a time
    # ------------------------------------------------------------------------

    global env

    chdir [file join $env(XCRYSDEN_TOPDIR) examples PWSCF_files]

    # PWscf input version < 1.2
    scripting::exec scripting::filter::pwscfInput EthAl001-2x2.inp 2 {
	1 13   2 6   3 1 
    }
    # PWscf input version >= 1.2
    scripting::exec scripting::filter::pwscfInput CH3Rh111.inp 2
    
    
    # PWscf output version < 1.2
    scripting::exec scripting::filter::pwscfOutput -lc EthAl001-2x2.out 2 {
    	1 13   2 6   3 1 
    }
    # PWscf output version >= 1.2
    scripting::exec scripting::filter::pwscfOutput -lc CH3Rh111.out 2

} {

    # ------------------------------------------------------------------------
    # This is the MINOR script
    #
    # It setup a given structure
    # ------------------------------------------------------------------------


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
    update

    # be alive for 5 seconds
    DisplayUpdateWidget dialog "Example will be alive for 5 seconds."
    after 5000 {exit 0}
}

#******