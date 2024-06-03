#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/check_package.tcl
# ------                                                                    #
# Copyright (c) 1996--2014 by Anton Kokalj                                  #
#############################################################################

proc check_package_awk {} {
    global system

    if { ![info exists system(awk)] } {
	# check first for GNU awk
	set system(awk) [auto_execok gawk]
	if { $system(awk) == "" } {
	    set system(awk) [auto_execok awk]
	    if { $system(awk) == "" } {
		WarningDialog "couldn't find \"awk\" program" "some features will not work !!!"
	    }
	}
    }
}


proc check_package_terminal {} {
    global system env

    if { ! [info exists system(term)] } {
	# check for env(TERM)	
	if { [info exists env(TERM)] } {
	    if { $env(TERM) != "dump" } {
		set system(term) [auto_execok $env(TERM)]
		if { $system(term) != ""  } {
		    return
		}
	    }
	}

	# now we check for these: 
	foreach term {xterm xvt rxvt eterm gnome-teminal konsole roxterm} {
	    set system(term) [auto_execok $term]
	    if { $system(term) != "" } {
		break
	    }
	}
    }
}


proc check_package_crystal {} {
    global system

    if { $system(c95_exist) } {
	if { [info exists system(c95_crystal)] } {
	    set system(c95_integrals) $system(c95_crystal)
	}

	#
	# CRYSTAL-95/98/03/06/09/14/...
	#
	set crystal_module $system(c95_integrals)
	set input "test\nMOLECULE\n1\n1\n1 0.0 0.0 0.0\nSTOP\nEND"
	WriteFile xc_inp.$system(PID) $input w
	
	set system(c95_integrals) [auto_execok $system(c95_integrals)]

	set OK 1
	if { ! [file exists $system(c95_integrals)] } {
	    set OK 0
	    ErrorDialog "File \"$crystal_module\" does not exist";
	}
	if { $OK && ! [file executable $system(c95_integrals)] } {
	    set OK 0
	    ErrorDialog "File \"$crystal_module\" is not executable";
	}
	
	if { $OK } {
	    set _status [catch {
		exec $system(c95_integrals) < xc_inp.$system(PID) >& xc_tmp.$system(PID)
	    } _errMsg]
	    
	    if { $_status } {
		ErrorDialog \
		    "Couldn't run CRYSTAL package.\n\nPlease check its definition in ~/.xcrysden/custom-definition file.\n\nError Message: $_errMsg"
		set system(c95_version) none
	    }
	    
	    if { [file exists xc_tmp.$system(PID)] } {
		set file [ReadFile xc_tmp.$system(PID)]
		foreach line [split $file \n] {
		    switch -glob -- $line {
			"*C R Y S T A L*" {
			    set vf [expr [llength $line] - 2]
			    set system(c95_version) [lindex $line $vf]
			    break
			}
			"*    CRYSTAL03    *" {
			    set system(c95_version) "03"
			    break
			}
			"*    CRYSTAL06    *" {
			    set system(c95_version) "06"
			    break
			}
			"*    CRYSTAL09    *" {
			    set system(c95_version) "09"
			    break
			}
			"*    CRYSTAL14    *" {
			    set system(c95_version) "14"
			    break
			}
			"*    CRYSTAL17    *" {
                            puts stderr "Package CRYSTAL: version-17 found"
                            puts stderr "        (BEWARE: adopting version-14 instead,"
                            puts stderr "                 support for version-17 may not be fully functional)"
			    set system(c95_version) "14"
			    break
			}
		    }	    
		}
	    }
	    if { ! [info exists system(c95_version)] } {
		puts stderr "Package CRYSTAL: unsupported version"
		puts stderr "                 adopting a crystal version 14 (BEWARE: this may not work properly)"
		set system(c95_version) "14"
	    }

	    puts stderr "Package CRYSTAL: $system(c95_integrals) (version: $system(c95_version))"

	    if { ! [info exists system(c95_scf)] } {
		set system(c95_scf) $system(c95_integrals)
	    }
	} else {
	    # !!! OK == 0 !!!
	    set system(c95_exist)   0
	    set system(c95_version) none
	}
    } else {
	set system(c95_version) none
    }
}
