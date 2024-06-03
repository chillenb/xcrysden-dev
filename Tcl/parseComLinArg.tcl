#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/parseComLinArg.tcl                               #
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc parseComLinArg comopt {
    global system xcMisc

    set load_attributes 0
    set i 0
    foreach option $comopt {
	incr i

	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option

	    switch -exact -- $tag {
		"-?"   { 
		    XCrySDenUsage
		}
		"--help" {
		    XCrySDenUsage
		}
		"--quiet" {
		    # quiet run -- without welcome window !!!
		    ViewMol .
		    return
		}
	    }

	} else {    	    
	    switch -glob -- $tag {
		"-r" -
		"--reduce*" {
		    set xcMisc(reduce_to) $option
		    continue
		}
                "--print" {
                    # print the structure to a file and exit
                    if { [winfo exists .mesa] } {
                        set filedir [GetAbsoluteFileName $option]
			raise .
                        update
                        update idletask
			printTogl .mesa $filedir
                        update
                        exit 0
                    }
                }
	    }	
	    
	    set filedir [GetAbsoluteFileName $option]

	    if { ! [file exists $filedir] && $tag != "--wien_struct" } {
		destroy .title
		ErrorDialog "file \"$filedir\" does not exist !!!"
		exit
	    }

            switch -glob -- $tag {
		"-a" -
		"--attrib*" {
		    set load_attributes 1
		    set file_attributes $filedir
		}
		"-l" -
		"--light*" {
		    glLight:load $filedir
		}
		"-c" -
		"--custom" {
		    uplevel #0 "source $filedir"
		}

                "--xsf" -
		"--xcr" -
		"--animxsf" -
		"--axsf" {
                    # reload-capable
                    reloadRegister $filedir XSF
		    ViewMol .
		    xsfOpen $filedir .mesa
		}
		"--xsf2" {
                    # reload-capable
                    reloadRegister $filedir XSF2
		    ViewMol .
		    xsf2Open $filedir .mesa
		}
		"--xmol" -
		"--xyz" {
                    reloadRegister $filedir XYZ
		    ViewMol .
		    xyzOpen $filedir
		}
		"--pdb" {
                    reloadRegister $filedir PDB
		    ViewMol .
		    OpenXYZPDB .mesa pdb $filedir
		}
		"--fermi" -
		"--bxsf" {
                    set xcMisc(titlefile) $filedir
		    if { [file isfile $filedir] } {
			genFSInit $filedir
		    } else {
			ErrorDialog "file $filedir doesn't exists\n(or is not a regular file) !!!"
			exit
		    }
		}

		"--pw_inp" -
                "--pwi" {
                    # reload-capable
		    if { [file exists $filedir] } {
                        reloadRegister $filedir openExtStruct:reload
			ViewMol .
			openExtStruct 3 crystal external \
			    [list sh $system(TOPDIR)/scripts/pwi2xsf.sh] \
			    pwi2xsf.xsf_out \
			    {PWSCF Input File} ANGS \
			    -file $filedir \
			    -preset pwInputPreset
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}

		"--pw_out" -
                "--pwo" {
                    # reload-capable
		    if { [file exists $filedir] } {
                        reloadRegister $filedir openExtStruct:reload
			ViewMol .
			openExtStruct 3 crystal external \
			    [list sh $system(TOPDIR)/scripts/dummy.sh] \
			    pwo2xsf.xsf \
			    {PWSCF Output File} ANGS \
			    -file $filedir \
			    -preset pwOutputPreset
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}

		    
		"--fhi_inpini" {
		    if { [file exists $filedir] } {
                        reloadRegister $filedir openExtStruct:reload
			ViewMol .
			openExtStruct 3 crystal external34 \
                            $system(BINDIR)/fhi_inpini2ftn34 \
                            $system(ftn_name).34 {FHI98MD "ini.inp"} \
                            BOHR -file $filedir \
                            -preset [list fhiPreset inpini]
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}
		"--fhi_coord" {
		    if { [file exists $filedir] } {
                        reloadRegister $filedir openExtStruct:reload
			ViewMol .
			openExtStruct 3 crystal external \
                            $system(BINDIR)/fhi_coord2xcr \
                            fhi_coord.xcr {FHI98MD "coord.out"} \
                            BOHR -file $filedir \
                            -preset [list fhiPreset coord]
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}

		# --------------------------------------------------
		#  --crystal_* options are intended as CRYSTAL-95/98/03 interfacing
		# --------------------------------------------------

		"--crystal_inp" {
		    if { [file exists $filedir] } {
			ViewMol .
			OpenCyrstalInput $filedir
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}
		"--crystal_f9" {
		    if { [file exists $filedir] } {
			ViewMol .
			PropC95 $filedir
		    } else {
			ErrorDialog "$filedir doesn't exists !!!"
			exit
		    }
		}

		# --------------------------------------------------
		#  --wien_* options are intended as WIEN2k interfacing
		# --------------------------------------------------

		"--wien_struct" {
		    # option can be directory, filehead or file
		    if { [file isdirectory $filedir] } {
			set filehead [file tail $filedir]
			set file $filedir/$filehead.struct
		    } elseif { [file exists $filedir.struct] } {
			set file $filedir.struct
		    } elseif { [file exists $filedir] } {
			set file $filedir
		    } else {
			ErrorDialog "$filedir is neither directory, neither filehead, niether file"
			exit
		    }
                    reloadRegister $file WIEN2k
		    ViewMol .
		    wnOpenSFile $file
		}
		"--wien_kpath" {
		    # struct file must be read, and converted
		    # pop-up just BZ Window
		    if { [file isdirectory $filedir] } {
			wnDetComOpt $filedir
			wnKPath $filedir
		    } elseif { [file isfile $filedir] } {
			set dir [file dirname $filedir]
			wnDetComOpt $dir
			wnKPath $filedir			
		    } else {
			ErrorDialog "$filedir directory|file does not exist !!!"
			exit
		    }
		}
		"--wien_renderdensity" {
		    # use wnDensity proc and pop-up Viewer and 2D isoControl
		    # toplevel
		    if { [file isdirectory $filedir] } {
			set filehead [file tail $filedir]
			set file $filedir/$filehead.struct
			ViewMol .
			wnOpenSFile $file
			wnDensity $filedir
		    } else {
			ErrorDialog "$filedir directory does not exist !!!"
			exit
		    }
		}
		"--wien_density" {
		    # popup Viewer & another toplevel for inquiring 2D or 3D
		    # plots; then grid window is displayed
		    if { [file isdirectory $filedir] } {
			set filehead [file tail $filedir]
			set file $filedir/$filehead.struct
			ViewMol .
			wnOpenSFile $file
			wnDetComOpt $filedir
			# pop-up decision window
			wnDensity2D_or_3D $filedir
		    } else {
			ErrorDialog "$filedir directory does not exist !!!"
			exit
		    }		    
		}
		"--wien_f*" {
		    # struct file must be read, and converted
		    if [file isdirectory $filedir] {
			wnDetComOpt $filedir
			wnFSInit $filedir
		    } else {
			ErrorDialog "$filedir directory does not exist !!!"
			exit
		    }
		}
		"--gzmat" {
                    reloadRegister $filedir GZMAT
                    ViewMol .
		    gzmat $filedir
		}
		"--g98_out" -
		"--gXX_out" - 
		"--gxx_out" -
		"--gaussian_out" {
                    # reload-capable
                    reloadRegister $filedir addOption:reload
		    # this option is done via new addOption interface ...
		    addOption:hardcoded \
			[list sh $system(TOPDIR)/scripts/g98toxsf.sh] $filedir
		}
		"--cube" -		
		"--gxx_cube" -
		"--gXX_cube" -
		"--g98_cube" {
		    g98Cube $filedir
		}
		"--orca_out" {
                    # reload-capable
                    reloadRegister $filedir addOption:reload
		    # this option is done via new addOption interface ...
		    addOption:hardcoded \
			[list sh $system(TOPDIR)/scripts/orca2xsf.sh] $filedir
		}
		"-s" -
		"--script" {
		    scripting::source $filedir		    
		}
                "--os" -
                "--oldscript" {
                    scripting::oldscript $filedir	
                }
		default {
		    #
		    # check for user-custom options
		    reloadRegister $filedir addOption:reload
                    addOption:parse $tag $filedir
		}
	    }
	}
    }
	
    if { $i%2 } {
	XCrySDenUsage
    }
    unset i

    if { $load_attributes && [winfo exists .mesa] } {
	xsfLoadAttributes $file_attributes
    }

    # open the xcrysden main window if it is not already opened

    update
    if { ! [winfo exists .mesa] && ! [winfo exists .xcrysden_logo] } {
	ViewMol .
    }
}


proc UpdateWMTitle {file} {
    global xcMisc system
    
    set xcMisc(titlefile) $file
    wm title . "XCrySDen: [file tail $file]"

    # if file exists (does not always) copy it    
    if { [file exists $file] } {
	# the following copy is used if someone wants to save 
	# file as WIEN2k STRUCT FILE    	
	file copy -force $file $system(SCRDIR)/xc_struc.$system(PID)
    }
}


proc XCrySDenUsage {} {
    global system

    wm withdraw .
    set message [ReadFile $system(TOPDIR)/usage]
    puts stderr "\n$message\n"
    addOption:printCustomUsage
    flush stderr
    exit 1
}

