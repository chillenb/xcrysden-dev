#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/kPath.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc kPath {{saveonly 0}} {
    global kpath prop BzOK system
   
    if { $saveonly == 0 } {
	set BzOK(done) 0
	if { ! [Bz_MakeToplevel generic] } {
	    return
	}
    }

    # so OK button was pressed
    # save prop(c95_BAND_script)
    set filetypes {
	{{XCRYSDEN k-Path File       }   {.kpf}   }
	{{PWscf K_POINTS File        }   {.pwscf} }
	{{CRYSTAL Property Input File}   {.d3 }   }
	{{All Files                  }   {*   }   }
    } 

    set sfile [tk_getSaveFile -initialdir $system(PWD) \
		   -title "Save k-Path" \
		   -defaultextension .kpf \
		   -filetypes $filetypes]

    # maybe Cancel button was pressed
    if { $sfile == {} } { 
	set kpath(active) 0
	return 
    }

    #
    # write the klist file in appropriate format 
    #
    set ext [file extension $sfile]
    if { $ext == ".kpf" } {
	#
	# XCRYSDEN's kPath file 
	#
	set content [kPath_XCrySDen]
	WriteFile $sfile $content w
    } elseif { [string match -nocase "*pwscf*" $ext] } {
	#
	# PWscf K_POINTS card
	#
	kPath_PWscf $sfile
    } elseif { $ext == ".d3" } {
	#
	# write CRYSTAL d3 file
	#
	set content $prop(c95_BAND_script)
	WriteFile $sfile $content w
    } else {
	#
	# supported format
	#
	WarningDialog "unknown k-list file type \"$ext\". Choose among supported file-types"
	kPath 1
	return
    }
    set kpath(active) 0
}


#
# Save selected k-path file in kind of xcrysden k-path info file
#
proc kPath_XCrySDen {} {
    global kpath BzOK

    set content {
#------------------------------------------------------------------------
# This is an XCRYSDEN k-Path File. Here you will find the
# data about the selected special k-point coordinates !!!
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# WARNING: 
#   special k-point coordinates in this file are in "CRYSTAL" 
#   (i.e. fractional) coordinates, which are expressed either 
#   in the basis of "primitive" or "conventional" reciprocal 
#   lattice vectors. See below !!!
#------------------------------------------------------------------------


k-point coordinates are expressed in the basis: }
    
    append content [format "%s\n" $kpath(basis)]	
    append content {
#
# #-------------------------------------#
# # INTEGER FORM of k-point COORDINATES #
# #-------------------------------------#
#
# some programs want an integer form of k-point coordinates,
# this looks like:
#
# k-point == (M*kx, M*ky, M*kz), where M is an integer multiplier
#
# Actual k-point coordinates are then calculated as:
#
# k-point == 1/M * (M*kx, M*ky, M*kz)
#------------------------------------------------------------------------

Multiplier M: }

    append content [format "%s\n" $BzOK(iss)]    
    append content "\nInteger form of k-point coordinates (M*kx,M*ky,M*kz,label):\n"
    append content [format "%s\n" $kpath(point_labels)]	
    append content {
#
# #----------------------------------#
# # REAL FORM of k-point COORDINATES #
# #----------------------------------#
#
    }

    append content "\nReal form of k-point coordinates (kx,ky,kz,label):\n"
    foreach line [split $kpath(point_labels) \n] {
	if { [llength $line] >= 3 } {
	    for {set i 0} {$i < 3} {incr i} {
		set k($i) \
		    [expr double([lindex $line $i]) / double($BzOK(iss))]
	    }
	    set label [lindex $line 3]
	    if { $label == {} } {
		set label A
	    }
	    append content [format "%15.10f  %15.10f  %15.10f     %s\n" \
				$k(0) $k(1) $k(2) $label]
	}
    }
    append content {
#
# #-------------#
# # END of FILE #
# #-------------#
#
    }
}



proc kPath_PWscf {sfile} {
    global system wnKP BzOK Bz

    set k_list_file $system(SCRDIR)/xc_pwklist.$system(PID)

    set    k_list "$wnKP(npoi) $BzOK(nK)\n"
    append k_list "$wnKP(type)\n"
    for {set i 1} {$i <= $wnKP(npoi)} {incr i} {
	append k_list "$wnKP(poi,$i,1) $wnKP(poi,$i,2) $wnKP(poi,$i,3)   $wnKP(label,$i)\n"
    }
    WriteFile $k_list_file $k_list w

    xcDebug -stderr "XCRYSDEN-to-PWscf k-list file:"
    xcDebug -stderr "-----------------------------"
    xcDebug -stderr $k_list

    xcCatchExecReturn \
	$system(FORDIR)/pwKPath \
	$system(SCRDIR)/xc_struc.$system(PID) $k_list_file > $sfile
    
    #
    # print some supporting information
    #
        if { $wnKP(type) == "prim" } {
	set type "RECIPROCAL-PRIMITIVE"
    } else {
	set type "RECIPROCAL-CONVENTIONAL"
    }
	
    set    supportInfo "================================================================\n"
    append supportInfo "    * * * S U P P O R T I N G    I N F O R M A T I O N * * *    \n"
    append supportInfo "================================================================\n"
    append supportInfo "\n"
    #append supportInfo "BEWARE: the conversion of selected k-point coordinates to format\n"
    #append supportInfo "        suitable for the PWscf is a new and untested feature.\n"
    #append supportInfo "\n"
    #append supportInfo "        !!! USE AT YOUR OWN RISK !!!\n"
    #append supportInfo "        !!! PLEASE CHECK BELOW DATA FOR CONSISTENCY !!!\n"
    #append supportInfo "\n"
    append supportInfo "Number of selected k-points: $wnKP(npoi)\n\n"
    append supportInfo "**ORIGINAL DATASET**\n"
    append supportInfo "--------------------\n"
    append supportInfo "THE FOLLOWING k-POINT WERE SELECTED:\n"
    append supportInfo "(crystal coordinates with respect to $type vectors)\n\n"
    for {set i 1} {$i <= $wnKP(npoi)} {incr i} {
	append supportInfo [format "   % 10.5f   % 10.5f   %10.5f     %s\n" \
				$wnKP(poi,$i,1) $wnKP(poi,$i,2) $wnKP(poi,$i,3)    $wnKP(label,$i)]
    }
    append supportInfo "\n**REPRODUCED(TRANSFORMED) DATASET for PWscf**\n"
    append supportInfo "----------------------------------------------\n"
    append supportInfo "IMPORTANT:\n"
    append supportInfo "Check if the data presented below are consistent with the data from ORIGINAL DATASET !!!\n\n"
    append supportInfo [ReadFile supportInfo-1.kpath]
    append supportInfo "\n"
    #append supportInfo "Transformed selected k-points for PWscf:\n"
    #append supportInfo "----------------------------------------\n"
    #append supportInfo [ReadFile supportInfo-2.kpath]
    #append supportInfo "\n"
    append supportInfo "The content of the generated PWscf K_POINTS file is the following:\n"
    append supportInfo "----------------------------------------------------------------\n"
    append supportInfo [ReadFile $sfile]

    WriteFile $system(PWD)/supportInfo.kpath $supportInfo w
    if { [file exists supportInfo-1.kpath] } { file delete supportInfo-1.kpath }
    if { [file exists supportInfo-2.kpath] } { file delete supportInfo-2.kpath }
    
    set t [xcDisplayVarText $supportInfo "k-path: Supporting Information"]
    if { [winfo exists $t.f1.t] } {
	catch {$t.f1.t configure -height 30 -width 90}
	catch {$t.f1.t configure -state disabled}
    }
    tkwait window $t
}
