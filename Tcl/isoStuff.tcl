#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/isoStuff.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


proc SetXC_Iso {{dim 3D}} {
    global isosurf dif_isosurf prop isosign isodata 
	    
    if { ! [info exists dif_isosurf(dif_map)] } {
	set dif_isosurf(dif_map) 0
    }

    ######################
    # ISOSIGN parameters #
    ######################
    set isosign {}
    append isosign "xc_isosign 3 0 +1;   "
    append isosign "xc_isosign 2 0 +1;   "
    append isosign "xc_isosign 1 all +1; "
    if { $dif_isosurf(dif_map) == 1 } {
	###################
	# DIFFERENCE MAPS #
	###################
	append isosign "xc_isosign 3 1 -1; "
    }

    if { $prop(spin_case) } {
	#############
	# SPIN CASE #
	#############
	if { $isosurf(spin) == "ALPHA" } {
	    append isosign "xc_isosign 0 0 +1; "
	    append isosign "xc_isosign 0 1 +1;  "
	} elseif { $isosurf(spin) == "BETA" } {
	    append isosign "xc_isosign 0 0 +1;  "
	    append isosign "xc_isosign 0 1 -1; "
	} elseif { $isosurf(spin) == "ALPHA+BETA" } {
	    append isosign "xc_isosign 0 0 +1; "
	    append isosign "xc_isosign 0 1  0; "
	} elseif { $isosurf(spin) == "ALPHA-BETA" } {
	    append isosign "xc_isosign 0 0 0; "
	    append isosign "xc_isosign 0 1 +1; "
	}	
    } else {
	#################
	# NON-SPIN CASE #
	#################
	append isosign "xc_isosign 0 0 +1; "
    }

    ######################
    # ISODATA parameters #
    ######################
    set NPZ_1 [expr $prop(NPZ) - 1]
    set isodata {}
    if { $dif_isosurf(dif_map) == 0 } {
	if $prop(spin_case) {
	    if { $isosurf(spin) == "ALPHA+BETA" } {
		if { $dim == "3D" } {
		    append isodata "xc_isodata 0 0 0-$NPZ_1 0; "
		} elseif { $dim == "2D" } {
		    append isodata "xc_isodata 0 0 0 0; "
		}
	    } elseif { $isosurf(spin) == "ALPHA-BETA" } {
		if { $dim == "3D" } {
		    append isodata "xc_isodata 0 0 0-$NPZ_1 1; "
		} elseif { $dim == "2D" } {
		    append isodata "xc_isodata 0 0 0 1; "
		}
	    } elseif { $isosurf(spin) == "ALPHA" || \
		    $isosurf(spin) == "BETA" } {		    
		MakeIsoData $NPZ_1 0 0
	    }
	} else {
	    if { $dim == "3D" } {
		append isodata "xc_isodata 0 0 0-$NPZ_1 0; "
	    } elseif { $dim == "2D" } {
		append isodata "xc_isodata 0 0 0 0; "
	    }		
	}
    } elseif { $dif_isosurf(dif_map) == 1 } {
	###################
	# DIFFERENCE MAPS #
	###################
	if $prop(spin_case) {
	    if { $isosurf(spin) == "ALPHA+BETA" } {
		if { $dim == "3D" } {
		    append isodata "xc_isodata 0 0 0-$NPZ_1 0; "
		    append isodata "xc_isodata 1 0 0-$NPZ_1 0; "
		} elseif { $dim == "2D" } {
		    append isodata "xc_isodata 0 0 0 0; "
		    append isodata "xc_isodata 1 0 0 0; "
		}
	    } elseif { $isosurf(spin) == "ALPHA-BETA" } {
		if { $dim == "3D" } {
		    append isodata "xc_isodata 0 0 0-$NPZ_1 1; "
		    append isodata "xc_isodata 1 0 0-$NPZ_1 1; "
		} elseif { $dim == "2D" } {
		    append isodata "xc_isodata 0 0 0 1; "
		    append isodata "xc_isodata 1 0 0 1; "
		}
	    } elseif { $isosurf(spin) == "ALPHA" || \
		    $isosurf(spin) == "BETA" } {
		MakeIsoData $NPZ_1 0 0
		MakeIsoData $NPZ_1 1 0
	    }
	} else { 
	    # RHF
	    if { $dim == "3D" } {
		append isodata "xc_isodata 0 0 0-$NPZ_1 0; "
		append isodata "xc_isodata 1 0 0-$NPZ_1 0; "
	    } elseif { $dim == "2D" } {
		append isodata "xc_isodata 0 0 0 0; "
		append isodata "xc_isodata 1 0 0 0; "
	    }
	}
    }
}



proc MakeIsoData {NPZ_1 stack3 stack2} {
    global isodata
    # because xc_isodata must be specified in sequential fashion 
    # care must be taken

    xcDebug "MakeIsoData:: $NPZ_1"
    for {set i 0} {$i <= $NPZ_1} {incr i} {
	append isodata "xc_isodata $stack3 $stack2 $i 0; "
	append isodata "xc_isodata $stack3 $stack2 $i 1; "
    }
}



proc GetCrystalVec {which_vec nvec} {
    global vec system
    # which_vec.......either prim, conv; lately converted to PRIMVEC/CONVVEC
    # nvec............number of vectors to read!!!!

    if { $which_vec != "prim" && $which_vec != "conv" } {
	tk_dialog [WidgetName] ERROR "ERROR: this is bug in program, please report to author !!!\nCode: GetCrystalVec {which_vec}" error 0 OK
	return
    }
    
#    if { $which_vec == "prim" } { set which_vec "PRIMVEC" }
#    if { $which_vec == "conv" } { set which_vec "CONVVEC" }
    if { $which_vec == "prim" } { 
	set which_vec "PRIMVEC" 
    } else {
	set which_vec "CONVVEC" 
    }
	
    # xc_gengeom must automatically exists when In properties mode 
    # (c95_properties_render)
    set fileID [open "$system(SCRDIR)/xc_gengeom.$system(PID)" r]

    # initialise vec($,$)
    for {set i 0} {$i < 3} {incr i} {
	for {set j 0} {$j < 3} {incr j} {
	    set vec($i,$j) 0
	}
    }

    set n 0
    set output [split [read $fileID] \n]
    foreach line $output {
	xcDebug "line: $line"
	if [string match " $which_vec*" $line] {
	    set n1 [expr $n + 1]
	    set n2 [expr $n + 2]
	    set n3 [expr $n + 3]
	    for {set j 0} {$j < $nvec} {incr j} {
		set vec(0,$j) [Angs2Bohr [lindex [lindex $output $n1] $j]]
		set vec(1,$j) [Angs2Bohr [lindex [lindex $output $n2] $j]]
		set vec(2,$j) [Angs2Bohr [lindex [lindex $output $n3] $j]]
		xcDebug "PRIMVEC: $vec(0,$j) $vec(1,$j) $vec(2,$j)\n"
	    }
	    break
	}
	incr n
    }

    close $fileID
}



# this is new routine (Mon May  3 13:02:04 CEST 1999), which is now used 
# instead of GetCrystalVec, since xc_isospacesel is used now !!!
proc GetCageVecOrig pov {
    global vec isosurf
    # which_vec.......either prim, conv; lately converted to PRIMVEC/CONVVEC
    # nvec............number of vectors to read!!!!
    
    set isosurf(origin,0) [Angs2Bohr [lindex $pov 0]]
    set isosurf(origin,1) [Angs2Bohr [lindex $pov 1]]
    set isosurf(origin,2) [Angs2Bohr [lindex $pov 2]]
    for {set i 0} {$i < 3} {incr i} {
	for {set j 0} {$j <3} {incr j} {
	    set index [expr 3 + $i*3 + $j]
	    set vec($i,$j) [Angs2Bohr [lindex $pov $index]]
	}
    }
}



proc GetNPY {com} {
    global vec prop isosurf
    # vec.........crystal vectors in Bohrs
    
    if { $isosurf(res_type) == "points" } {
	set prop(NPY) $isosurf(resol_poi)
    } else {
	# ANGSTROMS OR BOHRS; convert to NPY
	if { $com == "ECHD" } {
	    # take vec(1); look at "GetNumberOfPoints"
	    set distY [expr sqrt( \
		    $vec(1,0) * $vec(1,0) + \
		    $vec(1,1) * $vec(1,1) + \
		    $vec(1,2) * $vec(1,2))]
	} else {
	    # for "ECHG and rest" take vec(0); look at "GetNumberOfPoints"
	    set distY [expr sqrt( \
		    $vec(0,0) * $vec(0,0) + \
		    $vec(0,1) * $vec(0,1) + \
		    $vec(0,2) * $vec(0,2))]
	}
	
	if { $isosurf(mb_angs/bohr) == "Angstroms" } {
	    set res [Angs2Bohr $isosurf(resol_ang)]
	} else {
	    set res $isosurf(resol_ang)
	}

	set prop(NPY) [expr round( $distY / $res )]
    }
}



proc GetNumberOfPoints {command dim} {
    global prop vec
    # dim.............dimension of map to render 2/3
    
    if { $command == "ECHD" } {
	#######################################################################
	# the order of points submited must to "properties" program must be:
	# ---------------------------------------------------------------------
	# for ECHD:
	#    (  1,NPY)
	#    (  1,  1)
	#    (NPX,  1)
	#
	# so the order of points is:
	#    ECHD_order == point2, point0, point1 
	# and because point2 is inheriter of vec(1), it is vec(1) who is 
        # criteria for NPY
	# ---------------------------------------------------------------------
	set disY [expr sqrt( \
		$vec(1,0) * $vec(1,0) + \
		$vec(1,1) * $vec(1,1) + \
		$vec(1,2) * $vec(1,2))]
	set dY [expr $disY / [expr $prop(NPY) - 1.0]]
	xcDebug "distY = $disY;  dY = $dY"
	if { $dim == 3} {
	    set disZ [expr sqrt(\
		    $vec(2,0) * $vec(2,0) + \
		    $vec(2,1) * $vec(2,1) + \
		    $vec(2,2) * $vec(2,2))]
	    set prop(NPZ) [expr round( $disZ / $dY ) + 1]
	    set dZ [expr $disZ / [expr $prop(NPZ) - 1.0]]
	    xcDebug "distZ = $disZ;  NPZ = $prop(NPZ)";
	    xcDebug "NPZ = $prop(NPZ);  dZ = $dZ"
	}
    } else { 
	#######################################################################
	# for ECHG:
	#    (  1,  1)
	#    (NPY,  1)
	#    (NPY,NPX)
	#
	# so the order of points is:
	#    ECHG_order == point0, point1, (point2 + vecA)
	# and because point1 is inheriter of vec(0), it is vec(0) who is
	# criteria for NPY
	# ---------------------------------------------------------------------
	set disY [expr sqrt( \
		$vec(0,0) * $vec(0,0) + \
		$vec(0,1) * $vec(0,1) + \
		$vec(0,2) * $vec(0,2))]
	set dY [expr $disY / [expr $prop(NPY) - 1.0]]
	if { $dim == 3 } {
	    set disZ [expr sqrt(\
		    $vec(2,0) * $vec(2,0) + \
		    $vec(2,1) * $vec(2,1) + \
		    $vec(2,2) * $vec(2,2))]
	    set prop(NPZ) [expr round( $disZ / $dY ) + 1]
	    set dZ [expr $disZ / [expr $prop(NPZ) - 1.0]]
	}
    }
}



proc IsoCalc {n command {c95_output {}} {dir {}}} {
    global system vec prop isosurf dif_isosurf 
    # n          ........ this is used for diff. maps to know is it map A or 
    #                     map B
    # c95_output ........ c95's output file
    # dir        ........ working directory for RunC95

    #############################
    # vec(-1,$) is DUMMY & must be (0.0, 0.0, 0.0)!!!!!!!!
    set vec(-1,0) 0.0
    set vec(-1,1) 0.0
    set vec(-1,2) 0.0
    for {set i 0} {$i < $prop(NPZ)} {incr i} {	
	if { $prop(NPZ) > 1 } {
	    set f [expr $i / [expr $prop(NPZ) - 1.0]]
	} else {
	    set f 1.0
	}
	for {set j -1} {$j < 2} {incr j} {	    
	    set jj [expr $j + 1]
	    set v($i,$jj,0) [expr $isosurf(origin,0) + \
		    $vec($j,0) + $f * $vec(2,0)]
	    set v($i,$jj,1) [expr $isosurf(origin,1) + \
		    $vec($j,1) + $f * $vec(2,1)]
	    set v($i,$jj,2) [expr $isosurf(origin,2) + \
		    $vec($j,2) + $f * $vec(2,2)]	    
	    xcDebug "$f\n$v($i,$jj,0)   $v($i,$jj,1)   $v($i,$jj,2)\n"
	}
    }

    ##########################################################################
    # the order of points submited to "properties" program must be:
    #
    # -------------------------------------------------------------------------
    # for ECHD:
    #    (  1,NPY)
    #    (  1,  1)
    #    (NPX,  1)
    #
    # so the order of points is:
    #    ECHD_order == point2, point0, point1 
    # -------------------------------------------------------------------------
    # for ECHG:
    #    (  1,  1)
    #    (NPY,  1)
    #    (NPY,NPX)
    #
    # so the order of points is:
    #    ECHG_order == point0, point1, (point2 + vecA)
    # -------------------------------------------------------------------------

    ###############################
    # take care of density matrix #
    ###############################
    if { $dif_isosurf(dif_map) == 1 } {
	if { $n == 1 } {
	    if { $dif_isosurf(denmat_A) == "Density matrix as superposition of atomic densities" } {
		if { $prop(PATO_newbasis_A) == 0 } {
		    set    input "PATO\n"
		    append input "0 0\n"
		} else {
		    ###############################
		    # new basis set was specified #
		    ###############################
		    # insert code here
		}		
	    }
	} elseif { $n == 2 } {
	    if { $dif_isosurf(denmat_B) == "SCF density matrix" && \
		    $dif_isosurf(denmat_A) != "SCF density matrix" } {
		set input "PSCF\n" 
	    }
	    if { $dif_isosurf(denmat_B) == "Density matrix as superposition of atomic densities" } {
		if { $prop(PATO_newbasis_B) == 0 } {
		    set    input "PATO\n"
		    append input "0 0\n"
		} else {
		    ###############################
		    # new basis set was specified #
		    ###############################
		    # insert code here
		}
	    }
	}
    }
    
    ##########################################################################
    # t.k. - this is temporal, just to se the gamma point 5-th band of Rh band
    ##########################################################################
    #append input "NEWK\n0 0 0\n0 0\nPBAN\n1 0 1\n1 0 0 0\n10\n"
    ##########################################################################

    append input "BOHR\n"
    ########################################
    # prototype model for PBAN
    # append input "NEWK\n0 0 0\n1 0\n"
    # append input "PBAN\n2 0\n1 2\n"
    ########################################
    if { $command == "POTM" } {
	append input "POLI\n4 1\n0\n"
    }
    for {set i 0} {$i < $prop(NPZ)} {incr i} {	    
	if { $command == "ECHD" } {
	    # ECHD
	    append input "ECHD\n$prop(NPY) 25 0 1\n"
	    append input "$v($i,2,0) $v($i,2,1) $v($i,2,2)\n"
	    append input "$v($i,0,0) $v($i,0,1) $v($i,0,2)\n"
	    append input "$v($i,1,0) $v($i,1,1) $v($i,1,2)\n"
	} else {
	    # ECHG & REST
	    if { [string match ECHG* $command] } {
		append input "$command\n$prop(NPY)\nCOORDINA\n"
	    }
	    if { $command == "POTM" } {
		append input "POTM\n0 5\n$prop(NPY)\nCOORDINA\n"
	    }
	    set p1 [expr $v($i,2,0) + $vec(0,0)]
	    set p2 [expr $v($i,2,1) + $vec(0,1)]
	    set p3 [expr $v($i,2,2) + $vec(0,2)]
	    append input "$v($i,0,0) $v($i,0,1) $v($i,0,2)\n"
	    append input "$v($i,1,0) $v($i,1,1) $v($i,1,2)\n"
	    append input "$p1 $p2 $p3\n"
	}
	if { [string match ECHG* $command] || $command == "POTM" } { 
	    append input "END\n" 
	}
    }
    append input "END\n"

    xcDebug -stderr "CRYSTALxx INPUT:"
    xcDebug -stderr "----------------"
    xcDebug -stderr $input

    #------------------------------------------------------------------------
    # t.k.: here user can make some manuall editing of CRYSTAL d3 input !!!
    #------------------------------------------------------------------------
    if { $prop(editScript) == 1 } {
	cd $system(SCRDIR)
	
	set file CRYSTAL-map$n.d3
	WriteFile $file $input w
	update
	xcEditFile $file -foreground
	
	# can-not chack if editing was OK -> to long to wait -> take a risk !!!
	
	#
	# drop the comments from file
	#
	set input {}
	foreach line [split [ReadFile -nonewline $file] \n] {
	    if { [string match *\#* $line] == 0 } {
		if { $line != {} } {
		    append input [format "%s\n" $line]
		}
	    }
	}
    }
    
    ###################################################################
    # because calculating "Grid of Points" can take some time, it's needed to 
    # give some feed back to the user -> message mode of RunC95    
    if ![RunC95 $system(c95_properties) {message  "The CRYSTAL program is calculating the grid of points now.\nIt can take some time, so PLEASE WAIT!!!"} \
	    $input $c95_output {} $dir] {
	# error occure
	# tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	return 0
    }
    return 1
}


