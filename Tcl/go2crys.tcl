#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/go2crys.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# calculate & display STRUCTURE
proc CalStru {} {
    global nxdir nydir nzdir n_groupsel system periodic geng inp \
	    n_groupsel groupsel err species

    xcDebug -debug "IN CalStru"

    if { $species == "external" } {
	set groupsel   "P 1"
	set n_groupsel 1
	set inp(IFLAG) 0
    } else {
	if { ![info exists inp(IGR)] } { set inp(IGR) $n_groupsel }
	if { ![info exists inp(AGR)] } { set inp(AGR) $groupsel }
    }

    # update xcState
    xcUpdateState
    set geoInput [MakeInput]
    
    xcDebug -debug "CRYSTAL INPUT:\n $geoInput"
    xcDebug -debug "NXDIR NYDIR NZDIR: $nxdir $nydir $nzdir"
    
    if { ! [RunC95 $system(c95_integrals) {} $geoInput] } {
	return 0
    }

    # workaround of crystal-bug in EXTPRT (it prints number-of-atoms only up to 999 atoms, and then ***)    
    if { $species == "external" } {
	file copy -force $system(SCRDIR)/external_unit34   $system(SCRDIR)/$system(ftn_name).34
    }
    
    # usage of "gengeom" program:
    # 
    # gengeom  MODE1  MODE2  MODE3  IGRP  NXDIR  NYDIR  NZDIR  OUTPUT  INPUT
    #    0       1      2      3      4     5      6      7      8       9
    #
    # FIND FROM UNIT34 THE DIMENSIONALITY OF THE SYSTEM
    cd $system(SCRDIR)
    
    GenGeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) 1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID)

    #xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) 1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID)

    set fileID [open "$system(SCRDIR)/xc_gengeom.$system(PID)" r]
    GetDimGroup periodic(dim) periodic(igroup) $fileID
    close $fileID

    if { $periodic(dim) == 3 && $nzdir == 0 } { 
	# this can happen for "EXTERNAL"
	set nzdir 1 
    }

    # in periodic(dim) != 3 than updata (nxdir,nydir,nzdir)
    if { $periodic(dim) <= 2 } {
	set nzdir 0
    } 
    if { $periodic(dim) <= 1 } {
	set nydir 0
    }
    if { $periodic(dim) == 0 } {
	set nxdir 0
	# try this: due to problem reported by Klauss Doll
	GenGeom $geng(M1_PRIM) 1 $geng(M3_ARGUMENT) 1  1 1 1 $system(SCRDIR)/xc_struc.$system(PID)
    }
    
    #exec rm -f $system(SCRDIR)/xc_tmp.$system(PID)

    xcDebug -debug "GENGEOM>> dim igroup:: $periodic(dim) $periodic(igroup)"
    if { $periodic(dim) == 3 } {
	# take special care of hexagonal and trigonal systems
	if { $periodic(igroup) == 1 } {
	    # group can still be trigonal or hexagonal
	    if { $inp(IFLAG) == 0 } {
		if { $n_groupsel >= 143 && $n_groupsel < 168 } {
		    set periodic(igroup) $geng(IGRP_TRIG)
		} 
		if { $n_groupsel >= 168 && $n_groupsel < 195 } {
		    set periodic(igroup) $geng(IGRP_HEXA)
		} 
	    } else {
		set tt [lindex $inp(AGR) 1]
		if { [string match "3*" $tt] || [string match "-3*" $tt] } {
		    set periodic(igroup) $geng(IGRP_TRIG)
		}
		if { [string match "6*" $tt] || [string match "-6*" $tt] } {
		    set periodic(igroup) $geng(IGRP_HEXA)
		}
	    }
	}
    } 

    puts stdout "igroup:: $periodic(igroup)"
    flush stdout

    if { [xcIsActive render] } {
	# CellMode will do all the job
	GenGeomDisplay 1
    } else {
	GenGeomDisplay
	OpenXSFStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
    }

    CrysFrames

    return 1
}
    

proc MakeInput {} {
    global job_title n_groupsel groupsel species speciesName  \
	periodic inp AdvGeom H K L ISUP NL nxdir nydir nzdir mode \
	cxx

    # $save is used by SaveCrystalInput proc
            
    if { ! [info exists cxx(inputContent)] } {
	# JOB TITLE
	if { $job_title == "" } { set job_title "Input file constructed by XCrySDen" }
	set geoInput "$job_title\n"
	# SPECIES
	append geoInput [string toupper $species\n]

	if { $species != "external" } {

	    # GEOM-SYMM input for CRYSTAL
	    if { $species == "crystal" } {
		# RECORD: IFLAG IFHR IFSO
		append geoInput "$inp(IFLAG) $inp(IFHR) $inp(IFSO)\n"
		# RECORD: IGR
		if { $inp(IFLAG) == 0 } {
		    append geoInput "$n_groupsel\n"
		} else {
		    append geoInput "$groupsel\n"
		}
		puts stdout "#################GROUPSEL:: $groupsel"
		if { $inp(IFSO) > 1 } {
		    append geoInput "$inp(IX) $inp(IY) $inp(IZ)\n"
		}
	    }
	    
	    if { $species != "crystal" } {
		append geoInput "$n_groupsel\n"
	    }	
	    if { $species != "molecule" } {
		# RECORD: min. set of crystall. par.
		xcDebug "inp(A) = $inp(A)"
		SetUnDefPar
		append geoInput "$inp(A) $inp(B) $inp(C) \
		$inp(ALFA) $inp(BETA) $inp(GAMMA)\n"
	    }
	    
	    # NATR & Nst X Y Z
	    append geoInput "$inp(NATR)\n"
	    for {set i 1} {$i <= $inp(NATR)} {incr i} {
		append geoInput "$inp(NAT,$i) $inp(X,$i) $inp(Y,$i) $inp(Z,$i)\n"
	    }
	    
	    if [info exists speciesName] {
		set spec $speciesName
	    } else {
		set spec $species
	    }
	    # take care of (nxdir,nydir,nzdir)
	    if { $spec == "crystal" } {
		if { $nxdir == 0 } { set nxdir 1 }
		if { $nydir == 0 } { set nydir 1 }
		if { $nzdir == 0 } { set nzdir 1 }
		set periodic(dim) 3
	    }
	    
	    if { $spec == "slab" } {
		if { $nxdir == 0 } { set nxdir 1 }
		if { $nydir == 0 } { set nydir 1 }
		set periodic(dim) 2
		set nzdir 0
	    }
	    
	    if { $spec == "polymer" } {
		if { $nxdir == 0 } { set nxdir 1 }
		set periodic(dim) 1
		set nydir 0
		set nzdir 0
	    }
	    
	    if { $spec == "molecule" } {
		set periodic(dim) 0
		set nxdir 0
		set nydir 0
		set nzdir 0
	    }
	}
    } else {
	set geoInput $cxx(inputContent)
    }    

    ###############################
    # ADVANCE GEOMETRICAL OPTIONS #
    ###############################
    # check first if there is $AdvGeom(input) present
    if [info exists AdvGeom(input)] {
	append geoInput $AdvGeom(input)
    }
    # does user already make some AdvGeom
    set na [xcAdvGeomState current]
    xcDebug "\n\n ADV_GEOM ARRAY NAMES:: [array names AdvGeom *] $na\n\n"
    for {set i 1} {$i <= $na} {incr i} {
	if [info exists AdvGeom($i,rotate)] {
	    ##########
	    # ROTATE #
	    ##########
	    xcDebug -debug "ROTATE OPTION"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,rotate)"
	}    
	if [info exists AdvGeom($i,superCell)] {
	    #############
	    # SUPERCELL #
	    #############
	    xcDebug -debug "SUPERCELL OPTION"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,superCell)"
	}    
	if [info exists AdvGeom($i,elastic)] {
	    ###########
	    # ELASTIC #
	    ###########
	    xcDebug -debug "ELASTIC OPTION"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,elastic)"
	}    
	if [info exists AdvGeom($i,atomRemo)] {
	    ############
	    # ATOMREMO #
	    ############
	    xcDebug -debug "ATOMREMO OPTION"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,atomRemo)"
	}    
	if [info exists AdvGeom($i,atomSub)] {
	    ############
	    # ATOMSUBS #
	    ############
	    xcDebug -debug "ATOMSUBS OPTION:: $AdvGeom($i,atomSub)"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,atomSub)"
	}
	if [info exists AdvGeom($i,atomInse)] {
	    ############
	    # ATOMINSE #
	    ############
	    xcDebug -debug "ATOMINSE OPTION:: $AdvGeom($i,atomInse)"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,atomInse)"
	}
	if [info exists AdvGeom($i,atomDisp)] {
	    ############
	    # ATOMDISP #
	    ############
	    xcDebug -debug "ATOMDISP OPTION:: $AdvGeom($i,atomDisp)"
	    # fisrt puts symm flag
	    append geoInput "$AdvGeom($i,atomDisp)"
	}
	if [info exists AdvGeom($i,slab)] {
	    ########
	    # SLAB #
	    ########
	    xcDebug -debug "SLAB OPTION"
	    append geoInput "SLAB\n"
	    append geoInput "$AdvGeom($i,slab)"
	    # we are going from crystal --> slab so nzdir must be 0
	    set nzdir 0
	    set periodic(dim) 2	    
	}
	if [info exists AdvGeom($i,cluster)] {
	    ###########
	    # CLUSTER #
	    ###########
	    xcDebug -debug "CLUSTER OPTION"
	    append geoInput "CLUSTER\n"
	    append geoInput "$AdvGeom($i,cluster)"
	    set periodic(dim) 0
	}    
	if [info exists AdvGeom($i,option)] {
	    #################
	    # manual Option #
	    #################
	    xcDebug -debug "MANUAL OPTION"
	    append geoInput "$AdvGeom($i,option)"
	}
	if [info exists AdvGeom($i,edit)] {
	    ###############
	    # manual EDIT #
	    ###############
	    xcDebug -debug "MANUAL EDIT"
	    # this option is different, as the $AdvGeom($i,edit) holds
	    # the whole input !!!
	    set geoInput "$AdvGeom($i,edit)"
	}

    }
    # next line is experimental !!!
    xcUpdateState

    # handle correctly the EXTPRT/COORPRT/STOP keywords
    set geoInput [cxxHandleEXTPRT $geoInput]

    # debug::
    global periodic
    xcDebug -debug "dim= $periodic(dim), species= $species"    
    return $geoInput
}


proc SetUnDefPar {} {
    global inp 

    # proc 'set' undefined parameters
    if ![info exists inp(A)]        {set inp(A)     {}}
    if ![info exists inp(B)]        {set inp(B)     {}}
    if ![info exists inp(C)]        {set inp(C)     {}}
    if ![info exists inp(ALFA)]     {set inp(ALFA)  {}}
    if ![info exists inp(BETA)]     {set inp(BETA)  {}}
    if ![info exists inp(GAMMA)]    {set inp(GAMMA) {}}
}


