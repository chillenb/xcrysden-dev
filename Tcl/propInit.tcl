#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/propInit.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

###############################################
# prop(unit9) .. name of Crystal95's unit 9 - it's path
# prop(unit25_1).name of first Crystal95's unit25 on it's scratch directory
# prop(unit25_2).name of second Crystal95's unit25
# prop(dir) .... name of directory, where fortran units are stored
# prop(file) ... name of fortran unit (without possix)
# prop(n_band) . number of bands
# prop(band) ... selected bands 
# prop(firstband) ... first band for BWID
# prop(lastband)  ... last band for BWID
# prop(newk)      ... 0 = newk not yet set; 1 = newk already set
# prop(newk_script) . newk's script
# prop(doss_criteria) "band-interval criteria" or "energy-interval criteria" 
# prop(n_atom) ...... number of atoms per shell
# prop(type_of_run) . RHF/UHF
##########
# prop(NPY)    ...... number of point along B-A segment (used by C95 by 2D grid
#                     points evaluation
# prop(NPZ)    ...... number of point along third segment (used by C95 by 
#                     3D grid points evaluation
# prop(datagridDim) ..... dimensionality of datagrid
# prop(isolevel)..... isolevel to extract in 3D isosurface evaluation
# prop(pm_isolevel).. whether to render +/- of isolevel
##########
# prop(dif_prop3D_list) .... list of properties that can be rendered as 
#                            isosurface difference maps
# prop(spin_prop_list)  .... list of spin dependent properties
# prop(spin_case)       .... 0/1; we may deal with UHF & spin-undependent 
#                            property and so prop(spin_case) is zero
# prop(PATO_newbasis_A) .... when we select density matrix to be superposition
#                            of atomic densities, than we can specify new basis
#                            set 
# prop(PATO_newbasis_B) .... when using difference maps or properties on 
#                            surfaces we have to cards (A & B)
##########
# prop(c95_BAND_script) .... script for band structure for c95's properties
# prop(NLINE)           .... # of lines in reciprocal space to be explored
###############################################
###############################################
# dif_isosurf(dif_map)              0/1; 1 when diffrential map is on
# dif_isosurf(denmat_A)             what density matrix_A to take 
# dif_isosurf(file_textvar_A)       label_text for left frame; 
# dif_isosurf(prop_A)               properties A
# dif_isosurf(spin_A)               spin A
########## 
# dif_isosurf(load_another_unit9)   1 if another unit9 is reqursted; else 0
# dif_isosurf(unit9_B)              another unit9
# dif_isosurf(file_textvar_B)       label text for right frame
# dif_isosurf(prop_B)               properties B
# dif_isosurf(denmat_B)             what density matrix_B to take
# dif_isosurf(spin_B)               spin B
###############################################
# for definition of ISOSURF array luks in "grid.tcl" file
###############################################


##############################################################################
# THIS PROC IS CALLED WHEN WE OPEN CRYSTAL PROPERTIES
##############################################################################

proc PropC95 {{path {}}} {
    global prop system err periodic nxdir nydir nzdir \
	    dif_isosurf xcMisc

    # some initializations

    PropC95Init
    GraphInit

    if { $path == "" } {
	set path [fileselect {Open CRYSTAL's unit 9}]
    }

    if { $path == "" } {
	return
    }

    xcAppendState c95
    xcAppendState properties
    xcUpdateState

    #
    # update the title of toplevel
    #
    wm title . "XCrySDen: [file tail $path]"
    set xcMisc(titlefile) [file tail $path]

    set prop(unit9) $path
    append dif_isosurf(file_textvar_A) $path; #now we completed we element
    # this is just for now
    set prop(dir) $system(SCRDIR)
    
    ########################################
    # CD to $system(SCRDIR)
    cd $system(SCRDIR)
    
    #catch {file delete -force __fort_unit}
    #file mkdir __fort_unit
    #cd __fort_unit
    
    # find out what is the name of fortran units (without number)    
    set prop(file) [FtnName].
    puts stderr $prop(file)
    
    cd $system(SCRDIR)

    # now copy prop(unit9) to $prop(dir)/$prop(file)9
    file copy -force $prop(unit9) $system(SCRDIR)/$prop(file)9

    # now we will get some info about \"case\" and we will render a structure
    set input "BASE\n0\nEXTPRT\nCOORPRT\nEND"
    if { ! [RunC95 $system(c95_properties) {} $input {} {} $system(SCRDIR)] } {
	#tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	return
    }
    
    GetC95Info EXTPRT $system(SCRDIR)/xc_struc.$system(PID)
    GetC95Info TYPE_OF_RUN $system(SCRDIR)/xc_output.$system(PID)

    puts stderr "Dimension of system: $periodic(dim)"

    # maybe structure is not 3 dimensional
    if { $periodic(dim) < 3 } { set nzdir 0 }
    if { $periodic(dim) < 2 } { set nydir 0 }
    if { $periodic(dim) < 1 } { set nxdir 0 }

    xcAppendState render
    xcUpdateState
    GenGeomDisplay
    OpenXSFStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
}


proc PropC95Init {} {
    global prop dif_isosurf isosurf openGL

    if [array exists prop]        { unset prop }
    if [array exists dif_isosurf] { unset dif_isosurf }

    # initialize prop array
    set prop(newk) 0
    set prop(NPRO) 0
    set prop(PATO_newbasis_A) 0
    set prop(PATO_newbasis_B) 0
    set prop(datagridDim)     0

    #set prop(dif_prop3D_list) { \
    #	     {CHARGE DENSITY} \
    #	     {CHARGE DENSITY GRADIENT} \
    #	     {ELECTROSTATIC POTENTIAL} }

    set prop(dif_prop3D_list) { \
	    {CHARGE DENSITY} \
	    {ELECTROSTATIC POTENTIAL} }

    set prop(spin_prop_list) {ECHD "ECHG\n0"}
    set prop(isolevel)       {}
    set prop(pm_isolevel)    0

    # initialize dif_isosurf array
    set dif_isosurf(denmat_A)            "SCF density matrix"
    set dif_isosurf(load_another_unit9)  0
    set dif_isosurf(file_textvar_A)      "default UNIT 9:\n"
    set dif_isosurf(denmat_B) \
	    "Density matrix as superposition of atomic densities"
    set dif_isosurf(file_textvar_B)      "loaded UNIT 9:\nsame as in Map A"

    # initialization of dif_isosurf array
    SetIsoSurfArray
    ConvertTwoSideVar; # this is to switch isosurf(twoside_lighting) from 0/1 -> on/off
}


proc GetC95Info {whatlist file {dir {}}} {
    global system periodic geng prop

    if { $dir == {} } {
	########################################
	# CD to $system(SCRDIR)
	cd $system(SCRDIR)
	########################################
    } else {
	cd $dir
    }

    # whatlist:: lists of "whats" type of info needed
    foreach what $whatlist {
	if { $what == "EXTPRT" } {
    # usage of "gengeom" program:
    # 
    # gengeom  MODE1  MODE2  MODE3  IGRP  NXDIR  NYDIR  NZDIR  OUTPUT  INPUT
    #    0       1      2      3      4     5      6      7      8       9
    #
    # FIND FROM UNIT34 THE DIMENSIONALITY OF THE SYSTEM

	    if { $system(c95_version) == "95" } {
		set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
	    } elseif { $system(c95_version) == "98" } {
		set geng(M3_AGRUMENT) [GetGengM3Arg ANGS 98]
	    } elseif { $system(c95_version) == "03" } {
		set geng(M3_AGRUMENT) [GetGengM3Arg ANGS 03]
	    } elseif { $system(c95_version) == "06" } {
		set geng(M3_AGRUMENT) [GetGengM3Arg ANGS 06]
	    } elseif { $system(c95_version) == "09" } {
		set geng(M3_AGRUMENT) [GetGengM3Arg ANGS 09]
	    } elseif { $system(c95_version) == "14" } {
		set geng(M3_AGRUMENT) [GetGengM3Arg ANGS 14]
	    } 	    
		
	    set file $system(SCRDIR)/xc_gengeom.$system(PID)
	    xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) 1  1 1 1 $file

	    set fileID [open $file r]
	    GetDimGroup periodic(dim) periodic(igroup) $fileID
	    close $fileID	    
	    
	    puts stdout "dim igroup:: $periodic(dim) $periodic(igroup)"
	}
	if { $what == "BWID" || $what == "DOSS" || $what == "BAND" } {
	    # find out the number of bands
	    set fileID [open $file r]
	    foreach line [split [read $fileID] \n] {		
		if [string match "*NUMBER OF AO*" $line] {
		    set prop(n_band) [lindex $line 3]
		}
		if [string match "*N\. OF ATOMS PER CELL*" $line] {
		    set prop(n_atom) [lindex $line 5]
		}
	    }	
	    close $fileID
	}
	if { $what == "TYPE_OF_RUN" } {
	    # allowed prop(type_of_run) values: UHF/RHF
	    set fileID [open $file r]
	    foreach line [split [read $fileID] \n] {		
		if { [string match "*TYPE OF CALCULATION*" $line] } {
		    if { $system(c95_version) == "03" || $system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14" } {
			#
			# CRYSTAL03/06/09/14
			#
			if { [string match "*UNRESTRICTED*" $line] } {
			    set prop(type_of_run) UHF
			} else {
			    set prop(type_of_run) RHF
			}
		    } else {
			#
			# CRYSTAL-95/98
			#
			set prop(type_of_run) [lindex $line end]
		    }
		    xcDebug -stderr "Type of CRYSTAL run: $prop(type_of_run)"
		}
	    }
	}
    }
}
