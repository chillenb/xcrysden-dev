#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scriptingFilter.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################




# ------------------------------------------------------------------------
#****c* Scripting/scripting::filter
#
# NAME
# scripting::filter
#
# PURPOSE

# This namespace provide the filter functions for loading crystal
# (molecular) structures from various formats. The filter functions
# are similar to "scripting::open --format file.format". However for
# some formats using the latter form will query interactively some
# parameters, for instance, atomic numbers are usually not known in
# I/O files of pseudopotential-codes and are therefore
# queried. Instead by using the filter functions, the necessary
# information is supplied in function call, thus avoiding interactive
# querying. This makes the filter functions suitable for use in
# scripts. For example, if one wants to make 100 figures, the 100
# interactive queries would be very cumbersome, therefore one can use
# the filter functions instead.

#
# COMMANDS
#
# scripting::filter::g98cube      -- loads the Gaussian98 CUBE file 
# scripting::filter::crystalInput -- loads the CRYSTAL input file 
# scripting::filter::pwscfInput	  -- loads the PWSCF input file 
# scripting::filter::pwscfOutput  -- loads the PWSCF output file 
# scripting::filter::fhiInpini	  -- loads the FHI98MD inp.ini file 
# scripting::filter::fhiCoord     -- loads the FHI98MD coord.out file 
# 
#****
# ------------------------------------------------------------------------

namespace eval scripting::filter {}





# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::g98cube
#
# NAME
# scripting::filter::gXXcube
#
# USAGE
# scripting::filter::gXXcube cube_file mo_index
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display a particular
# molecular orbital (MO) from Gaussian cube-file by specifying the
# "mo_index", which is the ID of the MO as specified in the Gaussian
# cube-file. Using instead the "scripting::open --gXX_cube cube_file"
# proc will result in querying the mo_index. Hence this proc differs
# from "scripting::open --gXX_cube cube_file" in that the mo_index is
# already specified explicitly and is therefore not queried.

#
# ARGUMENTS
# cube_file -- name of Gaussian cube file.
# mo_index  -- which Molecular Orbital to display
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::gXXcube molecular_orbitals.cube 11
#****
# ------------------------------------------------------------------------
proc scripting::filter::gXXcube  {cube_file mo_index} {
    uplevel { scripting::filter::g98cube $cube_file $mo_index }
}
proc scripting::filter::g98cube {cube_file mo_index} {
    global g98

    set cube_file [_init $cube_file]
    
    set files [g98Cube:cube2xsf $cube_file]    
    set mo_file mo-[format %3.3d $mo_index].xsf
    
    set ind [lsearch $files $mo_file]
    if { $ind == -1 } {
	ErrorDialog "molecular orbital with index $mo_index does not exist in cube file: $cube_file"
	exit
    }
    set mo_file [lindex $files $ind]
    xsfOpen $g98(cube_dir)/$mo_file .mesa
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::crystalInput
#
# NAME
# scripting::filter::crystalInput
#
# USAGE
# scripting::filter::crystalInput crystalInput_file
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display a the structure
# from the CRYSTAL input file. Using instead the "scripting::open
# --crystal_inp file" call will result in interactive query of some
# action. Instead this proc immediately displays the structure and the
# interactive query is turned off.

#
# ARGUMENTS
# crystalInput_file -- name of CRYSTAL input file.
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::crystalInput urea.c98
#****
# ------------------------------------------------------------------------

proc scripting::filter::crystalInput {crystalInput_file} {
    global crystalInput

    set crystalInput_file [_init $crystalInput_file]
    
    OpenCyrstalInput $crystalInput_file
    if { [xcIsActive c95] } {
	CalStru
	foreach t $crystalInput(two_toplevels) {
	    destroy $t
	}
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::pwscfInput
#
# NAME
# scripting::filter::pwscfInput -- Filter for PWscf-Input files
#
# USAGE
# scripting::filter::pwscfInput pwscfInput_file reduce ?itypNatList?
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display crystal structure
# from PWSCF input file. It is similar to "scripting::open --pw_inp
# file" call. The latter call will result in querying the atomic
# numbers for PWscf version < 1.3, while this information is supplied
# in the function call by itypNatList argument. One can also specify a
# possible reduction of the structure dimensionality. For example, if
# the structure is a molecule, then one can get rid of 3D unit cell,
# by reducing the dimensionality to 0D.

#
# ARGUMENTS
# pwscfInput_file -- name of PWSCF input file.
# reduce          -- reduce dimenionality to "reduce"-D.
# itypNatList     -- used for PWSCF < 1.2 only.
#                    PWSCF's ityp --> atomic-number (i.e. nat) mapping list.
#                    The format of the list is the following: 
#                    {ityp1 nat1   ityp2 nat2   ...}

#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::pwscfInput water.inp 0 {1 8  2 1} ; # for PWSCF <  1.2
# scripting::filter::pwscfInput water.inp 0            ; # for PWSCF >= 1.2
#****
# ------------------------------------------------------------------------

proc scripting::filter::pwscfInput {pwscfInput_file reduce {itypNatList {}}} {
    global system xcMisc

    # the structure of itypNatList is the following:
    # ----------------------------------------------
    # set itypNatList [list ityp1 nat1   ityp2 nat2]

    # standard PWSCF pre-setting
    set pwscfInput_file [_absoluteFilename $pwscfInput_file]
    _pwscfInputCheckItypNatList $pwscfInput_file $itypNatList
    _pwscf scripting::filter::pwscfInput $pwscfInput_file $reduce $itypNatList

    # reload-capable
    
    reloadRegister $pwscfInput_file openExtStruct:reload

    # load the structure

    openExtStruct 3 crystal external \
	[list sh $system(TOPDIR)/scripts/pwi2xsf.sh] \
	$system(SCRDIR)/pwi2xsf.xsf_out \
	{PWSCF Input File} \
	ANGS \
	-file $pwscfInput_file
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::pwscfOutput
#
# NAME
# scripting::filter::pwscfOutput -- Filter for PWscf-Output files
#
# USAGE
# scripting::filter::pwscfOutput option pwscfOutput_file reduce ?itypNatList?
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display crystal structure
# from PWSCF output file. It is similar to "scripting::open --pw_out
# file" call. The latter call will result in querying the atomic
# numbers, ... Contrary, this information is supplied in the function
# call by itypNatList argument (used only for PWSCF outputs <
# 1.2). One can also specify possible reduction of the structure
# dimensionality. For example, if the structure is a molecule, then
# one can get rid of 3D unit cell, by reducing the dimensionality to
# 0D. Since PWSCF output file may contains several structures, one
# should specify which one to display using the option argument.

#
# ARGUMENTS
# option           -- which structure to render. Possibilities:
#                       --initcoor or -ic   ... render initial structure
#                       --latestcoor or -lc ... latest structure in the file
#                       --optcoor or -oc    ... optimized structure
#                       --animxsf or -a     ... extract all structure (animation)
# pwscfOutput_file -- name of PWSCF output file
# reduce           -- reduce dimenionality to "reduce"-D
# itypNatList      -- used for PWSCF < 1.2 only.
#                     PWSCF's ityp --> atomic-number (i.e. nat) mapping list.
#                     The format of the list is the following: 
#                     {ityp1 nat1   ityp2 nat2   ...}
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::pwscfOutput --optcoor water.out 0 {1 8  2 1} ; # for PWSCF <  1.2
# scripting::filter::pwscfOutput --optcoor water.out 0            ; # for PWSCF >= 1.2
#
#****
# ------------------------------------------------------------------------

proc scripting::filter::pwscfOutput {option pwscfOutput_file reduce {itypNatList {}}} {
    global system xcMisc pw
    # the structure of itypNatList is the following:
    # set itypNatList [list ityp1 nat1   ityp2 nat2]
    
    # standard PWSCF pre-setting

    set pwscfOutput_file [_absoluteFilename $pwscfOutput_file]
    _pwscfOutputCheckItypNatList $pwscfOutput_file $itypNatList
    _pwscf scripting::filter::pwscfOutput $pwscfOutput_file $reduce $itypNatList

    # reload-capable
    global openExtStruct
    set openExtStruct(state)   external
    set openExtStruct(program) [list sh $system(TOPDIR)/scripts/dummy.sh]
    set openExtStruct(program_output) pwo2xsf.xsf
    
    reloadRegister $pwscfOutput_file openExtStruct:reload
    
    # load the structure
    
    switch -exact -- $option {
	--inicoor    - -ic -
	--optcoor    - -oc -
	--latestcoor - -lc -
	--animxsf    - -a {
	    set pw(output_flag) $option
	}
	default {
	    ErrorDialog "wrong option $option, must be one of --initcor, --latestcoor, --optcoor, or --animxsf"
	    exit
	}
    }

    pwOutputPresetWhat dummy $pwscfOutput_file
    set xsf [ReadFile $system(SCRDIR)/pwo2xsf.xsf]
    xsfOpenParsed $system(SCRDIR)/pwo2xsf.xsf .mesa
}


proc scripting::filter::_pwscf {procname file reduce itypNatList} {
    global system
    # the structure of itypNatList is the following:
    # set itypNatList [list ityp1 nat1   ityp2 nat2]

    set len [_pp_common $procname $file $reduce $itypNatList]

    if { $itypNatList != "" } {
	
	# write the nuclei.charges file
        
	set out "[expr $len / 2]\n"
	foreach {ityp nat} $itypNatList {
	    append out "$ityp   $nat\n"
	}
	
    } else {
	# workaround for PWscf version 1.2
	set out 0
    }
    
    evalInScratch {WriteFile nuclei.charges $out w}
}


# ------------------------------------------------------------------------
# scripting::filter::_pwscfInputCheckItypNatList -- 
#
# Check the PWscf Input-file if the itypNatList argument is
# necessary !!!
# ------------------------------------------------------------------------

proc scripting::filter::_pwscfInputCheckItypNatList {file itypNatList} {

    set ntyp         0
    set is_new_input 0
    set is_old_input 0

    foreach line [split [ReadFile $file] \n] {
	foreach field [split $line ,] {
	    if { [string match -nocase *ntyp* $field] } {
		set ntyp [lindex [split $field =] 1]
		break
	    }
	}
	if { [string match -nocase *&input* $line] } {
	    set is_old_input 1
	}
	if { [string match -nocase *&system* $line] } {
	    set is_new_input 1
	}
    }
    puts stderr "File = file"
    if { $ntyp == 0 || (!$is_new_input && !$is_old_input) } {
	# the file is not PW-input file
	ErrorDialog "file \"$file\" is not a PWSCF Input file !!!\n\nTrace:\nntyp=$ntyp\nis_new_input=$is_new_input\nis_old_input=$is_old_input\n\nApplication will Exit."
	exit 1
    }

    if { $is_old_input && $itypNatList == "" } {
	ErrorDialog "please specify itypNatList argument for PWscf Input file < 1.3 !!!\n\nApplication will Exit."
	exit 1
    }
}


# ------------------------------------------------------------------------
# scripting::filter::_pwscfOutputCheckItypNatList -- 
#
# Check the PWscf Output-file if the itypNatList argument is
# necessary !!!
# ------------------------------------------------------------------------

proc scripting::filter::_pwscfOutputCheckItypNatList {file itypNatList} {

    set version 0.0
    set fid [open $file r]
    while { ! [eof $fid] } {
	gets $fid line
	if { [string match "* Program PWSCF *" $line] } {
	    #
	    # get the PWscf version of the output
	    #
	    set ver     [split [string trimleft [lindex $line 2] {v.}] .]
	    set version [lindex $ver 0].[lindex $ver 1]
	    if { [llength $ver] == 3 } {
		append version [lindex $ver 2]
	    }
	    break
	}
    }
    close $fid

    if { $version < 1.2 && $itypNatList == "" } {
	ErrorDialog "please specify itypNatList argument for PWscf Output file < 1.3 !!!\n\nApplication will Exit."
	exit 1
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::fhiInpini
#
# NAME
# scripting::filter::fhiInpini
#
# USAGE
# scripting::filter::fhiInpini fhiInpini_file reduce itypNatList
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display crystal structure
# from FHI98MD inp.ini file. It is similar to "scripting::open
# --fhi_inpini file" call. The latter call will result in querying the
# atomic numbers. Contrary, this information is supplied in function
# call of the scripting::filter::fhiInpini. Also one can specify
# possible reduction of the structure dimensionality. For example, if
# the structure is a molecule, then one can get rid of 3D unit cell,
# by reducing the dimensionality to 0D.

#
# ARGUMENTS
# fhiInpini_file  -- name of FHI98MD input file.
# reduce          -- reduce dimenionality to "reduce"-D
# itypNatList     -- FHI98MD's atomic-name --> atomic-number (i.e. nat) 
#                    mapping list. The format of the list is the following: 
#                    {name1 nat1   name2 nat2   ...}
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::fhiInpini inp.ini 0 {oxygen 8  hydrogen 1}
#****
# ------------------------------------------------------------------------

proc scripting::filter::fhiInpini {fhiInpini_file reduce nameNatList} {
    global system
    # the structure of nameNatList is the following:
    # set nameNatList [list name1 nat1   name2 nat2]
    
    set fhiInpini_file [_absoluteFilename $fhiInpini_file]
    _fhi scripting::filter::fhiInpini $fhiInpini_file $reduce $nameNatList

    # reload-capable    
    reloadRegister $fhiInpini_file openExtStruct:reload

    openExtStruct 3 crystal external34 \
	$system(BINDIR)/fhi_inpini2ftn34 \
	$system(ftn_name).34 {FHI98MD "ini.inp"} \
	BOHR \
	-file $fhiInpini_file
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::filter::fhiCoord
#
# NAME
# scripting::filter::fhiCoord
#
# USAGE
# scripting::filter::fhiCoord fhiCoord_file reduce itypNatList
#
# PURPOSE

# This proc is used in XCRYSDEN scripts to display crystal structure
# from FHI98MD coord.out file. It is similar to "scripting::open
# --fhi_coord file" call. The latter call will result in querying the
# atomic numbers. Contrary, this information is supplied in function
# call of the scripting::filter::fhiCoord. Also one can specify
# possible reduction of the structure dimensionality. For example, if
# the structure is a molecule, then one can get rid of 3D unit cell,
# by reducing the dimensionality to 0D.

#
# ARGUMENTS
# fhiCoord_file   -- name of FHI98MD input file.
# reduce          -- reduce dimenionality to "reduce"-D
# itypNatList     -- FHI98MD's atomic-name --> atomic-number (i.e. nat) 
#                    mapping list. The format of the list is the following: 
#                    {name1 nat1   name2 nat2   ...}
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::filter::fhiCoord coord.out 0 {oxygen 8  hydrogen 1}
#****
# ------------------------------------------------------------------------

proc scripting::filter::fhiCoord {fhiCoord_file reduce nameNatList} {
    global system
    # the structure of nameNatList is the following:
    # set nameNatList [list name1 nat1   name2 nat2]
    
    set fhiCoord_file [_absoluteFilename $fhiCoord_file]
    _fhi scripting::filter::fhiCoord $fhiCoord_file $reduce $nameNatList

    reloadRegister $fhiCoord_file openExtStruct:reload
    
    openExtStruct 3 crystal external \
	$system(BINDIR)/fhi_coord2xcr \
	fhi_coord.xcr {FHI98MD "coord.out"} \
	BOHR \
	-file $fhiCoord_file
}


proc scripting::filter::_fhi {procname file reduce nameNatList} {
    global system xcMisc
    # the structure of nameNatList is the following:
    # set nameNatList [list name1 nat1   name2 nat2]

    set len [_pp_common $procname $file $reduce $nameNatList]
    
    # get the list of atom names
    if { [namespace tail $procname] == "fhiInpini" } {
	xcCatchExecReturn $system(BINDIR)/fhi_inpini2ftn34 $file getlist
	#if { [catch {exec $system(BINDIR)/fhi_inpini2ftn34 $file getlist}] } {
	#    ErrorDialog "error while executing \"fhi_inpini2ftn34 $file getlist\" program"
	#    return 0
	#}
    } else {
	# type == fhiCoord
	xcCatchExecReturn $system(BINDIR)/fhi_coord2xcr $file getlist
	#if { [catch {exec $system(BINDIR)/fhi_coord2xcr $file getlist}] } {
	#    ErrorDialog "error while executing \"fhi_coord2xcr $file getlist\" program"
	#    return 0
	#}
    }

    # get the atom names
    set atomNames [lrange [ReadFile $system(SCRDIR)/fhi_species_name.list] 1 end]

    # write the nuclei.charges file
    
    set out "[expr $len / 2]\n"
    foreach {name nat} $nameNatList {
	set Nat($name) $nat
    }
    # the order in nameNatList and atomNames might be different, the
    # atomNames order is the right one
    foreach name $atomNames {
	append out $Nat($name)\n
    }

    # ------------------------------------------------------------------------
    # cd to SCRATCH directory and write the file
    cd $system(SCRDIR)
    WriteFile fhi_species_nat.list $out w
    # ------------------------------------------------------------------------
}


# pp_common == pseudo-potential common
proc scripting::filter::_pp_common {procname file reduce atomMappingList} {
    global xcMisc

    set file [_init $file]
    
    # REDUCE factor is used for reducing the periodic structure to
    # lower dimensions
    
    if { ! [string is digit $reduce] } {
	ErrorIn $procname "reduce-factor must be an integer, but got $reduce"
    } else {
	if { $reduce > -1 && $reduce < 4 } {
	    set xcMisc(reduce_to) $reduce
	}
    }
    
    # check the length itypNatList
    
    set len [llength $atomMappingList]
    if { $len % 2 } {
	ErrorIn $procname "odd number of elements in atomMappingList list"
	exit
    }
    return $len
}


proc scripting::filter::_init {file} {
    
    # pop-up the Viewer

    ViewMol .
    
    # destroy the welcome window
    destroyWelcomeWindow 
    
    set file [_absoluteFilename $file]

    if { ! [file exists $file] } {
	ErrorDialog "File $file does not exist"
	exit
    }
    return $file
}


proc scripting::filter::_absoluteFilename {file} {
    if { [file pathtype $file] != "absolute" } {
	global system
	evalInPWD {
	    set file [file join [pwd] $file]
	}
    }
    return $file
}
