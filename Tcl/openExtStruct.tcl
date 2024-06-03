#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/openExtStruct.tcl                                #
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################
proc openExtStruct:reload {filedir can {update 1}} {
    global openExtStruct system

    if { ! [info exists openExtStruct(state)] } {
        error "openExtStruct:reload called before openExtStruct"
    }    
    set file [gunzipFile $filedir]
    set program_output $openExtStruct(program_output)
        
    if { $program_output != "stdout" } {
        if { $program_output != "pwo2xsf.xsf" } {
            eval xcCatchExecReturn $openExtStruct(program) $file
        } else {
            # a bit special case for PWscf output files
            eval xcCatchExecReturn $openExtStruct(program) $file > $program_output
        }
    } else {
        set program_output $system(SCRDIR)/xc_out.$system(PID)
        eval xcCatchExecReturn $openExtStruct(program) < $file > $program_output
    }
    
    if { $openExtStruct(state) != "external34" } {
        # default
        xsfOpenParsed $program_output $can 1
    } else {
        global geng
        
        xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) 1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID)
        GetDimGroupXSF periodic(dim) periodic(igroup) $system(SCRDIR)/xc_gengeom.$system(PID)
        
        CellMode 1
        xcUpdateState
    }
}

proc openExtStruct {dim spec state program program_output file_type unit {args {}}} {
    global fileselect system periodic species geng xcMisc wn openExtStruct

    #
    # dim     
    #         ... dimension of the structure
    # spec    
    #         ... structure (i.e. species) name (molecule/polymer/slab/crystal)
    # state   
    #         ... could be either external or external34
    # program 
    #         ... name of the converting program (full-path name)
    # program_output 
    #         ... output file-name created by converting program; 
    #             
    #             if $program_output==stdout, then converting
    #             program reads from stdin and writes to stdout !!!
    #
    # file_type 
    #         ... type of file; it is just a name like "PWSCF Input File"
    # unit
    #         ... lenghth unit (ANGS/BOHR) for $program_oputput
    # args 
    #         ... additional arguments like:
    #
    #     -preset ... presetting procedure (used for pseudopotential codes,
    #                 where atomic numbers should be queried: see files
    #                 fhiPreset.tcl & pwPreset.tcl)
    #     -file   ... filename to open; used when launched as command line
    #                 option
    #

    set openExtStruct(state)   $state
    set openExtStruct(program) $program
    set openExtStruct(program_output) $program_output
    
    # 
    # internal checking of parameters  !!!
    #
    if { $unit != "ANGS" && $unit != "BOHR" } {
	xcDebug -stderr "openExtStruct: Syntax Error -- unit must be ANGS or BOHR"
	exit 1
    }
    if { $state != "external" && $state != "external34" } {
	xcDebug -stderr "openExtStruct: Syntax Error -- state must be external or external34"
	exit 1
    }
    if { $dim < 0 || $dim > 3 } {
	xcDebug -stderr "openExtStruct: Syntax Error -- dim must be within \[0,3\]"
	exit 1
    }

    ################
    # parsing args #
    ################
    set i 0
    set file   {}
    set preset {}    
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
	    xcDebug -debug "openExtStruct Options:: $tag $option"
            switch -- $tag {
                "-file"   {set file   $option}
                "-preset" {set preset $option}
		default { 
		    tk_dialog .mb_error Error "ERROR: Bad \"openExtStruct\" configure option $tag" error 0 OK 
		    return 0
		}
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error "ERROR: You called openExtStruct with an odd number of args !" error 0 OK
	return 0
    }

    #
    # check/set $file
    #
    if { $file == {}} {
        global system
        set file [tk_getOpenFile -defaultextension .in \
                      -filetypes { 
                          {{All Files}        {.*}}
                          {{Input Files}      {.in .inp}} 
                          {{Output Files}     {.out}} } \
                      -initialdir $system(PWD) \
                      -title "Open $file_type File"]
	if { $file == "" } {
	    return
	}

        # active "Reload button"
        reloadRegister $file openExtStruct:reload
    }
    
    if { ![file exists $file] } {
	tk_dialog .update \
            "WARNING !!!" "WARNING: File \"$file\" does not exist !!!" \
            warning 0 OK
	return
    }

    #
    # update the title of "."
    #
    set xcMisc(titlefile) $file
    wm title . "XCrySDen: [file tail $file]"
    
    set file [gunzipFile $file]

    # *****************
    # check/set $preset
    # *****************
    # RETURN status:
    #                0 ... CANCEL button was pressed OR an ERROR occured --
    #                      -- (should be taken care within $preset routine)
    #                1 ... $preset routine ended successfully
    if { $preset != {} } {
	if { ![eval $preset $file] } { 
	    return 0 
	}
    }

    # for GenGeom to know from which file to read XSF
    set xcMisc(external_xsf_name) $program_output

    set periodic(dim) $dim   
    set species       $spec
 
    file copy -force $file $system(SCRDIR)/xc_${state}.$system(PID) 

    ##############################
    cd $system(SCRDIR)
    ##############################
    # execute CONVERTING program #
    ##############################
    if { $program_output != "stdout" } {
	eval xcCatchExecReturn $program $file
    } else {
	set program_output $system(SCRDIR)/xc_out.$system(PID)
	eval xcCatchExecReturn $program < $file > $program_output
    }

    # check if $program_output was created ?!!
    if { ! [file exists $program_output] } {
	tk_dialog [WidgetName] ERROR "openExtStruct ERROR:: $program_output file does not exist !!!" error 0 OK
	return
    }

    set geng(M3_ARGUMENT) [GetGengM3Arg $unit 95]
    
    #
    # determine periodic(igroup)
    #
    if { $state == "external34" } {
	xcAppendState $state
	xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) \
		1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID)
	
	set fileID [open "$system(SCRDIR)/xc_gengeom.$system(PID)" r]
	GetDimGroup periodic(dim) periodic(igroup) $fileID
	close $fileID
	
	ResetDispModes
	CellMode
	OpenXSFStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
	set light On
	Lighting On
	# append 'render' to XCState(state) if "render" is not defines yet
	xcAppendState render
	xcUpdateState
	return	
    } else {
	# state == external (means the format of $program_output is XSF)
	xsfOpenParsed $program_output .mesa
	return
	
	#xcCatchExecReturn $system(BINDIR)/gengeom $geng(M1_INFO) 1 $geng(M3_ARGUMENT) \
	#	1  1 1 1 $system(SCRDIR)/xc_gengeom.$system(PID) \
	#	$program_output
    }
}
