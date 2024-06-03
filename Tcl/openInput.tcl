#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/openInput.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc OpenCyrstalInput {{file {}}} {
    global fileselect distext Alist species speciesName \
	    type_group type_group1 job_title \
	    inp n_groupsel groupsel crdatom XCState system \
	    AdvGeom xcMisc crystalInput

    if { $system(c95_exist) == 0 } {
	ErrorDialog "can't open CRYSTAL-95/98/03 input File. CRYSTAL package is not installed !!!"
	return
    } 

    # distext .... here go informations to be displayed in 
    #              information text widget

    
    if { $file == "" } { 
	#fileselect "Open CRYSTAL Input"
        set fileselect(path) [tk_getOpenFile -defaultextension .r1 \
                                  -filetypes { 
                                      {{All Files}           {.*}}
                                      {{CRYSTAL Input Files} {.r1}}
                                  } -initialdir $system(PWD) -title "Open CRYSTAL Input"]
        
    } elseif { [file isdirectory $file]} {
	set fileselect(path) [tk_getOpenFile -defaultextension .r1 \
				  -filetypes { 
				      {{All Files}           {.*} }
				      {{CRYSTAL Input Files} {.r1}}
				  } -initialdir $file \
				  -title "Open CRYSTAL Input"]
	if { $fileselect(path) == "" } {
	    return
	}
    } else {
	set fileselect(path) $file
    }

    # maybe CANCEL button was pressed
    if { $fileselect(path) == "" } {
	xcDeleteState c95
        xcDeleteState openinput
        xcUpdateState
	return
    }

    #################
    # initialisation
    set XCState(state) c95_openinput
    xcUpdateState
    xcAdvGeomState reset

    #
    # reset the title of "." 
    #
    wm title . "XCrySDen: [file tail $fileselect(path)]"
    set xcMisc(titlefile) [file tail $fileselect(path)]


    # OK button was pressed
    # check if selected file is Crystal95 file;
    # the best way for doing it is to go and check it with Crystal95
    
    ########################################
    # CD to $system(SCRDIR)
    cd $system(SCRDIR)
    ########################################

    # test only geom part of the input file (this is quick)
    # xc_inp.$system(PID)...here just geom input will be stored
    xcCatchExecReturn $system(AWKDIR)/ginp.awk \
	    $fileselect(path) > $system(SCRDIR)/xc_inp.$system(PID)

    puts stdout "FILE: $fileselect(path)"
    puts stdout "FILTERED INPUT"
    catch {ReadFile $system(SCRDIR)/xc_inp.$system(PID)}

    # BEGIN t.k.
    # for EXTERNAL: copy also fort.34 
    if { [file exists $system(PWD)/$system(ftn_name).34 ] } {
	file copy -force $system(PWD)/$system(ftn_name).34 $system(SCRDIR)/external_unit34
	file copy -force $system(SCRDIR)/external_unit34   $system(SCRDIR)/$system(ftn_name).34
    }
    # END t.k.

    # if we catch error than selected file is of the right type,
    # but is corrupted
    if { [catch {exec $system(c95_integrals) < \
		     $system(SCRDIR)/xc_inp.$system(PID) >& \
		     $system(SCRDIR)/xc_tmp.$system(PID)} errmsg] } {
	set idx [tk_dialog .idx1 ERROR "Selected file seems to be \
		Crystal95 input file, but is corrupted" error 0 OK Details]
	if { $idx == 1 } {
	    #user want's to see details
	    tk_dialog .errm Details "ERROR MESSAGE:\n$errmsg" {} 0 OK
	}
	catch {file delete $system(SCRDIR)/xc_tmp.$system(PID)}
	OpenCyrstalInput
	return	
    } else {
	# if Crystal95 find out that file is "bad", it has exited
	# nicely, but with ERROR message; grep ERROR --> if grep
	# doesn't find anything, we must catch the grep error

	#eval [list exec grep ERROR $system(SCRDIR)/xc_tmp.$system(PID)] 	
	if { ! [catch {exec grep ERROR $system(SCRDIR)/xc_tmp.$system(PID)}] } {
	    set is_error 1
	    global system
	    if { $system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14" } {
		# if the "STOP KEYWORD - EXECUTION STOPS" is found then everything is OK"
		if { ! [catch {exec grep "STOP KEYWORD - EXECUTION STOPS" $system(SCRDIR)/xc_tmp.$system(PID)}] } {
		    # OK, we don't have error
		    set is_error 0
		}	
	    }
	    if { $is_error } {
		puts stdout "grep ERROR catched"
		tk_dialog .idx2 ERROR "Selected file is bad !!" error 0 OK
		file delete $system(SCRDIR)/xc_tmp.$system(PID)
		OpenCyrstalInput
		return
	    }
	}
	
    }
    # it Looks that Selected file is good !!!

    # READ THE FILE; "distext" variable collects a information to be displayed
    set input [open "$system(SCRDIR)/xc_inp.$system(PID)"]
    set job_title [gets $input]
    set distext "> TITLE::\n$job_title\n"
    append distext "--------------------------------------------------\n\n" 
    # just in any case (lindex 0) if there is anything bisides the "$species"
    set species [string tolower [lindex [gets $input] 0]]    
    append distext "> SPECIES::  $species\n"
    append distext "--------------------------------------------------\n\n"

    if { $species == "external" } {
	# EXTERNAL OPTION
	CalStru
	return
    } elseif { $species == "crystal" } {
	# ==========================
	# SPECIES == CRYSTAL
	# ==========================

	set type_group "space"
	set type_group1 "Space"
	set ifl [gets $input]
	set inp(IFLAG) [lindex $ifl 0]
	set inp(IFHR) [lindex $ifl 1]
	set inp(IFSO) [lindex $ifl 2]
	append distext "> CRYSTAL FLAGS::\n"
	append distext "IFLAG:  $inp(IFLAG),    IFHR:  $inp(IFHR),    IFSO:  $inp(IFSO)\n"
	append distext "--------------------------------------------------\n\n"
	# ===================================================
	# WHAT ABOUT sequ. number or alfanum. code for "group"
	if { $inp(IFLAG) == 0 } {
	    # n_groupsel is synonym for IGR; (lindex 0) is just in any case
	    set n_groupsel [lindex [gets $input] 0]
	    set groupsel [Igr2Agr $n_groupsel space_group]
	    append distext "> SPACE GROUP::\n"
	    append distext "IGR = $n_groupsel  -->  AGR = $groupsel\n"
	    append distext "--------------------------------------------------\n\n"
	
	} else {
	    # $groupsel & $AGR are synonyms
	    set inp(AGR) [gets $input]
	    set n_groupsel [Agr2Igr $inp(AGR)]
	    set groupsel $inp(AGR); #gropusel & AGR are synonyms
	    append distext "> SPACE GROUP::\n"
	    append distext "AGR = $groupsel  -->  IGR = $n_groupsel\n\n"
	    append distext "--------------------------------------------------\n\n"
	}
	if { $inp(IFSO) > 1} {
	    # non-standard shift of the ORIGIN
	    set ixyz [gets $input]
	    set inp(IX) [lindex $ixyz 0]
	    set inp(IY) [lindex $ixyz 1]
	    set inp(IZ) [lindex $ixyz 2]
	    append distext "> NON-STANDARD ORIGIN SHIFT::\n"
	    append distext "IX = $inp(IX),    IY = $inp(IY),    IZ = $inp(IZ)\n"
	    append distext "--------------------------------------------------\n\n"	    
	}
	# verify which unit-cell parameter must be read & read it !!!
	WhichPar2Read $input 
    } elseif { $species == "slab" } {
	set type_group "plane"
	set type_group1 "Plane"
	# n_groupsel is synonym for IGR; (lindex 0) is just in any case
	set n_groupsel [lindex [gets $input] 0]
	set groupsel [Igr2Agr $n_groupsel plane_group]
	append distext "> LAYER GROUP::\n"
	append distext "IGR = $n_groupsel  -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
	# verify which unit-cell parameter must be read & read it !!!
	WhichPar2Read $input 
    } elseif { $species == "polymer" } {
	set type_group "line"
	set type_group1 "Line"
	# n_groupsel is synonym for IGR; (lindex 0) is just in any case
	set n_groupsel [lindex [gets $input] 0]
	set groupsel [Igr2Agr $n_groupsel line_group]
	append distext "> ROD GROUP::\n"
	append distext "IGR = $n_groupsel -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
	# for all polymers we must read just A parameter
	set inp(A) [gets $input]
	append distext "> UNIT CELL PARAMETER::\nA:  $inp(A)\n"
	append distext "--------------------------------------------------\n\n"
    } elseif { $species == "molecule" } {
	set type_group "point"
	set type_group1 "Point"
	# n_groupsel is synonym for IGR; (lindex 0) is just in any case
	set n_groupsel [lindex [gets $input] 0]
	set groupsel [Igr2Agr $n_groupsel point_group]
	append distext "> POINT GROUP::\n"
	append distext "IGR = $n_groupsel  -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
    }
	
    # ===============================================
    # THIS IS COMMON FOR ALL SPECIES
    # inp(NATR); (lindex 0) is just in any case to be save
    set inp(NATR) [lindex [gets $input] 0]
    # crdatom is used for checking the variables
    set crdatom 1
    append distext "> NUMBER OF NON-EQUIVALENT ATOMS::\nNATR = $inp(NATR)\n\n"
    append distext "> ATOMIC NUMBERS & COORDINATES OF NON-EQUIVALENT ATOMS::\n"
    # read Nat,X,Y,Z
    for {set i 1} {$i <= $inp(NATR)} {incr i} {
	set natrx [gets $input]
	set inp(NAT,$i) [lindex $natrx 0]
	set inp(X,$i) [lindex $natrx 1]
	set inp(Y,$i) [lindex $natrx 2]
	set inp(Z,$i) [lindex $natrx 3]
	# this is to load atom names
	AtomNames
	append distext [format "%-3d %-4s %10.5f %10.5f %10.5f\n" \
	$inp(NAT,$i) [Nat2Aname $inp(NAT,$i)] $inp(X,$i) $inp(Y,$i) $inp(Z,$i)]
    }
    append distext "--------------------------------------------------\n\n"

    set speciesName $species

    #######################################
    # check for ADVANCE GEOMETRICAL INPUT #
    #######################################
    set line [gets $input]
    if { $line != "STOP" } {
	#############################
	# ADVANCE GEOMETRICAL INPUT #
	#############################	
	set AdvGeom(input) "$line\n"
	while { [set line [gets $input]] != "STOP" } {
	    append AdvGeom(input) "$line\n"
	}
	xcDebug "\nADVANCE GEOMETRICAL INPUT FOUND::\n$AdvGeom(input)\n"
	append distext "     #########################################\n"
	append distext "     # GEOMETRY MANIPULATION INPUT was found #\n"
	append distext "     #########################################\n\n"
	append distext "Geometry manipulation input::\n"	
	append distext "$AdvGeom(input)\n"
	# chack if dimensionality of the system has changed during 
        # the geometry manipulation
	foreach word $AdvGeom(input) {
	    switch -glob -- $word {
		*SLAB*       { set speciesName "slab" }
		*MOLECULE*   { set speciesName "molecule" }
		*CLUSTER*    { set speciesName "cluster" }
	    }
	}
    }

    # ============================================================
    #   INPUT FILE HAS BEEN READ OUT !!!!!!!!
    #
    #   make a toplevel where some information about selected file 
    #   will be printed out
    #
    #   there will be text widget & OK button
    #   for text widged there is a DispText procedure; 
    #   all will be displayed after the input file will be read out
    #
    #   produce some toplevel where it will
    #   be possible to modify different parameters
    # =============================================================

    # StatusWidget creates Status Widget & return path of toplevel
    set tx [StatusWidget]

    # TOPLEVEL FOR DECISION (VIEWER,MODIFY)
    set td [xcToplevel .opfd "Open Crystal Input" "Open Crystal Input" \
	    . 50 100]
    set crystalInput(two_toplevels) {.opfd .opftx}
    AlwaysOnTopON . $crystalInput(two_toplevels)
    focus $td
 
    set l [label $td.lbl -text "What to do?" -relief raised -bd 2]
    set f [frame $td.frm]
    set b1 [button $f.but1 -text "Modify File" -command \
	    [list OpenCyrstalInputModify $tx $td]]
    puts stdout "tx> $tx"
    set b2 [button $f.but2 -text "View $speciesName" -command \
	    [list OpFile2ViewMol $tx $td]]
    pack $l $f -side top -expand 1 -fill both -ipadx 10 -ipady 10
    pack $b1 $b2 -side left -expand 1 -padx 7
}


proc OpenCyrstalInputModify {tx td} {
    global fileselect distext species type_group type_group1 \
	    inp n_groupsel groupsel XCState XCTrace

    CancelProc $td
    # MODIFY TOPLEVEL
    set t [xcToplevel .openfile "Modify" "Modify" \
	    . 50 100]
    # .opftx -- status widget
    AlwaysOnTopON . {.openfile .opftx}
    focus $t

    # because what ever I choose there will apper some new toplevel, that will
    # override AlwaysOnTop flag for .openfile & .opftx toplevels, 
    # so we must set
    #set XCState(toplevel) {.openfile .opftx}

    puts stdout "species: $species"
    flush stdout
    # label goes on the top
    set l [label $t.lbl -text "MODIFY/CHANGE:" -relief groove -bd 2]
    pack $l -side top -expand 1 -fill x -padx 7 -pady 7 -ipady 7 -ipadx 10
    # for every option/parameter make a button
    set spe [button $t.b1 -text "Species" \
	    -command ChooseSpecies]
    # if you will change $species you must change name of this button
    # this is possible throuh buutn-entry combination
    set igr [button $t.b4 -text "Group" \
	    -command [list ModGroup $t]]
    # Only for rhombohedral group is not meanningless to specify IHFR
    # so doit by CheckGroup proc
    set ifhr [button $t.b2 -text "Type of Cell for \n\
	    Rhombohedral Groups" -command CheckGroup]
    set ifso [button $t.b3 -text "Origin Setting" \
	    -command PreSetOrigin]
    set par [button $t.b5 -text "Cell Parameters" \
	    -command [list PreGeom_sym_input .opflgeom $t]]
    set coor [button $t.b6 -text "Atomic Coordinates &\n\
	    Atomic Numbers" -command atom_num_coord]
    set view [button $t.b8 -text "View Structure" \
	    -command [list OpFile2ViewMol $tx $t]]
    set close [button $t.b9 -text "Close" \
	    -command [list DestroyOpfl $t $tx]]
    
    if { $species == "crystal" } {
	pack $spe $igr $ifhr $ifso $par $coor $view $close \
		    -fill x -expand 1 -padx 5 -pady 3 -ipadx 0 -ipady 0
	if { [lindex $groupsel 0] != "R" } {
	    $ifhr config -state disabled
	    set XCTrace(RHOMBO_TYPE_BUTTON) $ifhr
	    trace variable groupsel w xcTrace
	}
    } elseif { $species == "slab" || $species == "polymer"} {
	pack $spe $igr $par $coor $view $close -fill x -expand 1 \
		-padx 5 -pady 3 -ipadx 0 -ipady 0
    } else {
	pack $spe $igr $coor $view $close -fill x -expand 1 \
		-padx 5 -pady 3 -ipadx 0 -ipady 0
    }
}


proc cxxManageExternal {} {
    global species system

    # if the species is external we need to copy
    # $system(SCRDIR)/external_unit34 to
    # $system(SCRDIR)/$system(ftn_name).34
    # (see proc OpenCyrstalInput)
    
    xcDebug -stderr "cxxManageExternal:: species = $species"

    if { $species == "external" } {
	file copy -force $system(SCRDIR)/external_unit34  $system(SCRDIR)/$system(ftn_name).34
    }
}


proc Agr2Igr {agr} {
    global group_list

    # "load" a $group_list
    space_group
    set n 1
    # assign a sequ. number to $igr that correspond to $agr
    foreach word $group_list {
	# pure the $word
	set last [ string length $word ]
	set word [ string range $word 5 $last ]
	regexp {(([A-Z0-9] )|[A-Z0-9\/\-])+[A-Z0-9]} $word word	
	if { $agr == $word } {
	    set igr $n
	}
	incr n
    }
    
    puts stdout "Agr2Igr> agr = $agr"
    flush stdout
    # puts stdout "         igr = $igr"
    # maybe "agr"symbol is not standard one
    if ![info exists igr] { return "\"$agr\" is not a standard space group" }
    return $igr
}


proc Igr2Agr {igr comm} {
    global group_list

    # "load" a $group_list
    eval $comm
    set n 1
    # assign asymbol to $agr that correspond to $igr
    foreach word $group_list {
	# purify the $word
	if { $igr == $n } {
	    set last [ string length $word ]
	    set word [ string range $word 5 $last ]
	    regexp {(([A-Z0-9] )|[A-Z0-9\/\-])+[A-Z0-9]} $word word	
	    set agr $word
	}
	incr n
    }
    return $agr
}


proc OpFile2ViewMol {tx td} {
    # there maybe trace on "groupsel" variable
    xcTraceDelete groupsel
    CancelProc $tx
    CancelProc $td

    # CHECK THE VARIABLES --> variables must be checked everytimes- 
    # we do sometning in Modify; so if we are here everything is OK
    CalStru
    return
}


proc ModGroup {t} {
    global species

    if { $species == "molecule" } {
	# load point groups
	point_group
	geom_sym_input .opflmod $t
    } elseif { $species == "polymer" } {
	# load groups
	line_group
	geom_sym_input .opflmod $t
    } elseif { $species == "slab" } {
	# load groups
	plane_group
	crys_slab_sym .opflmod $t
    } elseif {$species == "crystal" } {
	# load groups
	space_group
	crys_slab_sym .opflmod $t
    }
}


proc UpdateStatus {t} {
    global fileselect distext Alist species type_group type_group1 job_title \
	    inp n_groupsel groupsel crdatom
    #t....text-widget-path

    # put updated information in $distext
    set distext "> TITLE::\n$job_title\n"
    append distext "--------------------------------------------------\n\n" 
    append distext "> SPECIES::  $species\n"
    append distext "--------------------------------------------------\n\n"
    # ==========================
    # SPECIES == CRYSTAL
    # ==========================
    if { $species == "crystal" } {
	append distext "> CRYSTAL FLAGS::\n"
	append distext "IFLAG:  $inp(IFLAG),    IFHR:  $inp(IFHR),    IFSO:  $inp(IFSO)\n"
	append distext "--------------------------------------------------\n\n"
	if { $n_groupsel == 999 } {
	    set n_groupsel "\"$groupsel\" is not a standard space group"
	}
	if { $inp(IFLAG) == 0 } {
	    append distext "> SPACE GROUP::\n"
	    append distext "IGR = $n_groupsel  -->  AGR = $groupsel\n"
	    append distext "--------------------------------------------------\n\n"
	} else {
	    append distext "> SPACE GROUP::\n"
	    append distext "AGR = $groupsel  -->  IGR = $n_groupsel\n\n"
	    append distext "--------------------------------------------------\n\n"
	}
	if { $inp(IFSO) > 1} {
	    append distext "> NON-STANDARD ORIGIN SHIFT::\n"
	    append distext "IX = $inp(IX),    IY = $inp(IY),    IZ = $inp(IZ)\n"
	    append distext "--------------------------------------------------\n\n"	    
	}
	# verify which unit-cell parameter must be read & read it !!!
	WhichPar2Print  
    } elseif { $species == "slab" } {
	append distext "> LAYER GROUP::\n"
	append distext "IGR = $n_groupsel  -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
	# verify which unit-cell parameter must be read & read it !!!
	WhichPar2Print  
    } elseif { $species == "polymer" } {
	append distext "> ROD GROUP::\n"
	append distext "IGR = $n_groupsel -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
	append distext "> UNIT CELL PARAMETER::\nA:  $inp(A)\n"
	append distext "--------------------------------------------------\n\n"
    } elseif { $species == "molecule" } {
	append distext "> POINT GROUP::\n"
	append distext "IGR = $n_groupsel  -->  $groupsel\n"
	append distext "--------------------------------------------------\n\n"
    }
	
    append distext "> NUMBER OF NON-EQUIVALENT ATOMS::\nNATR = $inp(NATR)\n\n"
    append distext "> ATOMIC NUMBERS & COORDINATES OF NON-EQUIVALENT ATOMS::\n"
    # read Nat,X,Y,Z
    for {set i 1} {$i <= $inp(NATR)} {incr i} {
	AtomNames
	append distext [format "%-3d %-4s %10.5f %10.5f %10.5f\n" \
	$inp(NAT,$i) [Nat2Aname $inp(NAT,$i)] $inp(X,$i) $inp(Y,$i) $inp(Z,$i)]
    }
    append distext "--------------------------------------------------"
    # now display the updated information
        
    # first delete old text if nessecary
    if { [winfo exists $t] } {
	set text [DispText $t $distext 50 25 1]
	$text config -state disabled
    } else {
	set t [StatusWidget]
    }
}


proc StatusWidget {} {
    global distext
    
    # TOPLEVEL WITH TEXT
    set tx [xcToplevel .opftx "Crystal Input:  status" "Crystal Input" \
	    . 330 100]
    AlwaysOnTopON . .opftx

    set text [DispText $tx.frm1 $distext 50 25]
    $text config -state disabled
    
    set frm [frame $tx.f -height 10]
    pack $frm -side bottom -before $tx.frm1 -fill x
    set ok [button $frm.ok -text "Close" -command "destroy $tx"]
    set updat [button $frm.upd -text "Update Status" \
	    -command [list UpdateStatus $tx.frm1]]    
    pack $ok $updat -side left -expand 1 -ipadx 2 -ipady 2 -pady 10
    return $tx
}


proc PreSetOrigin {} {
    global n_groupsel species inp

    # Origin settings is just for crystals & rhombohedral one;
    # if user has changed $species or group --> if origin settings become 
    # meaningless, make a note to user
    if { $species != "crystal" } {
	set b [tk_dialog .preorgset WARNING "Species has been changed \
		and is not any more a CRYSTAL, so Origin Setting is \
		meaningless !!" warning 0 OK]
	return
    }
    
    set_origin        
}
    

proc PreGeom_sym_input {w t} {
    global species

    set b 1
    if { $species == "molecule" } {
	set b [tk_dialog .pregeom WARNING "Species has been changed to \
		Molecule & for Molecules there is no Cell Parameters!!" \
		warning 0 OK]
    }
    if { $b == 1 } {
	geom_sym_input .opflgeom $t
    }
    return
}


proc DestroyOpfl {t1 t2} {
    
    # delete trace on groupsel variable
    xcDeleteState c95
    xcDeleteState openinput
    xcUpdateState
    xcTraceDelete groupsel
    CancelProc $t1
    CancelProc $t2
    return
}
