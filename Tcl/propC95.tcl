#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/propC95.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc PropC95Cmd {cmd} {
    global prop system grafdata periodic err
    
    xcDebug "\nIn PropC95Cmd; prop(newk) == $prop(newk);   cmd == $cmd\n"
    
    # there are some properties, that are caracteristic for periodic systems
    set perlist [list BWID DOSS BAND]
    foreach com $perlist {
	if { $cmd == $com && $periodic(dim) == 0} {
	    #return silently
	    return
	}
    }

    # delete old unit-25
    if { [file exists $prop(dir)/$prop(file)25] } {
	file delete $prop(dir)/$prop(file)25
    }

    # INFO
    set info 0
    if { $cmd == "INFO" } {
	set input "BASE\n0\n"
	if { $periodic(dim) > 0 } {
	    set button [tk_dialog [WidgetName] Question \
		    "Do You want \"band widths\" to be included in INFO record" question 0 Yes No]
	    if { $button == 0 } {
		set info 1
		if ![RunC95 $system(c95_properties) err "END\n"] {
		    return
		}
		append input [NEWK_Init]
		GetC95Info BWID $system(SCRDIR)/xc_output.$system(PID)
		append input "BWID\n1 $prop(n_band)\n"
	    }
	}
	append input "END\n"
	if ![RunC95 $system(c95_properties) err $input] {
	    return
	}
	DispC95Output $system(SCRDIR) {} "INFO: $prop(unit9)" 1
	return
    }

    # 
    # do user want to run PropC95Cmd that require NEWK ????
    # NOTE: band do not require NEWK, but if we want correct Fermi Energy
    #       to be displayed on the spaghetti-graphs, this is needed !!!
    #
    set newk_cmd_list [list BWID DOSS BAND] 
    foreach newkcmd $newk_cmd_list {
    	if { $cmd == $newkcmd } {
	    set prop(newk) 1
	    set prop(newk_script) [NEWK_Init]	
	    if { $prop(newk_script) == "" } {
		set prop(newk) 0
		return
	    }
	    break
	}
    }

    if { $cmd == "BWID" } {		
	GetC95Info BWID $system(SCRDIR)/xc_output.$system(PID)

	# SINTAX selband will be {{alfa band} {betaband}}
	# SelBandIntv $prop(n_band)
	set input "$prop(newk_script)BWID\n1 $prop(n_band)\nEND\n"

	# run properties program
	if { ![RunC95 $system(c95_properties) {message "CRYSTAL program is calculating band widths.\nIt can take some time, so PLEASE WAIT !!!"} $input \
		   {} "BAND LIMITS" $system(SCRDIR)] } {
	    # error occure
	    #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	    return
	}

	# update crystal95's output
	DispC95Output $system(SCRDIR) {} {} 1
	
	# now load a Graph
	set id [NextGrapherID]
	set grafdata($id,Xtitle) "Band Widths: $prop(unit9)"
	set grafdata($id,Xicon) "Band Widths"
	BWIDGraph $system(SCRDIR)/xc_output.$system(PID)
	Grapher BARGraph
	set cmd ""
	return
    }

    if { $cmd == "DOSS" } {
	global doneDOSS
	if ![info exists prop(NPRO)] { set prop(NPRO) 0 }
	GetC95Info DOSS $system(SCRDIR)/xc_output.$system(PID)

	DOSS_Init
	if { $doneDOSS == 0 } {
	    return
	}
	
	xcDebug "after DOSS_Init"
	# T.K.: now I must set prop(N,$i) for each $i
	for {set i 1} {$i <= $prop(NPRO)} {incr i} {
	    set prop(N,$i) [llength $prop(NDM,$i)]
	}
	set    input "$prop(newk_script)"
	append input "DOSS\n"
	append input "$prop(NPRO) $prop(NPT) $prop(firstband) \
		$prop(lastband) 1 $prop(NPOL) $prop(NPR)\n"
	if { $prop(firstband) < 0 && $prop(lastband) < 0 } {
	    append input "$prop(BMI) $prop(BMA)\n"
	}
	for {set i 1} {$i <= $prop(NPRO)} {incr i} {
	    if { $prop(proj_NPRO,$i) != "set of N AOs" && $prop(N,$i) > 0} {
		set prop(N,$i) -$prop(N,$i)
	    }
	    append input "$prop(N,$i) $prop(NDM,$i)\n"
	}
	append input "END\n"
	# run properties program
	if ![RunC95 $system(c95_properties) err $input {} {} $system(SCRDIR)] {
	    # error occure
	    #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	    return
	}
	xcDebug "BEFORE DOSSGraph"

	# update crystal95's output
	DispC95Output $system(SCRDIR)

	# now load a Graph
	set id [NextGrapherID]
	set grafdata($id,Xtitle) "Density of States: $prop(unit9)"
	set grafdata($id,Xicon) "Density of States"
	DOSSGraph $prop(NPRO)
	Grapher XYGraph
	return
    }

    if { $cmd == "BAND" } {
	global BzOK
	if { ! [info exists prop(NPRO)] } { 
	    set prop(NPRO) 0 
	}
	GetC95Info DOSS $system(SCRDIR)/xc_output.$system(PID)
	
	if { [Bz_MakeToplevel] == 0 } {
	    return
	}
	#if { $BzOK(done) == 0 } {
	#    # cancel button was pressed; RETURN
	#    return
	#}
	append prop(c95_BAND_script) "END\n"
	xcDebug "BAND_Init: 1"
	if { ![RunC95 $system(c95_properties) err $prop(c95_BAND_script) \
		   {} {} $system(SCRDIR)] } {
	    # error occure
	    #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	    return
	}
	xcDebug "BEFORE BANDGraph"

	# update crystal95's output
	DispC95Output $system(SCRDIR)

	# now load a Graph
	set id [NextGrapherID]
	set grafdata($id,Xtitle) "Band Structure: $prop(unit9)"
	set grafdata($id,Xicon)  "Band Structure"
	BANDGraph $prop(NLINE)
	Grapher XYGraph    
	return
    }
}


##############################################################################
# NEWK
# proc return NEWK-block for runing "properties" program
proc NEWK_Init {} {    
    global periodic done newk prop system

    if { ! [info exists prop(IS)]  }  { set prop(IS) 0 }
    if { ! [info exists prop(ISHF)] } { set prop(ISHF) 0 }
    if { ! [info exists prop(ISP)]  } { set prop(ISP) 0 }

    if { $periodic(dim) > 0 } {
	set t [xcToplevel .bwid "NEWK Options" "NEWK" . 100 50 1]
	catch { grab $t }
	set f1 [frame $t.f1 -relief raised -bd 2]
	pack $f1 -side top -fill both -ipady 10
	frame $f1.f1 -relief flat
	frame $f1.f2 -relief flat

	if ![info exists newk(option)] {
	    set newk(option) "Same as in SCF"
	}
	if ![info exists newk(IFE)] {
	    set newk(IFE) "Yes"
	}
	pack $f1.f1 $f1.f2 -side top -padx 10 -fill both -expand 1
	
	if { $system(c95_version) != "06" && $system(c95_version) != "09" && $system(c95_version) != "14" } {
	    #
	    # CRYSTAL prior to version-06
	    #
	    RadioButVarCmd $f1.f1 "Which (IS,ISHF,ISP) values to take:" \
		newk(option) NEWK_RadioButCmd top left 1 1 \
		"Same as in SCF" "New (IS,ISHF,ISP) values" 
	} else {
	    #
	    # CRYSTAL06 or later
	    #
	    RadioButVarCmd $f1.f1 "Which (IS,ISP) values to take:" \
		newk(option) NEWK_RadioButCmd top left 1 1 \
		"Same as in SCF" "New (IS,ISP) values" 
	}

	RadioBut $f1.f2 "Calculate Fermi Energy:" \
	    newk(IFE) left left 1 1 \
	    "Yes" "No"
	
	message $f1.msg -justify center -aspect 400 -text \
		"\nWARNING: In order to obtain correct Fermi energy and \
		eigenvalues spectra when a shift of eigenvalues was \
		requested in SCF (LEVSHIFT, SPINLOCK)\n it is required to\
		re-calculate Fermi energy" 
	pack $f1.msg -side bottom -padx 10 -pady 10

	set f2 [frame $t.f2 -relief raised -bd 2]	
	set ok [button $f2.ok -text "OK" -command [list NEWK_OK $t]]
	set can [button $f2.can -text "Cancel" -command [list NEWK_Can $t]]
	pack $f2 -side bottom -fill x 
	pack $ok $can -side left -expand 1 -padx 10 -pady 10
    }

    tkwait variable done    
    # if cancel button was pressed, return 0
    if $newk(cancel) { return }

    if { $newk(IFE) == "Yes" } {
	set prop(IFE) 1
    } else {
	set prop(IFE) 0
    }
    set prop(NPR) 0

    #
    # NEWK
    # is (ishf) isp
    # if (is==0 && dim>0) is1 is2 is3
    # ife npr
    #

    set    input "NEWK\n"
    if { $system(c95_version) != "06" && $system(c95_version) != "09" && $system(c95_version) != "14" } {
	#
	# CRYSTAL prior to version-06
	#
	append input "$prop(IS) $prop(ISHF) $prop(ISP)\n"
    } else {
	#
	# CRYSTAL06 or later
	#
	append input "$prop(IS) $prop(ISP)\n"
    }	
    if { $periodic(dim) > 0 } {
	# add IS1 IS2 IS3 support
	#append input "$prop(IS1) $prop(IS2) $prop(IS3)\n"
    }
    append input "$prop(IFE) $prop(NPR)\n"
    
    return $input
}


proc NEWK_RadioButCmd {what} {
    global prop

    if { $what == "Same as in SCF" } {
	set prop(IS)   0
	set prop(ISHF) 0
	set prop(ISP)  0
    } else {
	NEWK_IS_ISHF_ISP
    }
}


proc NEWK_IS_ISHF_ISP {} {
    global prop ok system periodic

    set oldgrab [grab current]
    set t [xcToplevel .newkis "Enter values" "NEWK" . 120 70 1]    
    set f1  [frame $t.f1 -relief raised -bd 2]
    set f11 [frame $t.f11 -relief raised -bd 2]
    set f2  [frame $t.f2 -relief raised -bd 2]

    if { $system(c95_version) != "06" && $system(c95_version) != "09" && $system(c95_version) != "14" } {
	#
	# CRYSTAL prior to version-06
	#
	Entries $f1 {IS: ISHF: ISP:} {prop(IS) prop(ISHF) prop(ISP)} 3
	set varlist [list {prop(IS) posint} {prop(ISHF) posint} {prop(ISP) posint}]
	set foclist [list $f1.frame.entry1 $f1.frame.entry2 $f1.frame.entry3]
    } else {
	#
	# CRYSTAL06 or later
	#
	Entries $f1 {IS: ISP:} {prop(IS) prop(ISP)} 2
	set varlist [list {prop(IS) posint} {prop(ISP) posint}]
	set foclist [list $f1.frame.entry1 $f1.frame.entry2]
    }
    
    # if { $periodic(dim) > 0 } {
    # 	# add IS1 IS2 IS3 support
    # 	Entries $f11 {IS1 IS2 IS3} {prop(IS1) prop(IS2) prop(IS3)} 3
    # 	append varlist " {prop(IS1) posint} {prop(IS2) posint} {prop(IS3) posint}"
    # 	append foclist " $f11.frame.entry1 $f11.frame.entry2 $f11.frame.entry3"
    # }
    focus $f1.frame.entry1

    set ok  [button $f2.ok -text OK \
	    -command [list check_var $varlist $foclist]]
    set can [button $f2.can -text Cancel -command NEWK_IS_Can]

    pack $f1 $f11 $f2 -side top -fill both
    pack $ok $can -side left -expand 1 -padx 10 -pady 10

    update
    catch { grab $t }

    tkwait variable ok
    if { [winfo exists $t] } { 
	catch { grab release $t }
	destroy $t 
    }
    if { $oldgrab != "" } {
	catch { grab $oldgrab }
    }
}


proc NEWK_IS_Can {} {
    global prop ok newk
    
    set prop(IS)     0
    set prop(ISHF)   0
    set prop(ISP)    0
    set newk(option) "Same as in SCF"
    set ok           1
}


proc NEWK_Can {t} {
    global done newk

    set done 1
    set newk(cancel) 1
    if { [winfo exists $t] } { 
	catch { grab release $t }
	destroy $t 
    }
}


proc NEWK_OK {t} {
    global done newk

    puts stdout "NEWK_OK"
    set done 1
    set newk(cancel) 0
    if { [winfo exists $t] } { 
	catch { grab release $t }
	destroy $t 
    }
}


proc SelBandIntv w {
    global prop done
    
    set l [label $w.l \
	    -text "Number of Bands:   $prop(n_band)" \
	    -relief flat -anchor w]
    set font [$l cget -font]
    set font [ModifyFont $font $l -underline 1 -weight bold]
    $l config -font $font
    set l1 [label $w.l1 -text \
	    "Please choose the range of bands to consider:" \
	    -relief flat -anchor w]
    pack $l $l1 -side top -padx 5 -fill x -expand 1
    set foclist [OneEntries $w {"First band to consider:" \
	    "Last band to consider:"} \
	    {prop(firstband) prop(lastband)} 37 8 5 -fill x -expand 1]
    return $foclist
}


proc SelEnerIntv {w} {
    global prop ok

    # in future I will have to make some mechanism to check in energy 
    # interval is correctly entered
    xcDebug "inSelEnerIntv"

    set l1 [label $w.l1 -text \
	    "Please enter energy interval for DOSS calculation:\n(Boundaries of energy interval must be in a band gap !!!)" \
	    -relief flat -justify left -anchor w]

    pack $l1 -side top -pady 5 -fill x -expand 1
    OneEntries $w {"Minimum Energy:" "Maximum Energy:"} \
	    {prop(BMI) prop(BMA)} 37 8
}


#proc checkIntv {varlist foclist lastvar last focus} {
#    global ok done
#    upvar #0 $lastvar var
#
#    puts stdout "lastvar = $lastvar"
#    check_var $varlist $foclist
# 
#    if [info exist ok] {
#	if $ok {
#	    if { $last >= $var } {
#		set done 1
#	    } else {
#		tk_dialog .num ERROR "ERROR !\nYou have specified to large \
#			number for \"$lastvar\" variable. It should be \
#			lower or equal to $last" error 0 OK
#		set focus $focus
#	    }
#	}
#    }
#}


##############################################################################
# DENSITY OF STATES
proc DOSS_Init {} {
    global prop doneDOSS

    set doneDOSS 0

    set t [xcToplevel [WidgetName] "Density od States" "DOSS" . 100 50 1]
    catch { grab $t }
    set f [frame $t.f -relief raised -bd 2]
    set f2 [frame $t.f2 -relief raised -bd 2]    
    set f3 [frame $t.f3 -relief raised -bd 2]
    set f4 [frame $t.f4 -relief raised -bd 2]
    set f5 [frame $t.f5 -relief raised -bd 2]
    pack $f $f2 $f3 $f4 $f5 -expand 1 -fill both -side top   

    set ff1 [frame $f.1 -relief groove -borderwidth 2]
    pack $ff1 -side top -padx 2 -pady 2 -fill x -expand 1
    set sc [scale $ff1.scale -from 0 -to 10 -length 100 -variable prop(NPRO) \
	    -orient horizontal -label "Number of DOS projections:" \
	    -tickinterval 2 \
	    -digits 2 -resolution 1 -showvalue true -width 10]
    pack $sc -side top -fill x -padx 5

    if ![info exists prop(doss_criteria)] {
	set prop(doss_criteria) "band-interval criteria"
    } 

    set ff2 [frame $f.2 -relief flat]
    RadioBut $ff2 "Select criteria for spanning DOS:" \
	    prop(doss_criteria) top left 1 1 \
	    "band-interval criteria" "energy-interval criteria"     
    trace variable prop(doss_criteria) w xcTrace
    ##########################################################################
    OneEntries $f2 {"Number of uniformly spaced energy points:" \
	    "Number of Legendre polynomials:"} \
	    {prop(NPT) prop(NPOL)} 37 8
    focus $f2.frame1.entry1
    set varlist [list {prop(NPT) posint} {prop(NPOL) posint}]
    set foclist [list $f2.frame1.entry1 $f2.frame2.frame1]
    pack $ff2 -side top -expand 1 -fill x 
    ###########################################################################
    # BAND INTERVAL
    SelBandIntv $f3
    set prop(bandIntv_labels)  \
	    [list $f3.l $f3.l1 $f3.frame1.lab1 $f3.frame2.lab2]
    set prop(bandIntv_entries) [list $f3.frame1.entry1 $f3.frame2.entry2]

    append varlist " {prop(firstband) posint} {prop(lastband) posint}"
    append foclist " $f3.frame1.entry1 $f3.frame2.entry2"
	
    ###########################################################################
    # ENERGY INTERVAL
    SelEnerIntv $f4
    set prop(enerIntv_labels)  [list $f4.l1 $f4.frame1.lab1 $f4.frame2.lab2]
    set prop(enerIntv_entries) [list $f4.frame1.entry1 $f4.frame2.entry2]

    append varlist " {prop(BMI) real} {prop(BMA) real}"
    append foclist " $f4.frame1.entry1 $f4.frame2.entry2"

    proc DOSS_Init_BandIntv {} {
	global prop
	set dis_c [GetWidgetConfig button -disabledforeground]
	set ena_c [GetWidgetConfig button -foreground]
	# disable "ENERGY-INTERVAL"
	foreach lab $prop(enerIntv_labels) {
	    $lab config -fg $dis_c
	}
	foreach entry $prop(enerIntv_entries) {
	    $entry config -relief flat -state disabled
	}
	# enable "BAND-INTERVAL"
	foreach lab $prop(bandIntv_labels) {
	    $lab config -fg $ena_c
	}
	foreach entry $prop(bandIntv_entries) {
	    $entry config -relief sunken -state normal
	}
    }
    
    proc DOSS_Init_EnerIntv {} {
	global prop
	set dis_c [GetWidgetConfig button -disabledforeground]
	set ena_c [GetWidgetConfig button -foreground]
	# enable "ENERGY-INTERVAL"
	foreach lab $prop(enerIntv_labels) {
	    $lab config -fg $ena_c
	}
	foreach entry $prop(enerIntv_entries) {
	    $entry config -relief sunken -state normal
	}
	foreach lab $prop(bandIntv_labels) {
	    $lab config -fg $dis_c
	}
	foreach entry $prop(bandIntv_entries) {
	    $entry config -relief flat -state disabled
	}
    }
    #
    # query the prop(doss_criteria) state
    #
    trace variable prop(doss_criteria) w xcTrace
    xcTrace prop doss_criteria w
        
    set ok [button $f5.ok -text "OK" \
	    -command [list DOSS_Init2 $t $varlist $foclist]]
    set can [button $f5.can -text "Cancel" -command [list DOSS_InitCan $t]]
    pack $can $ok -side left -expand 1 -padx 10 -pady 10
    
    tkwait variable doneDOSS
    # delete prop(doss_criteria)'s trace
    xcTraceDelete prop(doss_criteria)
    
    xcDebug "end of DOSS_Init"
    return 1
}
    

proc DOSS_InitCan tplw {
    global doneDOSS
    xcTraceDelete prop(doss_criteria)
    CancelProc $tplw doneDOSS 
}

proc DOSS_Init2 {t varlist foclist} {
    global prop doneDOSS

    if { $prop(doss_criteria) == "band-interval criteria" } {
	set varl [lrange $varlist 0 3]
    } else {
	set varl [lrange $varlist 0 1]
	append varl " [lrange $varlist 4 5]"
	set prop(firstband) -1
	set prop(lastband)  -1
    }
    if ![check_var $varl $foclist] {
	return
    }
    # check if band interval is specified correctly:
    if { $prop(doss_criteria) == "band-interval criteria" } {
	if { $prop(firstband) < 1 } {
	    tk_dialog [WidgetName] ERROR \
		    "ERROR: Lowest boundary of band-interval \"$prop(firstband)\" is out of range, should be greater than 0. Try Again !!!" \
		    error 0 OK
	    focus [lindex $foclist 2]
	    return
	}
	if { $prop(lastband) > $prop(n_band) } {
	    tk_dialog [WidgetName] ERROR \
		    "ERROR: Upper boundary of band-interval \"$prop(firstband)\" is out of range, should be lower or equal than $prop(n_band). Try Again !!!" error 0 OK
	    focus [lindex $foclist 3]
	    return
	}
	if { $prop(firstband) > $prop(lastband) } {
	    tk_dialog [WidgetName] ERROR \
		    "ERROR: Lowest boundary of band-interval \"$prop(firstband)\" greater than upper boundary \"$prop(lastband)\". Try Again !!!" \
		    error 0 OK
	    focus [lindex $foclist 2]
	    return
	}
    }

    xcDebug "DOSS_Init2:: $prop(NPRO)    $prop(doss_criteria)"
    if { [winfo exists $t] } {
	catch { grab release $t }
	destroy $t
    }

    if { $prop(NPRO) > 0 } {
	global donePDOSS
	PDOSS_Init
	set doneDOSS $donePDOSS
    } else {
	set doneDOSS 1
    }
}
    
proc PDOSS_Init {} {
    global prop donePDOSS

    set t [xcToplevel .npro "Projected DOSS" "PDOSS" . 100 50 1]
    catch { grab $t }
    set fb [frame $t.fb]
    pack $fb -side bottom -expand true -fill both 
    # and one frame where canvas&scrollbar will be!!
    set ft [frame $t.ft -relief sunken -bd 2]
    pack $ft -side top -expand true -fill both 
	
    set c [canvas $ft.canv -yscrollcommand [list $ft.yscroll set]]
    set scb [scrollbar $ft.yscroll -orient vertical -command [list $c yview]]
    pack $scb -side right -fill y
    pack $c -side left -fill both -expand true
	
    # create FRAME to hold every LABEL&ENTRY
    set f [frame $c.f -bd 0]
    $c create window 0 0 -anchor nw -window $f
    set varlist ""
    set foclist ""
    for {set i 1} {$i <= $prop(NPRO)} {incr i 1} {	    
	set gro [frame $f.fr$i -relief groove -bd 2]
	pack $f.fr$i -padx 5 -pady 5
	if ![info exists prop(proj_NPRO,$i)] {
	    set prop(proj_NPRO,$i) "set of N AOs"
	}
	RadioBut $gro "Density of states projected onto:" \
		prop(proj_NPRO,$i) top left \
		0 1 "set of N AOs" "set of all AOs of the N atoms" 
	set e [Entries $gro \
		{"Specify sequence numbers of AOs/Atoms for projection:"} \
		prop(NDM,$i) 40 1 top -fill x]
	button $gro.b -text "Select sequence" \
		-command [list SelectSequence $i prop(NDM,$i)]
	pack $gro.b -side left -padx 10 -pady 10
	lappend varlist [list prop(NDM,$i) posint]
	append foclist " $e"
    }

    focus $f.fr1.frame.entry1

    puts stdout "FOCLIST: $foclist\n\n"
    puts stdout "VARLIST: $varlist"
    set child [lindex [pack slaves $f] 0]
    
    # set the focus to first entry that upper FOR-LOOP create
    
    tkwait visibility $child
    set width [winfo width $f]
    set height [winfo height $f]
    if { $prop(NPRO) < 3 } {
	$c config -width $width -height $height 
    } else {
	$c config -width $width -height \
		[expr $height / $prop(NPRO) * 3] \
		-scrollregion "0 0 $width $height"
    }
    
    button $fb.butok -text OK -command \
	    [list PDOSS_InitOK $t $varlist $foclist]
    
    button $fb.butcan -text Cancel -command \
	    [list CancelProc $t donePDOSS]
    
    pack $fb.butok $fb.butcan -side left \
	    -expand 1 -padx 10 -pady 10
    
    tkwait variable donePDOSS
    
    catch { grab release $t }
    destroy $t
}


proc PDOSS_InitOK {t varlist foclist} {
    global donePDOSS

    # varlist has the following shape {varname posint} {varname posint}
    set n 0
    foreach var $varlist {
	xcDebug "PDOSS_InitOK:: var=$var"
	set varn [lindex $var 0]
	upvar #0 $varn val 
	set type [lindex $var 1]
	xcDebug "PDOSS_InitOK:: varn=$varn ; val=$val ; type=$type"
	if { $val == "" } {
	    tk_dialog .number2 ERROR "ERROR !\nYou forget to specify \
		    the \"$varn\" variable. Please do so !" error 0 OK
	    focus [lindex $foclist $n]
	    return
	}
	foreach value $val {
	    if { [catch {expr abs($value)}] } {
		# this CATCH specify if $var is a number;
		# if we get 1 --> not number, else number
		# string is not a number
		dialog .number1 ERROR "ERROR !\nYou have specified a character instead of number for \"$varn\" variable.\
			TRY AGAIN \!" error 0 OK
		focus [lindex $foclist $n]
		return
	    }		
	    if { $value != int($value) || [string match *.* $value] } {
		tk_dialog .number2 ERROR "ERROR !\nYou have specified \
			a non-integer number instead of \
			positive integer number \
			for \"$varn\" variable. TRY AGAIN \!" error 0 OK
		focus [lindex $foclist $n]
		return
	    }
	    if { $value < 0 } {
		tk_dialog .number2 ERROR "ERROR !\nYou have specified \
			a negative integer number instead of positive \
			integer number \
			for \"$varn\" variable. TRY AGAIN \!" error 0 OK
		focus [lindex $foclist $n]
		return
	    }
	}
	incr n
    }
    set donePDOSS 1
}
		

proc SelectSequence {i var} {
    global prop

    set what $prop(proj_NPRO,$i)
    xcDebug "SelectSequence $what $var"
    if { $what == "set of N AOs" } {
	set $var [SelectItems "Select the AOs" $prop(n_band)]
    } elseif { $what == "set of all AOs of the N atoms" } {
	set $var [SelectItems "Select Atoms" $prop(n_atom)]
    }
}


proc SelectItems {title num} {
    global prop done value

    set oldgrab [grab current]
    # if .sband already exists, return silenlty
    if [winfo exists .sband] { return }
    set t [xcToplevel .sband $title $title . 100 50 1]
    catch { grab $t }
    # in the future I will have to know whether I'm dealing with 
    # close or open shell
    set bandlist ""
    for {set i 1} {$i <= $num} {incr i} {
	append bandlist "   $i"
    }

    if { $num > 15 } {
	set h 15
    } else {
	set h $num
    }

    set title "$title\n(to select multiple items:\nCTRL-key + mouse-click)"
    label $t.l -text $title -relief flat
    pack $t.l -side top
    set scr [ScrolledListbox2 $t.f -width 10 -height $h -setgrid true \
	    -selectmode extended]
    eval {$scr insert 0} $bandlist
    button $t.ok -text OK -command [list SelectItemsOK $scr]
    pack $t.ok -side left -padx 10 -pady 10
    
    tkwait variable done
    destroy $t
    if { $oldgrab != "" } { 
	catch { grab $oldgrab }
    }
    return $value
}


proc SelectItemsOK {scr} {
    global done value
    
    set done 1
    set value ""
    set vals [$scr curselection]
    foreach va $vals {
	set v [expr $va + 1]
	append value "$v "
    }
}


proc DispC95Output {output_dir {file {}} {title {}} {newWin 0}} {
    global system prop dispC95out unmapWin

    if ![info exists dispC95out(offset)] {
	set dispC95out(offset) 0
    } else {
	incr dispC95out(offset) 15
	if { $dispC95out(offset) > 200 } {
	   set dispC95out(offset) 0
	} 
    }
    set x [expr -50 + $dispC95out(offset)]
    set y [expr  70 + $dispC95out(offset)] 
    # newWin  ...... if newWin == 1 --> always make a new Toplevel
    if { $file == {} } { 
	set file xc_output.$system(PID)
    }
    set fileID [open $output_dir/$file r]
    set output [read $fileID]
    xcDebug "CRYSTAL OUTPUT output::\n\n$output"
    # if .dispC95 already exists, we will just update the display
    if { $title == {} } {
	set title "Crystal95: Output"
    }

    set update 1
    if $newWin {	
	set dispC95out(tplw) [xcToplevel [WidgetName] $title $title . $x $y 1]
	set update 0
	tkwait visibility $dispC95out(tplw)
	set t $dispC95out(tplw)
	xcRegisterUnmapWindow $t $unmapWin(frame,main) \
		C95output_$t -image unmap-C95output
	bind $t <Unmap> [list xcUnmapWindow unmap %W $t \
		$unmapWin(frame,main) C95output_$t]
	bind $t <Map>   [list xcUnmapWindow map %W $t \
		$unmapWin(frame,main) C95output_$t]
    } elseif ![winfo exists .dispC95] {	
	set dispC95out(tplw) [xcToplevel .dispC95 $title $title . $x $y 1]
	set update 0
	tkwait visibility $dispC95out(tplw)
	set t $dispC95out(tplw)
	xcRegisterUnmapWindow $t $unmapWin(frame,main) \
		C95output_$t -image unmap-C95output
	bind $t <Unmap> [list xcUnmapWindow unmap %W $t \
		$unmapWin(frame,main) C95output_$t]
	bind $t <Map>   [list xcUnmapWindow map %W $t \
		$unmapWin(frame,main) C95output_$t]
    } else {
	# wm title $t $title
	;
    }

    if !$update {
	set f2 [frame $dispC95out(tplw).f2 -height 10]
	pack $f2 -side bottom -fill x
	set close [button $f2.cl -text "Close" \
		 -command [list destroy $dispC95out(tplw)]]
	set hid [button $f2.hid -text "Hide" \
		 -command [list DispC95Output_Hide $t]]
	
	pack $hid $close -side left -expand 1 -ipadx 2 -ipady 2 -pady 10
	tkwait visibility $hid
    }
    set text [DispText $dispC95out(tplw).f1 $output 80 35 1]
    $text configure -state disabled

    close $fileID
    return $dispC95out(tplw) 
}


proc DispC95Output_Hide t {
    global unmapWin
    wm withdraw $t
    xcUnmapWindow unmap $t $t $unmapWin(frame,main) C95output_$t
}
    
#set system(PWD) [pwd]
#set system(PID) [pid]
#set periodic(dim) 3
#set prop(n_band)  20
#set prop(n_atom)  4
#lappend auto_path /home/tone/prog/XCrys/Mesa
#DOSS_Init

    
