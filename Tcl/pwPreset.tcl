#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/pwPreset.tcl                                 #
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc pwInputPreset {file} {
    global system pw
    
    # just for any case
    cd $system(SCRDIR)

    if { [file exists $system(SCRDIR)/nuclei.charges] } {
	file delete $system(SCRDIR)/nuclei.charges
    }

    #------------------------------------------------------------------------
    # first scan the "file" and found the "ntyp" (i.e. number
    # of different atomic species)
    #------------------------------------------------------------------------
    set f_content [ReadFile $file]
    set pw(ntyp) 0
    set is_new_input 0
    set is_old_input 0
    foreach line [split $f_content \n] {
	foreach field [split $line ,] {
	    if [string match -nocase *ntyp* $field] {
		set pw(ntyp) [lindex [split $field =] 1]
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
    if { $pw(ntyp) == 0 || (!$is_new_input && !$is_old_input) } {
	# the file is not PW-input file
	ErrorDialog "file \"[file tail $file]\" is not a PWSCF Input File !!!"
	CloseCase
	return 0
    }

    if { $is_old_input } {
	set pw(input) old
    } elseif { $is_new_input } {
	set pw(input) new
    }
    
    uplevel 1 { pwPreset $file Input }

    return $pw(status)
}


proc pwOutputPreset {file} {
    global system pw

    # just for any case
    cd $system(SCRDIR)

    if { [file exists $system(SCRDIR)/nuclei.charges] } {
	file delete $system(SCRDIR)/nuclei.charges
    }

    #------------------------------------------------------------------------
    # first scan the "file" and found the "number of atomic types" 
    # (i.e. number of different atomic species)
    #------------------------------------------------------------------------
    set pw(ntyp) 0
    set pw(version) {}

    set fid [open $file r]
    while { ! [eof $fid] } {
	gets $fid line
	if { [string match "* Program PWSCF *" $line] } {
	    #
	    # get the PWscf version of the output
	    #
	    set ver [split [string trimleft [lindex $line 2] {v.}] .]
	    set pw(version) [lindex $ver 0].[lindex $ver 1]
	    if { [llength $ver] == 3 } {
		append pw(version) [lindex $ver 2]
	    }
	    xcDebug -stderr "PWSCF version of the output: $pw(version)"
	}
	if [string match {*number of atomic types*} $line] {
	    set pw(ntyp) [lindex $line end]
	    break
	}
    }
    close $fid

    if { $pw(ntyp) == 0 } {
	set pw(status) 0

	# some older PWSCF versions didn't have "number of atomic
	# types" printout -> user will have to make nuclei.charges
	# file by himself/herself !!!

	# todo: "exec xterm" should in future be replaced by some Tcl Console
	if { [catch {exec \
			 xterm -geometry 80x15 \
			 -e $system(TOPDIR)/scripts/pwGetNC.sh output $file} \
		  errMsg] } {
	    # the file is not PW-input file
	    set id [tk_dialog [WidgetName] Error \
			"ERROR: File \"[file tail $file]\" is not a PWSCF Output File !" error 0 OK {Error Info}]
	    if { $id == 1 } {
		set t1 [xcDisplayVarText $errMsg {Error Info}]
		tkwait window $t1
	    }
	    CloseCase
	    return
	}
	set pw(status) 1
    } else {
	uplevel 1 { pwPreset $file Output }
    }

    #
    # Query what user would like to do: -lc|-oc|-a
    # 
    
    if { $pw(status) == 1 } {
	set t [xcToplevel [WidgetName] "Question" "Question"]
	set f [frame $t.f -class StressText]
	set l [label $t.f.l -text "What would you like to do:" \
		   -relief groove -bd 2]
	pack $f $l -side top -expand 1 -fill both \
	    -padx 5 -pady 5 -ipadx 3 -ipady 3

	set pw(output_what) {Display Optimized Coordinates}
	RadioButtons $f pw(output_what) top \
	    {Display Initial Coordinates} \
	    {Display Optimized Coordinates} \
	    {Display Latest Coordinates} \
	    {Display All Coordinates as Animation}
	
	button $f.continue -text "Continue" \
	    -command [list pwOutputPresetWhat $t $file] \
	    -default active
	pack $f.continue -side top -expand 1 -padx 5 -pady 5


	#
	# check which type of coordinates are present in the output
	# file (and enable only the corresponding radiobuttons)
	#
	# radiobuttons = $f.choices."digit":
	# 0 = initial
	# 1 = optimized
	# 2 = latest
	# 3 = anim
	set ctypes "initial"
	set ctypes [xcCatchExecReturn $system(awk) -f $system(AWKDIR)/pwo_coortype.awk $file]
	puts stderr "ctypes = $ctypes"

	if { [catch {set ctypes [exec $system(awk) -f  $system(AWKDIR)/pwo_coortype.awk $file]}] } {
	    ErrorDialog "error while executing \"pwo_coortype.awk\" program"
	    return
	}

	if { ! [string match *init* $ctypes] } {
	    $f.choices.0 configure -state disabled
	} else {
	    set pw(output_what) {Display Initial Coordinates}
	}

	if { ! [string match *inter* $ctypes] } {
	    if { ! [string match *opt* $ctypes] } {
		$f.choices.3 configure -state disabled
	    }
	    $f.choices.2 configure -state disabled
	} else {
	    set pw(output_what) {Display Latest Coordinates}
	}
	
	if { [string match *dyna* $ctypes] } {
	    set pw(output_what) {Display All Coordinates as Animation}
	}

	if { ! [string match *opt* $ctypes] } {
	    $f.choices.1 configure -state disabled
	} else {
	    set pw(output_what) {Display Optimized Coordinates}
	}       
	
	tkwait window $t
    }
    return $pw(status)
}


proc pwOutputPresetWhat {t file} {
    global system pw openExtStruct
    
    if { [info exists pw(output_flag)] } {
	# calling pwOutputPresetWhat from scripting filter
	set flag $pw(output_flag)
    } else {
	if { $pw(output_what) == {Display Initial Coordinates} } {
	    set flag {--inicoor}
	} elseif { $pw(output_what) == {Display Optimized Coordinates} } {
	    set flag {--optcoor}
	} elseif [string match $pw(output_what) {Display Latest Coordinates}] {
	    set flag {--latestcoor}
	} else {
	    set flag {--animxsf}
	}
    }

    cd $system(SCRDIR)

    xcDebug -stderr "pwOutputPresetWhat: $system(TOPDIR)/scripts/pwo2xsf.sh $flag $file"
    
    set openExtStruct(program) [list sh $system(TOPDIR)/scripts/pwo2xsf.sh $flag]
    if { [xcCatchExec sh $system(TOPDIR)/scripts/pwo2xsf.sh $flag $file > pwo2xsf.xsf] } {
	CloseCase
	return
    }
    if { [winfo exists $t] } {
	destroy $t
    }
}


proc pwPreset {file {filetype {Input}}} {
    global system pw xcMisc
    
    #------------------------------------------------------------------------
    # now ask user for ityp->nat replacement
    #------------------------------------------------------------------------
    
    set grab [grab current]
    catch { grab release $grab }
    
    set t [xcToplevel [WidgetName] \
	       "PWscf $filetype: \"[file tail $file]\"" "PWSCF $filetype" . 0 0 1]
    
    set pw(disptext) [xcDisplayFileText $file \
			  "PWscf $filetype File: \"[file tail $file]\"" . 200 200 1]
    
    set f1 [frame $t.1]
    set f2 [frame $t.2]
    set f3 [frame $t.3]
    pack $f1 $f2 $f3 -side top -padx 10 -pady 10 -fill x
    
    if { ! [info exists pw(input)] } {
	set pw(input) none
    }
    if { ! [info exists pw(version)] } {
	set pw(version) 0.0
    }
    if { $pw(input) != "new" && $pw(version) < 1.2 } {
	frame $f1.title -class StressText
	label $f1.title.l -text "Specify Atomic Numbers of Species" \
	    -relief groove -bd 2
	$f1.title.l config -font [ModifyFontSize $f1.title.l 18]
	pack $f1.title -side top -expand 1
	pack $f1.title.l -ipadx 10 -ipady 3 -padx 10 -pady 10 -expand 1
	
	for {set i 0} {$i < $pw(ntyp)} {incr i} {
	    set ii [expr $i + 1]
	    set pw(ityp,$ii) $ii
	}
	set labellist  [list "Atom type:" "Atomic Number:"]    
	set arraylist  [list ityp nat]
	set arraytype  [list posint posint]
	set buttonlist [list 1 \
			    [list {Periodic Table} scroll_ptableSelectNAT $t pw nat]]	    
	
	ScrollEntries \
	    $f1 \
	    $pw(ntyp) \
	    "Set Atomic Number for Species #" \
	    $labellist \
	    $arraylist \
	    $arraytype \
	    15 \
	    pw \
	    $buttonlist \
	    3    
	update
	set width [winfo width $f1]
    } else {
	set width 400
    }

    # old + new input

    set mf [frame $f2.f -class StressText -relief groove -bd 2]
    set m  [message $mf.m -justify left -anchor e \
		-width [expr $width - 25] -text {Multi-slab and Molecule-in-a-Box Note:

If your structure is a so-called multi-slab model of the surface, then a nicer display is yielded if only a single slab is rendered. If you want to render multi-slab as a single slab, then select the "reduce dimension to 2D" radiobutton below.

In case you have molecule-in-a-box then select the "reduce dimension to 0D" radiobutton below.}]
    
    set xcMisc(reduce_to) 3
    set rb3 [radiobutton $mf.rb3 -text "do not reduce dimensionality" \
		 -variable xcMisc(reduce_to) -value 3]
    set rb2 [radiobutton $mf.rb2 -text "reduce dimension to 2D" \
		 -variable xcMisc(reduce_to) -value 2]
    set rb1 [radiobutton $mf.rb1 -text "reduce dimension to 1D" \
		 -variable xcMisc(reduce_to) -value 1]
    set rb0 [radiobutton $mf.rb0 -text "reduce dimension to 0D" \
		 -variable xcMisc(reduce_to) -value 0]
    pack $mf $m $rb3 $rb2 $rb1 $rb0 -side top -expand 1 -fill both \
	-padx 5 -pady 5 -ipadx 3 -ipady 3
    
    # OK and cancel button should be cerated as well
    set ok  [button $f3.ok  -text OK     -command [list pwPresetOK $t] \
		 -default active]    
    set can [button $f3.can -text Cancel -command [list pwPresetCan $t]]
    pack $can $ok -side left -expand 1 -padx 10
    update
    if { [winfo exists $pw(disptext)] } {
	raise $pw(disptext) .
	raise $t $pw(disptext)
    }
    
    tkwait window $t

    if { $grab != "" } {
	catch { grab $grab }
    }
}

    
#
# When OK button is preset, then "nuclei.charges" file should be written
# to $system(SCRDIR) directory
#
proc pwPresetOK t {
    global varlist foclist pw system
    
    if { $pw(input) != "new" && $pw(version) < 1.2 } {
	# PWscf < 1.2
	if { ![check_var $varlist $foclist] } {
	    return 0
	}
	
	# its seems OK; write the "nuclei.charges" file and destroy 
	# the toplevel window
	
	set out "$pw(ntyp)\n"
	for {set i 1} {$i <= $pw(ntyp)} {incr i} {
	    append out "$pw(ityp,$i) $pw(nat,$i)\n"
	}
    } else {
	# PWscf >= 1.2: we need 0 ityp->nat conversion
	set out 0
    }
    
    # always write the nuclei.charges file (some workaround for PWscf
    # version 1.2)
    evalInScratch {
	WriteFile nuclei.charges $out w
    }

    destroy $t
    if { [winfo exists $pw(disptext)] } { destroy $pw(disptext) }
    set pw(status) 1
}
	    

proc pwPresetCan t {
    global pw system

    ####################
    # just in any case #
    cd $system(SCRDIR)
    ####################
    if [file exists nuclei.charges] {
	file delete nuclei.charges
    }
    destroy $t
    if { [winfo exists $pw(disptext)] } { destroy $pw(disptext) }
    set pw(status) 0
}
