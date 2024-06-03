#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnFSInit.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnFSInit filedir {
    global wn xcMisc
        
    XCRYSDEN_Logo $filedir

    #####################################
    set wn(dir)      $filedir
    set wn(filehead) [file tail $filedir]
    set head         $wn(filehead)

    wnReadStruct $wn(filehead)
    ######################################

    if [winfo exists .title] { 
	destroy .title 
    }
    set t [xcToplevel .fs_init \
	    "*** XCrySDen: Fermi Surface Preparation: $head" \
	    "Fermi Surface" . 290 000 1]

    set cf .fs_init.f
    if [winfo exists $cf] {
	destroy $cf
    }
    frame $cf
    pack $cf

    #
    # WIEN flags 
    #    
    set wn(c) {}
    if $wn(complex) { 
	set wn(c) c 
    }

    #
    # check if *.klist exists and assign # of k-points !!!
    #
    if { ![info exists wn(fs_nkp)] } {
	set klist $wn(dir)/$wn(filehead).klist
	if ![file exists $klist] {
	    ErrorDialog "File \"$wn(dir)/$wn(filehead).klist\" does not exists !!!"
	    exit_pr
	}
	# this is bad test
	set wn(fs_nkp) [exec head -1 $klist | awk "{print \$9}"]

	if { ! [string is integer $wn(fs_nkp)] || $wn(fs_nkp) == "" } {
	    set wn(fs_nkp) 2000
	}
    }
    if { ![info exists wn(fs_shift)] } {
	set wn(fs_shift) 0
    }
    
    #-----------------------------------------
    # 1st part of WIEN fermi TASK-LIST
    #-----------------------------------------
    # This is the definition of the widgets::
    #-----------------------------------------
    # LABEL:       "Number of k-points: "
    # ENTRY:       wn(fs_nkp)
    # BUTTON:      "Generate k-mesh"            COM:: wnGenKMesh
    
    #set widlist [list \
    #	    [list label -text "Number of k-points: "] \
    #	    [list entry -textvariable wn(fs_nkp) -width 20 -relief sunken] \
    #	    [list checkbutton -text "Allow the shift of k-mesh" \
    #	    -variable wn(fs_shift)] \
    #	    [list button -text "Generate k-mesh" -command wnGenKMesh] \
    #	    [list button -text "Change to unit 5 in $wn(filehead).in1$wn(c)" \
    #	    -command "wnTo5In1 $wn(filehead).in1$wn(c)" ] \
    #	    [list button -text "Edit $wn(filehead).in1$wn(c) file (optional)" \
    #	    -command "xcEditFile $wn(dir)/$wn(filehead).in1$wn(c)"]]
    
    set widlist [subst { 
	{label -text "Number of k-points: "}
	{entry -textvariable wn(fs_nkp) -width 20 -relief sunken}
	{label -text "NOTE: shift of k-mesh is not supported !!!\n(Fermi surfaces calculated with the shited k-mesh\nwill be erratically displayed !!!)"}
	{button -text "Generate k-mesh" -command wnGenKMesh}
    }]
	    
    #
    # Widget section !!!
    #
    foreach wid $widlist {
	set w [WidgetName $cf] 
	set com [concat [lindex $wid 0] $w [lrange $wid 1 end]]
	eval $com
	pack $w -side top -expand 1 -fill x -padx 2m -pady 1m	
    }

    #
    # if SPIN-POLARIZED make a NoteBook for up- & down-spin
    #
    if $wn(spin_polarized) {
	set nb [NoteBook [WidgetName $cf]]
	pack $nb -side top -expand 1 -fill x -padx 2m -pady 3m
	$nb insert 0 spin_up -text "Spin type: up"
	$nb insert 1 spin_dn -text "Spin type: dn"
	$nb raise spin_up

	set fup [$nb getframe spin_up]
	wnFSSetSpinFlag $fup up
	set fdn [$nb getframe spin_dn]
	wnFSSetSpinFlag $fdn dn 

    } else {
	wnFSSetSpinFlag $cf
    }

    #-----------------------------------------
    # 3rd part of WIEN fermi TASK-LIST
    #-----------------------------------------
    # This is the definition of the widgets::
    #-----------------------------------------
    # BUTTON:      "Exit"

    set e [button [WidgetName $cf] -text Exit -command exit_pr]
    pack $e -side top -expand 0 -padx 10m -pady 3m    
}


#
# we need this proc to handle correctly the spin polarized case
# (up- & down-spin each in its NoteBook tab --> proc called 
#  for each spin)
#
proc wnFSSetSpinFlag {cf {spin {}}} {
    global wn 
    
    set wn(options) {}
    if {$wn(c) != {}} {
	append wn(options) " -$wn(c)"
    } 
    if {$spin != {}} {
	append wn(options)   " -$spin"
    }

    set lapw2 [concat $wn(options) -fermi]

    set wn(updn)      $spin
    set wn(exe_flag)  [concat x lapw1 $wn(options)]
    set wn(lapw2)     [concat x lapw2 $lapw2]

    #-----------------------------------------
    # 2nd part of WIEN fermi TASK-LIST
    #-----------------------------------------

    #-----------------------------------------
    # This is the definition of the widgets::
    #-----------------------------------------
    # CHECKBUTTON: "Parallel execution of lapw1/lapw2"
    # BUTTON:      "Calculate Eigenvalues (x lapw1)"  COM:: wnRunWIEN
    # BUTTON:      "Calculate Fermi energy (x lapw2)" COM:: wnRunWIEN
    # BUTTON:      "Render Fermi Surface"             COM:: wnFSGo

    if { ! [info exists wn(checkbutton_parallel)] } {
	set wn(checkbutton_parallel) 0
    }
    set cb [checkbutton [WidgetName $cf] -text "Parallel execution of lapw1/lapw2" -variable wn(checkbutton_parallel)]
    pack $cb -side top -expand 1 -fill x -padx 2m -pady 1m

    set textlist [list \
		      "Calculate Eigenvalues \[x lapw1 $wn(options)\]" \
		      "Calculate Fermi energy \[x lapw2 $lapw2\]" \
		      "Render Fermi Surface"]
    
    set kgen $wn(dir)/$wn(filehead).outputkgen
    set comlist \
	[list \
	     [list wnFS_ExecLapw $wn(exe_flag) {WIEN program is calculating the eigenvalues. It can take a lot of CPU time. PLEASE WAIT!!!} \
		  $wn(dir)/$wn(filehead).output1$wn(updn)] \
	     [list wnFS_ExecLapw $wn(lapw2) {WIEN program is calculating the Fermi energy. PLEASE WAIT!!!} $wn(dir)/$wn(filehead).output2$wn(updn)] \
	     [list wnFSGo $kgen $spin]]
    
    foreach text $textlist com $comlist {
	set b [button [WidgetName $cf] -text $text -command $com]
	pack $b -side top -expand 1 -fill x -padx 2m -pady 1m
    }
}

proc wnFS_ExecLapw {cmd display_text outputf} {
    global wn

    if { $wn(checkbutton_parallel) } {
	set wn(parallel) 1
	append cmd " -p"
    } else {
	set wn(parallel) 0
    }
    wnRunWIEN $cmd $display_text $outputf
}

#set textlist [list \
#	     "Generate k-mash
#	 "Run SCF cycle first!"\
#	 "Generate k-list with XCrysDen"\
#	 "Insert k-Points and change unit to 5 in $case.in1$c" \
#	 "Calculate Eigenvalues (x lapw1 $cmdopt)"\
#	 "Edit $case.insp (insert correct EF)"\
#	 "Calculate Bandstructure (x spaghetti $cmdopt)"\
#	 "Preview Bandstructure (ghostview $case.spaghetti${fspin}_ps)"\
#	 "Save Bandstructure"\
#	 "Reset unit to 4 in $case.in1$c" \
#	 ]
#
#set tasklist [list  \
#	 "title" \
#	 "xcrysden --wien_kpath $workdir" \
#	 "edit $case.in1$c" \
#	 "exe x lapw1 $cmdopt"\
#	 "edit $case.insp"\
#	 "exe x spaghetti $cmdopt"\
#	 "exe gv $case.spaghetti${fspin}_ps"\
#	 "savetask $task"\
#	 "edit $case.in1$c"\
#	 ]


