#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/advGeom2.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ---------------------------------------------------------------------------
# PART#1>>
#          ATOMSUBS
#          ^^^^^^^^
proc AtomSubs {} {
    global atomSub XCState doneAL
    
    #if XCState is not c95, we can not substitute an atom
    puts stdout $XCState(state)
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { return }
    
    if [winfo exists .atsubs] { return }
    # how many atoms to substitute
    OneEntryToplevel .atsubs "Substitute Atoms" "AtomSubs" \
	    "Number of atoms to be substituted" 3 atomSub(NSOST) posint 310 60
    # make a toplevel with enries for AtomSubs
    if [winfo exists .atmsub] { return }
    set t [xcToplevel .atmsub "Substitute Atoms" "AtomSubs" . 180 0 1]

    set fb [frame $t.fb1]
    set atomSub(sym) "BREAKSYM"
    RadioBut $fb "Symmetry option::" atomSub(sym) \
	    left left 1 1 "BREAKSYM" "KEEPSYMM"
    pack $fb -side bottom -expand true -padx 10 -pady 10 -anchor center
 
    # -------------------------------
    # SCROLLENTRIES --- SCROLLENTRIES
    ScrollEntries .atmsub $atomSub(NSOST) "Atom N.:" \
	{{Atomic Label:} {New Atomic Number:}} {LB NA} {posint nat} \
	    3 atomSub [list 2 \
	    [list {Periodic Table} scroll_ptableSelectNAT .atmsub atomSub NA] \
	    [list {Select AtomSubs} SelAtomLabel .atmsub atomSub LB \
	    {Select an atom to substitute it} \
	    {To select the atom click on the atom}]] 3

    button $fb.butok -text OK -command \
	    [list AtomLabOK $t atomSub]
    
    button $fb.butcan -text Cancel -command [list CancelProc $t]
    
    pack $fb.butok $fb.butcan -side left -expand 1 -padx 10 -pady 10
    tkwait window $t
}


proc SelAtomLabel {t type elem title text i} {
    global SelLabel atomSub atomRemo atomDisp
    
    # array..... name of global array
    # elem...... array's element
    # i......... global array is of type globvarname($elem,$i)

    set grab 0
    if { [grab current] == $t } {
	set grab 1
	catch { grab release $t }
    }

    set selwin [PreSel .selatom .mesa $title $text LabelSel 1]
    tkwait variable SelLabel(done)

    if { $type == "atomSub" } {
	# ATOMSUB
	set atomSub(${elem},${i}) $SelLabel(label)
    } elseif { $type == "atomRemo" } {
	# ATOMREMO
	set atomRemo(${elem},${i}) $SelLabel(label)
    } elseif { $type == "atomDisp" } {
	# ATOMDISP
	set atomDisp(${elem},${i}) $SelLabel(label)
    }
    if { $grab } {
	catch { grab $t }
    }
}    


proc LabelSel {{w {}}} {
    global select maxsel system
    
    WriteFile "$system(SCRDIR)/xc_tmp.$system(PID)" \
	"$select(Nat1) $select(X1) $select(Y1) $select(Z1)"

    set com [list exec $system(FORDIR)/atomlab 1 \
		 $system(SCRDIR)/xc_struc.$system(PID) \
		 $system(SCRDIR)/xc_tmp.$system(PID)]
    if { $system(c95_version) != 95 } {
	set com [concat $com cr98]
    }

    # now exec a atomsel program; catch the error if it occure

    $select(textWid) config -state normal

    xcDebug -stderr "Executing: $com"
    if { [catch {set label [string trim [eval $com]]} errmsg] } {
	error "error executing atomlab: $errmsg"
	$select(textWid) delete [expr $maxsel + 4].0 [expr $maxsel + 4].end
	$select(textWid) insert [expr $maxsel + 4].0 $errmsg
	file delete $system(SCRDIR)/xc_tmp.$system(PID)
	$select(textWid) config -state disabled
	return "Error"
    } else {
	$select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end
	$select(textWid) insert [expr $maxsel + 5].0 "Selected atom has label: $label"
	file delete $system(SCRDIR)/xc_tmp.$system(PID)
	$select(textWid) config -state disabled
	return [string trim $label]
    }
}


# ----------------------------------------------------------------------------
# PART#2 >>
#          ATOMREMO
#          ^^^^^^^^
proc AtomRemo {} {
    global atomRemo XCState doneAL
    
    #if XCState is not c95, we can not cut a cluster
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { return }
    
    if [winfo exists .atremo] { return}
    # how many atoms to remove
    OneEntryToplevel .atremo "Remove Atoms" "AtomRemo" \
	    "Number of atoms to be removed" 3 atomRemo(NL) posint 310 60
    # make a toplevel with enries for AtomRemo
    if [winfo exists .atmrem] { return }
    set t [xcToplevel .atmrem "Remove Atoms" "AtomRemo" . 180 0 1]
 
    set fb [frame $t.fb1]
    set fb2 [frame $t.fb2]
    set atomRemo(sym) "BREAKSYM"
    RadioBut $fb "Symmetry option::" atomRemo(sym) \
	    left left 1 1 "BREAKSYM" "KEEPSYMM"
    pack $fb2 $fb  -side bottom -expand true -padx 10 -pady 10 -anchor center \
	    -fill x

    # -------------------------------------
    # SCROLLENTRIES --- SCROLLENTRIES #####
    ScrollEntries .atmrem $atomRemo(NL) "Atom N.:" \
	    {{Atomic Label:}} LB posint \
	    3 atomRemo [list 1 \
	    [list {Select AtomRemo} SelAtomLabel $t atomRemo LB \
	    {Select an atom to remove it} \
	    {To select the atom click on the atom}]] 3

    button $fb2.butok -text OK -command \
	    [list AtomLabOK $t atomRemo]
    
    button $fb2.butcan -text Cancel \
	    -command [list CancelProc $t]
    
    pack $fb2.butok $fb2.butcan -side left \
	    -expand 1 -padx 10

    tkwait window $t
}


proc AtomLabOK {t type} {
    global doneAL varlist foclist err atomSub atomRemo atomDisp \
	    AdvGeom system nxdir nydir nzdir n_groupsel periodic

    puts stdout "ARRAY ELEMENTS:: [array names atomRemo]"
    
    check_var $varlist $foclist
    if $err {return}
    
    set n [xcAdvGeomState new]

    if { $type == "atomSub" } {
	puts stdout "ATOMSUBS:: $atomSub(NSOST) $atomSub(LB,1) $atomSub(NA,1)"
	set AdvGeom($n,atomSub)    "$atomSub(sym)\n"
	append AdvGeom($n,atomSub) "ATOMSUBS\n"
	append AdvGeom($n,atomSub) "$atomSub(NSOST)\n"
	for {set i 1} {$i <= $atomSub(NSOST)} {incr i} {
	    append AdvGeom($n,atomSub) "$atomSub(LB,$i) $atomSub(NA,$i)\n"
	}
    } elseif { $type == "atomRemo" } {
	puts stdout "SYMM:: $atomRemo(sym)"
	set AdvGeom($n,atomRemo) "$atomRemo(sym)\n"
	append AdvGeom($n,atomRemo) "ATOMREMO\n"
	append AdvGeom($n,atomRemo) "$atomRemo(NL)\n"
	for {set i 1} {$i <= $atomRemo(NL)} {incr i} { 
	    append AdvGeom($n,atomRemo) "$atomRemo(LB,$i) "
	}
	puts stdout "ATOMREMO:: $AdvGeom($n,atomRemo)"
	append AdvGeom($n,atomRemo) "\n"
    } elseif { $type == "atomDisp" } {
	set AdvGeom($n,atomDisp)    "$atomDisp(sym)\n"
	append AdvGeom($n,atomDisp) "ATOMDISP\n"
	append AdvGeom($n,atomDisp) "$atomDisp(NDISP)\n"
	for {set i 1} {$i <= $atomDisp(NDISP)} {incr i} {
	    append AdvGeom($n,atomDisp) "$atomDisp(LB,$i) $atomDisp(DX,$i) $atomDisp(DY,$i) $atomDisp(DZ,$i)\n"
	}
    }

    # run crystal95
    set input [MakeInput]
    xcDebug -debug "CRYSTAL INPUT:\n$input"	
    if { ![RunC95 $system(c95_integrals) {} $input] } {
	#tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	return
    }

    # this  is for UNDO/REDO
    switch -exact -- $type {
	atomSub  {GenCommUndoRedo "Substitute Atom"}
	atomRemo {GenCommUndoRedo "Remove Atom"}
	atomDisp {GenCommUndoRedo "Displace Atom"}
    }
    
    GenGeomDisplay 1
    destroy $t
}


# this is old routine (Thu Sep  2 11:11:18 CEST 1999)
#proc AtomLabOK {t type} {
#    global doneAL varlist foclist err atomSub atomRemo AdvGeom \
#	     system nxdir nydir nzdir n_groupsel periodic
#
#    puts stdout "ARRAY ELEMENTS:: [array names atomRemo]"
#    
#    check_var $varlist $foclist
#    if $err {return}
#    
#    set n [xcAdvGeomState new]
#
#    if { $type == "atomSub" } {
#	 puts stdout "ATOMSUBS:: $atomSub(NSOST) $atomSub(LB,1) $atomSub(NA,1)"
#	 set AdvGeom($n,atomSub)    "$atomSub(sym)\n"
#	 append AdvGeom($n,atomSub) "ATOMSUBS\n"
#	 append AdvGeom($n,atomSub) "$atomSub(NSOST)\n"
#	 for {set i 1} {$i <= $atomSub(NSOST)} {incr i} {
#	     append AdvGeom($n,atomSub) "$atomSub(LB,$i) $atomSub(NA,$i)\n"
#	 }
#    } elseif { $type == "atomRemo" } {
#	 puts stdout "SYMM:: $atomRemo(sym)"
#	 set AdvGeom($n,atomRemo) "$atomRemo(sym)\n"
#	 append AdvGeom($n,atomRemo) "ATOMREMO\n"
#	 append AdvGeom($n,atomRemo) "$atomRemo(NL)\n"
#	 for {set i 1} {$i <= $atomRemo(NL)} {incr i} { 
#	     append AdvGeom($n,atomRemo) "$atomRemo(LB,$i) "
#	 }
#	 puts stdout "ATOMREMO:: $AdvGeom($n,atomRemo)"
#	 append AdvGeom($n,atomRemo) "\n"
#    }
#    # run crystal95
#    set input [MakeInput]
#    puts stdout "CRYSTAL INPUT:\n $input"	
#    if ![RunC95 $system(c95_integrals) {} $input] {
#	 #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
#	 return
#    }
#    
#    GenGeomDisplay 1
#    destroy $t
#}


# ----------------------------------------------------------------------------
# PART#3 >>
#          ATOMINSE
#          ^^^^^^^^
proc AtomInse {} {
    global atomInse XCState doneAI
    
    #if XCState is not c95, we can not cut a cluster
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0} { return }

    if [winfo exists .atinse] { return }
    # how many atoms to add
    OneEntryToplevel .atinse "Insert Atoms" "AtomInse" \
	    "Number of atoms to be added" 3 atomInse(NINS) posint 310 60
    # make a toplevel with enries for AtomInse
    if [winfo exists .atins] { return }
    set t [xcToplevel .atmins "Insert Atoms" "AtomInse" . 180 0 1]
 
    # two bottom frames, one for "select AtomInse" & KEEPSYMM,
    # second for OK button 
    set fb [frame $t.fb1]
    set atomInse(how2) "Hole-Adding"
    RadioBut $fb "How to add the atom:" atomInse(how2) \
	    left left 1 1 "Hole-Adding" "Cell-Adding"  "Line-Adding"
    pack $fb -side top -expand true -padx 10 -pady 10 -anchor center
    set fb2 [frame $t.fb2]
    set fb3 [frame $t.fb3]
    pack $fb2 $fb3 -side top -expand true -fill both
    # -------------------------------
    # SCROLLENTRIES --- SCROLLENTRIES
    ScrollEntries $fb2 $atomInse(NINS) "Atom N.:" \
	    {{Atm Num:} X: Y: Z:} {NA X Y Z} {nat real real real} \
	    7 atomInse [list 2 \
	    [list {Select AtomInse} AddAtom $t] \
	    [list {Periodic Table} scroll_ptableSelectNAT $t atomInse NA]] 2

    set fb3f [frame $fb3.f]
    set fb3a [frame $fb3.f.a]
    set fb3b [frame $fb3.f.b]
    pack $fb3f -side left -expand true -fill both
    pack $fb3a $fb3b -side top -expand true -fill both

    set atomInse(sym)       "BREAKSYM"
    set atomInse(coor_type) "FRACTION"
    RadioBut $fb3a "Symmetry option::" atomInse(sym) \
	    left left 1 1 "BREAKSYM" "KEEPSYMM"
    RadioBut $fb3b "Coordinates unit::" atomInse(coor_type) \
	    left left 1 1 "FRACTION" "ANGSTROM"

    button $fb3.butok -text "OK" -command \
	    [list AddAtomOKCan ok $t]
    
    button $fb3.butcan -text Cancel -command \
	    [list AddAtomOKCan can $t]
    
    pack $fb3.butok $fb3.butcan -side left \
	    -expand 1 -padx 10 -pady 10

    tkwait window $t
}


proc AddAtom {t i} {
    global periodic addAtom SelHoleCL atomInse SelLine \
	    atompos system af bf cf
    
    set SelHoleCL(transl) 1

    # atomInse(how2) can be Hole-Adding or Cell-Adding
    if { $atomInse(how2) == "Hole-Adding" } {
	# here I will use "SelCentreHoleCL" proc that was designd 
        # for CutCluster
	PreSel .addatom .mesa "Coordinates of the atom to be added" \
	    "A \"hole\" is geometrical centre of several atoms!!!\n\
		To select the hole click on desired number of atoms" \
	    SelCentreHoleCL 15; # 15 is maximum allowed number of 
	                        # selected atoms
	tkwait variable SelHoleCL(done)
	if { [info exists SelHoleCL(transl)] } {unset SelHoleCL(transl)}

	if { $atomInse(coor_type) == "ANGSTROM" } {
	    set coor $SelHoleCL(centre)
	} else {
	    # factional units
	    set coor [GetFracCoor $SelHoleCL(centre)]
	}
	set atomInse(X,$i) [lindex $coor 0]
	set atomInse(Y,$i) [lindex $coor 1]
	set atomInse(Z,$i) [lindex $coor 2]	
    } elseif { $atomInse(how2) == "Cell-Adding" } {
	# if periodic(dim) == 0 --> return silently
	if { $periodic(dim) == 0 } { return }
	
	set tplw [xcToplevel .cell "Coordinates of the atom to be added" \
		"AtomInse" . 30 0]
	set f1 [frame $tplw.f1 -relief raised -bd 2]
	set f2 [frame $tplw.f2 -relief raised -bd 2]
	pack $f1 $f2 -side top -fill both -expand 1
	set af 0.3
	set bf 0.3
	set cf 0.3
	if { $periodic(dim) == 3 } {
	    Entries $f1 {fract-A: fract-B: fract-C:} {af bf cf} 8 1
	} elseif { $periodic(dim) == 2 } {
	    Entries $f1 {fract-A: fract-B: Z:} {af bf cf} 8 1
	} elseif { $periodic(dim) == 1 } {
	    Entries $f1 {fract-A: Y: Z:} {af bf cf} 8 1
	} else {
	    tk_dialog .dialog ERROR \
		    "Wrong dimensionality of structure. Goodbye!!!" \
		    error 0 OK
	    return
	}
	
	# initialise atomadd procedure

	set atompos [xc_atomadd .mesa begin]

	set f11 [frame $f1.1 -relief flat]
	set l1 [label $f11.l1 -text "Current ATOMINSE position:" -relief flat]
	set l2 [label $f11.l2 -textvar atompos -relief sunken -bd 2 -width 45]
	pack $f11 -side top -fill both -expand 1
	pack $l1 $l2 -side left -padx 5 -pady 5
	
	set upd [button $f2.upd -text "Update" -command [list CellAdding upd]]
	set ok  [button $f2.ok -text "OK" -command [list CellAdding ok]]
	set can [button $f2.can -text "Cancel" -command [list CellAdding can]]
	pack $upd $ok $can -side left -expand 1 -padx 5 -pady 5
	
	tkwait variable addAtom(done)

	if { $atomInse(coor_type) == "ANGSTROM" } {
	    set atomInse(X,$i) [lindex $atompos 0]
	    set atomInse(Y,$i) [lindex $atompos 1]
	    set atomInse(Z,$i) [lindex $atompos 2]
	} else {
	    # factional units
	    set atomInse(X,$i) $af
	    set atomInse(Y,$i) $bf
	    set atomInse(Z,$i) $cf
	}
	destroy $tplw
    } elseif { $atomInse(how2) == "Line-Adding" } { 
	# there must also be some entry for specifying the fraction of 
	# line from first atom --> that is the adding atom position
	PreSel .addatom .mesa "Coordinates of the atom to be added" \
		"For Line-Adding it is necessary to select two atoms\n \
		& to specify the fraction of line-length from first atom \
		in FRACTIONAL units" \
		LineSel 2
	tkwait variable SelLine(done)

	if { $atomInse(coor_type) == "ANGSTROM" } {
	    set coor $SelLine(coor)
	} else {
	    # factional units
	    set coor [GetFracCoor $SelLine(coor)]
	}
	set atomInse(X,$i) [lindex $coor 0]
	set atomInse(Y,$i) [lindex $coor 1]
	set atomInse(Z,$i) [lindex $coor 2]	
    }
}


proc AddAtomOKCan {type w} {
    global varlist foclist err atomInse AdvGeom \
	    system nxdir nydir nzdir n_groupsel periodic
    
    if { $type == "ok" } { 
	puts stdout "ARRAY ELEMENTS:: [array names atomRemo]"
	
	check_var $varlist $foclist
	if $err {return}
	
	set n [xcAdvGeomState new]

	puts stdout "ATOMINSE:: $atomInse(NINS) \
		$atomInse(NA,1) $atomInse(X,1) $atomInse(Y,1) $atomInse(Z,1)"
	set AdvGeom($n,atomInse) "$atomInse(coor_type)\n$atomInse(sym)\n"
	append AdvGeom($n,atomInse) "ATOMINSE\n"
	append AdvGeom($n,atomInse) "$atomInse(NINS)\n"
	for {set i 1} {$i <= $atomInse(NINS)} {incr i} {
	    append AdvGeom($n,atomInse) \
		    "$atomInse(NA,$i) $atomInse(X,$i) $atomInse(Y,$i) \
		    $atomInse(Z,$i)\n"
	}
	
	# run crystal95
	set input [MakeInput]
	xcDebug -debug "CRYSTAL INPUT:\n$input"
	if { ![RunC95 $system(c95_integrals) {} $input] } {
	    #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	    return
	}
	
	# this  is for UNDO/REDO
	GenCommUndoRedo "Insert Atom"
	
	GenGeomDisplay 1
    } 
    
    # this is for type == ok & can
    if [winfo exists $w] { destroy $w } 
}


proc CellAdding {type} {
    global addAtom atompos af bf cf
    
    # updating::
    if { $type == "upd" } {
	puts stdout "CellAdding upd ATOMPOS: $af $bf $cf"
	set atompos [xc_atomadd .mesa update $af $bf $cf]
	puts stdout "CellAdding upd ATOMPOS: $atompos"
    } elseif { $type == "ok" } {
	# OK::
	if [info exists addAtom(done)] { unset addAtom(done) }    
	set atompos [xc_atomadd .mesa update $af $bf $cf]
	puts stdout "CellAdding ok ATOMPOS: $atompos"
	xc_atomadd .mesa clean
	set addAtom(done) 1
    } elseif { $type == "can" } {
	# Cancel::
	xc_atomadd .mesa clean
	destroy .cell
    }
}



# ----------------------------------------------------------------------------
# PART#4 >>
#          SUPERCELL
#          ^^^^^^^^^
proc SuperCell {} {
    global AdvGeom superCell periodic

    # if XCState is not c95, we can not generate a supercell
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { 
	return 
    }
    
    # SuperCell option is available just for periodic systems
    if { $periodic(dim) == 0 } {
	return
    }

    #
    # initialize SuperCell Option
    #
    SetWatchCursor
    xc_supercell .mesa init
    ResetCursor

    set t [xcToplevel [WidgetName] "Generation of a SuperCell" \
	    "SuperCell" . 180 0 1]

    set f1 [frame $t.f1 -highlightthickness 0]
    set f2 [frame $t.f2 -relief raised -bd 2 -highlightthickness 0]
    pack $f1 $f2 -side top -fill x

    ############
    # FRAME #1 #
    ############
    set f21 [frame $f2.1 -relief raised -bd 1]
    set f22 [frame $f2.2 -relief raised -bd 1]
    set f2a [frame $f2.a]
    set f2b [frame $f2.b]
    pack $f21 -side left -fill both -expand 1
    pack $f22 -side left -fill both
    pack $f2a -side left -fill both -in $f21 -padx 5 -pady 5

    set b1 [button $f1.b1 \
	    -text "Specify an Expansion Matrix" \
	    -bd 3 \
	    -highlightthickness 0 \
	    -command [list SuperCell_Show matrix $f2a $f2b $f21 $f1.b1 $f1.b2]]
    set b2 [button $f1.b2 \
	    -text "Elongate Unit-Cell Vectors" \
	    -bd 1 \
	    -highlightthickness 0 \
	    -command [list SuperCell_Show elongate $f2a $f2b $f21 $b1 $f1.b2]]
    pack $b1 $b2 -side left -fill y

    ############
    # FRAME #2 #
    ############
    #
    # frame $f2a
    set superCell(type) matrix
    set superCell(matrix_varlist) {}
    set superCell(matrix_foclist) {}
    label $f2a.lt -text "New_Vectors == Expansion_Matrix x Vectors" -anchor w
    grid $f2a.lt -column 0 -row 0 -columnspan 4 -sticky w -pady 10
    label $f2a.ll -text "Expansion Matrix =="
    if { $periodic(dim) == 3 } {
	set row 2
    } else {
	set row 1
    }
    grid $f2a.ll -column 0 -row $row
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	for {set j 1} {$j <= $periodic(dim)} {incr j} {	
	    set e [entry $f2a.e${i}_$j \
		    -textvariable superCell(e$i,e$j) -width 5]
	    lappend superCell(matrix_varlist) [list superCell(e$i,e$j) real]
	    lappend superCell(matrix_foclist) $e
	    grid $e -column $j -row $i -padx 3 -pady 3
	}
    }
    focus [lindex $superCell(matrix_foclist) 0]
    tkwait visibility $f2a.e1_1
    set w [winfo width  $f21]
    set h [winfo height $f21]
    pack propagate $f21 false 

    #
    # frame $f2b
    set superCell(elongate_varlist) {}
    set superCell(elongate_foclist) {}
    label $f2b.lt -text "New_Vector(i) == factor(i) Vector(i)" -anchor w
    grid $f2b.lt -column 0 -row 0 -columnspan 4 -sticky w -pady 10
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	label $f2b.l$i -text "factor($i)   =="
	grid $f2b.l$i -column 0 -row $i -sticky w
	set e [entry $f2b.e$i \
		-textvariable superCell(e$i) -width 5]
	lappend superCell(elongate_varlist) [list superCell(e$i) real]
	lappend superCell(elongate_foclist) $e
	grid $e -column 1 -row $i -padx 3 -pady 3 -sticky w
    }    

    #
    # frame $f22 --> Test It/OK/Cancel
    #
    proc SuperCell_Can t {
	global superCell
	unset superCell
	#
	# clear the SuperCell OPTION
	#
	SetWatchCursor
	xc_supercell .mesa clear
	# go back to previous render state
	CalStru
	ResetCursor
	xcSwapBuffers
	CancelProc $t
    }
    set testit [button $f22.test -text "Test It" \
	    -command [list SuperCell_TestIt $f22.ok]]
    set ok     [button $f22.ok   -text "OK" -command [list SuperCell_OK $t]]
    set can    [button $f22.can  -text "Cancel" \
	    -command [list SuperCell_Can $t]]
    pack $testit $ok $can -side top -padx 5 -pady 5 -fill x
    $ok configure -state disabled

    xcSwapBuffers
}


proc SuperCell_Show {type fmatrix felongate parent bmatrix belongate} {
    global superCell

    if { $type == "matrix" } {
	set superCell(type) matrix
	$bmatrix   config -bd 3
	$belongate config -bd 1
	pack forget $felongate
	pack $fmatrix -side left -fill both -in $parent -padx 5 -pady 5
	focus [lindex $superCell(matrix_foclist) 0]
    } else {
	set superCell(type) elongate
	$bmatrix   config -bd 1
	$belongate config -bd 3
	pack forget $fmatrix
	pack $felongate -side left -fill both -in $parent -padx 5 -pady 5
	focus [lindex $superCell(elongate_foclist) 0]
    }
}


proc SuperCell_TestIt okbut {

    $okbut configure -state normal

    SetWatchCursor
    update
    if { ![SuperCell_ValidEMatrix] } {
    	ResetCursor
    	xcSwapBuffers
    	return 0
    }
    
    if { ![SuperCell_Run95] } {
    	ResetCursor
    	xcSwapBuffers
    	return 0
    }
    
    # everything looks OK
    GenGeomDisplay 1
    xc_supercell .mesa testit
    xcAdvGeomState delete
    ResetCursor
    xcSwapBuffers

    return 1
}


proc SuperCell_OK tplw {
    global superCell

    SetWatchCursor
    update
    if ![SuperCell_ValidEMatrix] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    if ![SuperCell_Run95] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    # this  is for UNDO/REDO
    GenCommUndoRedo "SuperCell"

    xc_supercell .mesa clear
    ResetCursor
    xcSwapBuffers
    destroy $tplw
}



proc SuperCell_ValidEMatrix {} {
    global superCell periodic
    
    set type $superCell(type)
    if ![check_var $superCell(${type}_varlist) $superCell(${type}_foclist)] {
	return 0
    }
    
    xcDebug -debug "Periodic(dim) == $periodic(dim)"
    if { $superCell(type) == "matrix" } {
	for {set i 1} {$i <= $periodic(dim)} {incr i} {
	    for {set j 1} {$j <= $periodic(dim)} {incr j} {
		set superCell($i,$j) $superCell(e$i,e$j)
	    }
	}
    } else {
	for {set i 1} {$i <= $periodic(dim)} {incr i} {
	    for {set j 1} {$j <= $periodic(dim)} {incr j} {
		if { $i != $j } {
		    set superCell($i,$j) 0.0
		} else {
		    if { $superCell(e$i) == 0 } {
			tk_dialog [WidgetName] ERROR \
				"ERROR !\nYou have badly specified the $i factor. Please try again !" error 0 OK
			focus [lindex $superCell(elongate_foclist) \
				[expr $i - 1]]
			return 0
		    }
		    set superCell($i,$i) $superCell(e$i)
		}
	    }
	}
    }

    return 1
}


proc SuperCell_Run95 {} {
    global superCell periodic system AdvGeom

    # now set the AdvGeom variable
    set n [xcAdvGeomState new]
    set AdvGeom($n,superCell) "SUPERCELL\n"
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	for {set j 1} {$j <= $periodic(dim)} {incr j} {
	    append AdvGeom($n,superCell) "$superCell($i,$j) "
	}
	append AdvGeom($n,superCell) "  "
    }
    append AdvGeom($n,superCell) \n

    set   input    [MakeInput]
    set   output   $system(SCRDIR)/xc_output.$system(PID)
    set   inp      $system(SCRDIR)/xc_tmp.$system(PID)
    set   fileID   [open $inp w]
    puts  $fileID  $input
    flush $fileID
    close $fileID

    xcDebug -debug "CRYSTAL INPUT:\n$input"

    return [RunC95_advGeom $inp SUPERCELL \
		"super-cell was badly choosen or the upper limit for the number of atoms in SUPERCELL was exceeded."]
}


# ---------------------------------------------------------------------------
# PART#5>>
#          ATOMDISP
#          ^^^^^^^^
proc AtomDisp {} {
    global atomDisp XCState doneAL
    
    # if XCState is not c95, we can not substitute an atom
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { return }
    
    # how many atoms to substitute
    if [winfo exists .atdisp] { return }
    OneEntryToplevel .atdisp "Displace Atoms" "atomDisp" \
	    "Number of atoms to be displaced" 3 atomDisp(NDISP) posint 310 60
    # make a toplevel with enries for atomDisp
    if [winfo exists .atmdis] { return }
    set t [xcToplevel .atmdis "Displace Atoms" "atomDisp" . 180 0 1]

    set fb [frame $t.fb1]
    set atomDisp(sym) "BREAKSYM"
    RadioBut $fb "Symmetry option::" atomDisp(sym) \
	    left left 1 1 "BREAKSYM" "KEEPSYMM"
    pack $fb -side bottom -expand true -padx 10 -pady 10 -anchor center
 
    # -------------------------------------
    # SCROLLENTRIES --- SCROLLENTRIES #####
    ScrollEntries .atmdis $atomDisp(NDISP) "Atom N.:" \
	    {{Atomic Label:} DX: DY: DZ:} {LB DX DY DZ} \
	    {posint real real real} 3 atomDisp [list 1 \
	    [list {Select AtomDisp} SelAtomLabel $t atomDisp LB \
	    {Select an atom to substitute} \
	    {To select an atom click on the atom}]] 3

    button $fb.butok -text OK -command \
	    [list AtomLabOK $t atomDisp]
    
    button $fb.butcan -text Cancel -command [list CancelProc $t]
    
    pack $fb.butok $fb.butcan -side left \
	    -expand 1 -padx 10 -pady 10

    tkwait window $t
}

# ----------------------------------------------------------------------------
# PART#6 >>
#          ELASTIC
#          ^^^^^^^
proc Elastic {} {
    global AdvGeom elastic periodic

    # if XCState is not c95, we can not do ELATSIC
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { 
	return 
    }
    
    # Elastic option is available just for periodic systems
    # I think just for CRYSTALS ???
    if { $periodic(dim) == 0 } {
	return
    }

    #
    # initialize SuperCell Option 
    #    (I will do the same feedback as for SuperCell)
    SetWatchCursor
    xc_supercell .mesa init
    ResetCursor

    set t [xcToplevel [WidgetName] "Elastic Deformation" \
	    "Elastic" . 180 0 1]

    set elastic(type) Z_vol
    set elastic(varlist) {}
    set elastic(foclist) {}

    set f1  [frame $t.1 -relief raised -bd 1]
    set f2  [frame $t.2 -relief raised -bd 1]
    set f1a [frame $f1.a -relief groove  -bd 2]
    set f1b [frame $f1.b -relief groove  -bd 2]
    pack $f1 -side left -fill both -expand 1
    pack $f2 -side left -fill both
    pack $f1a $f1b -side top -fill both -padx 1m -pady 1m

    label $f1a.1 -text \
	    "Deformation through \"Z\" matrix is defined as:\nL' == L + Z"
    label $f1a.2 -text \
	    "Deformation through \"e\" matrix is defined as:\nL' == (I + e)L"
    pack $f1a.1 $f1a.2 -side top -fill both

    set textlist {
	{Z matrix deformation -- volume CONSERVING}
	{Z matrix deformation -- NOT volume conserving}
	{e matrix deformation -- volume CONSERVING}
	{e matrix deformation -- NOT volume conserving}
    }
    set valuelist { Z_vol Z e_vol e }
    foreach text $textlist value $valuelist i {1 2 3 4} {
	radiobutton $f1a.r$i -variable elastic(type) \
		-text $text -value $value -anchor sw
	pack $f1a.r$i -side top -fill both -padx 10 -pady 5
    }

    label $f1b.l -text "Specify Matrix:"
    grid $f1b.l -row 0 -padx 3 -pady 3 -columnspan $periodic(dim) -sticky w
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	for {set j 1} {$j <= $periodic(dim)} {incr j} {	
	    set e [entry $f1b.e${i}_$j \
		    -textvariable elastic(e$i,e$j) -width 5]
	    lappend elastic(varlist) [list elastic(e$i,e$j) real]
	    lappend elastic(foclist) $e
	    grid $e -column $j -row $i -padx 3 -pady 3
	}
    }
    focus [lindex $elastic(foclist) 0]
    tkwait visibility $f1b.e1_1

    #
    # frame $f2 --> Test It/OK/Cancel
    #
    proc Elastic_Can t {
	global elastic
	unset elastic
	#
	# clear the SuperCell OPTION
	#
	SetWatchCursor
	xc_supercell .mesa clear
	# go back to previous render state
	CalStru
	ResetCursor
	xcSwapBuffers
	CancelProc $t
    }
    set testit [button $f2.test -text "Test It" \
	    -command [list Elastic_TestIt $f2.ok]]
    set ok     [button $f2.ok   -text "OK" -command [list Elastic_OK $t]]
    set can    [button $f2.can  -text "Cancel" \
	    -command [list Elastic_Can $t]]
    pack $testit $ok $can -side top -padx 5 -pady 5 -fill x
    $ok configure -state disabled

    xcSwapBuffers
}

proc Elastic_TestIt okbut {

    $okbut configure -state normal

    SetWatchCursor
    update
    if {![Elastic_ValidEMatrix] } {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    if { ![Elastic_Run95] } {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    # everything looks OK
    GenGeomDisplay 1
    xc_supercell .mesa testit
    xcAdvGeomState delete
    ResetCursor
    xcSwapBuffers

    return 1
}


proc Elastic_OK tplw {
    global elastic

    SetWatchCursor
    update
    if { ![Elastic_ValidEMatrix] } {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    if { ![Elastic_Run95] } {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    # this  is for UNDO/REDO
    GenCommUndoRedo "Elastic Deformation"

    xc_supercell .mesa clear
    ResetCursor
    xcSwapBuffers
    destroy $tplw
}



proc Elastic_ValidEMatrix {} {
    global elastic periodic
    
    set type $elastic(type)
    if ![check_var $elastic(varlist) $elastic(foclist)] {
	return 0
    }
    
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	for {set j 1} {$j <= $periodic(dim)} {incr j} {
	    set elastic($i,$j) $elastic(e$i,e$j)
	}
    }
    
    return 1
}


proc Elastic_Run95 {} {
    global elastic periodic system AdvGeom

    switch -exact -- $elastic(type) {
	Z_vol { set idef +1 }
	Z     { set idef -1 }
	e_vol { set idef +2 }
	e     { set idef -2 }
    }

    # now set the AdvGeom variable
    set n [xcAdvGeomState new]
    set AdvGeom($n,elastic) "ELASTIC\n$idef\n"
    for {set i 1} {$i <= $periodic(dim)} {incr i} {
	for {set j 1} {$j <= $periodic(dim)} {incr j} {
	    append AdvGeom($n,elastic) "$elastic($i,$j) "
	}
	if { $i < $periodic(dim) } {
	    append AdvGeom($n,elastic) "\n"
	}
    }
    append AdvGeom($n,elastic) \n

    set   input    [MakeInput]
    set   output   $system(SCRDIR)/xc_output.$system(PID)
    set   inp      $system(SCRDIR)/xc_tmp.$system(PID)
    set   fileID   [open $inp w]
    puts  $fileID  $input
    flush $fileID
    close $fileID

    xcDebug -debug "CRYSTAL INPUT:"
    xcDebug -debug "--------------\n$input"

    return [RunC95_advGeom $inp ELASTIC "error occured while executing ELASTIC option."]
}


# ----------------------------------------------------------------------------
# PART#7 >>
#          ROTATE
#          ^^^^^^
proc RotFrame {} {
    global AdvGeom rotFrame periodic

    # if XCState is not c95, we can not do ELATSIC
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { 
	return 
    }
    
    # ROTATE option is available just for periodic systems
    # I think just for CRYSTALS ???
    if { $periodic(dim) == 0 } {
	return
    }

    #
    # initialize SuperCell Option 
    #    (I will do the same feedback as for SuperCell)
    SetWatchCursor
    xc_supercell .mesa init
    ResetCursor

    set t [xcToplevel [WidgetName] "Rotate Cartesian Frame" \
	    "Rotate" . 180 0 1]

    set rotFrame(varlist) {}
    set rotFrame(foclist) {}

    set f1  [frame $t.1 -relief raised -bd 1]
    set f2  [frame $t.2 -relief raised -bd 1]
    pack $f1 -side left -fill both -expand 1
    pack $f2 -side left -fill both

    label $f1.1 -text \
	    "Specify Miller indices of the basal layer of the new 3D unit cell"
    pack $f1.1 -side top -fill both

    set rotFrame(foclist) \
	    [Entries $f1 {I: J: K:} {rotFrame(I) rotFrame(J) rotFrame(K)} 4 1]
    set whead [string trimright $rotFrame(foclist) 1]
    lappend rotFrame(foclist) ${whead}2 ${whead}3
    set rotFrame(varlist) {
	{rotFrame(I) int}
	{rotFrame(J) int}
	{rotFrame(K) int}
    }

    #
    # frame $f2 --> Test It/OK/Cancel
    #
    proc RotFrame_Can t {
	global rotFrame
	unset rotFrame
	#
	# clear the SuperCell OPTION
	#
	SetWatchCursor
	xc_supercell .mesa clear
	# go back to previous render state
	CalStru
	ResetCursor
	xcSwapBuffers
	CancelProc $t
    }
    set testit [button $f2.test -text "Test It" \
	    -command [list RotFrame_TestIt $f2.ok]]
    set ok     [button $f2.ok   -text "OK" -command [list RotFrame_OK $t]]
    set can    [button $f2.can  -text "Cancel" \
	    -command [list RotFrame_Can $t]]
    pack $testit $ok $can -side top -padx 5 -pady 5 -fill x
    $ok configure -state disabled

    xcSwapBuffers
}


proc RotFrame_TestIt okbut {

    $okbut configure -state normal

    SetWatchCursor
    update
    if ![RotFrame_Valid_IJK] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    if ![RotFrame_Run95] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    # everything looks OK
    GenGeomDisplay 1
    xc_supercell .mesa testit
    xcAdvGeomState delete
    ResetCursor
    xcSwapBuffers

    return 1
}


proc RotFrame_OK tplw {
    global rotFrame

    SetWatchCursor
    update
    if ![RotFrame_Valid_IJK] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    if ![RotFrame_Run95] {
	ResetCursor
	xcSwapBuffers
	return 0
    }

    # this  is for UNDO/REDO
    GenCommUndoRedo "Rotate Cartesian Frame"

    xc_supercell .mesa clear
    ResetCursor
    xcSwapBuffers
    destroy $tplw
}



proc RotFrame_Valid_IJK {} {
    global rotFrame periodic
    
    if ![check_var $rotFrame(varlist) $rotFrame(foclist)] {
	return 0
    }
    return 1
}


proc RotFrame_Run95 {} {
    global rotFrame periodic system AdvGeom

    set n [xcAdvGeomState new]
    set AdvGeom($n,rotate) "ROTATE\n$rotFrame(I) $rotFrame(J) $rotFrame(K)\n"

    set   input    [MakeInput]
    set   output   $system(SCRDIR)/xc_output.$system(PID)
    set   inp      $system(SCRDIR)/xc_tmp.$system(PID)
    set   fileID   [open $inp w]
    puts  $fileID  $input
    flush $fileID
    close $fileID

    xcDebug -debug "CRYSTAL INPUT:"
    xcDebug -debug "--------------\n$input"

    return [RunC95_advGeom $inp ROTATE "an error occured while executing ROTATE option"]
}
