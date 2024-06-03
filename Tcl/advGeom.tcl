#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/advGeom.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# --------------------------------------------------------------------------- #
#                 ADV_GEOM OPTIONS --- ADV_GEOM OPTIONS                       #
# in this file: SLAB
#               CLUSTER options are implemented
# --------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------
# PART#1>>
#          CAT A SLAB
#          ^^^^^^^^^^
proc CutSlab {} {
    global H K L ISUP NL cutslabfoclist species nxdir nydir XCState periodic

    # there must be manual & 'graphical' way of doing it
    # slab could be cut only if $species == "Crystal"
 
    #if XCState is not c95, we can not cut a slab
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 } { return }
    
    puts stderr "SPECIES: $species"
    puts stderr "DIM:     $periodic(dim)"

    set ok 1
    if { $species == "external" } {
	if { $periodic(dim) != 3 } {
	    set ok 0
	}
    } else {
	if { $species != "crystal" } {
	    set ok 0
	}
    }
    if { ! $ok } {
	tk_dialog .cutslabwarning "Warning" "WARNING: \"Cut a Slab\" \
		option can be applied only for CRYSTALS !!!" warning 0 OK
	return
    }

    if [winfo exists .cutslab] { return }
    set t [xcToplevel .cutslab "Cut a Slab" "Cut a Slab" . 100 0 1]

    set cutslabfoclist {}

    set f2 [frame $t.f2 -relief raised -bd 1]
    set f3 [frame $t.f3 -relief raised -bd 1]
    pack $f2 $f3 -side top -fill both -expand true

    # FRAME #1 has been deleted
    # in FRAME #2 goes three entries & button
    set fr2 [frame $f2.fr -relief groove -bd 2]
    pack $fr2 -side top -padx 7 -pady 7 -ipadx 7 -ipady 7 -fill x
    set l1 [label $fr2.l1 -text "Miller idices of a plane:"]
    pack $l1 -side top -expand 1
    Entries $fr2 "H: K: L:" "H K L" 3
    append cutslabfoclist " $fr2.frame.entry1 $fr2.frame.entry2 \
	    $fr2.frame.entry3"
    focus $fr2.frame.entry1
    set b [button $fr2.b -text "Select the Plane" \
	    -command [list SelPlane $t]]
    pack $b -expand 1
    set frm2 [frame $f2.frm -relief groove -bd 2]
    pack $frm2 -side top -padx 7 -pady 7 -ipadx 7 -ipady 7 -fill x
    set f21 [frame $frm2.1]
    pack $f21 -side top -fill both -expand 1
    Entries $f21 {{Label of Surface layer:}} ISUP 3
    
    set f22 [frame $frm2.2]
    pack $f22 -side top -fill both -expand 1
    Entries $f22 {{Number of layers:}} NL 3

    # here will be specification of number of cells to be drawn in X&Y dir
    set frm21 [frame $f2.frm1 -relief groove -bd 2]
    pack $frm21 -side top -padx 7 -pady 7 -ipadx 7 -ipady 7 -fill x
    set l2 [label $frm21.lbl -text "Specify number of cells to be drawn in"]
    pack $l2 -side top -padx 3 -pady 7
    set f23 [frame $frm21.3]
    pack $f23 -side top -fill both -expand 1
    Entries $f23 {{in X-dir:} {in Y-dir:}} {nxdir nydir} 3
    
    append cutslabfoclist " $f21.frame.entry1 $f22.frame.entry1\
	    $f23.frame.entry1 $f23.frame.entry2"
    
    # in FRAME #3 there will be OK & CANCEL BUTTON
    # make a ring around OK button
    set okf [frame $f3.okf -relief sunken -bd 2]
    set ok  [button $okf.ok -text OK -command [list CutSlabOK $t]]
    set can [button $f3.can -text Cancel -command [list Cancel $t]]
    pack $okf -side left -expand 1 -padx 5m -pady 2m 
    pack $can -side left -expand 1 -padx 5m -pady 2m -ipadx 1m -ipady 1m
    pack $ok -padx 1m -pady 1m -ipadx 1m -ipady 1m
}


proc CutSlabOK {t} {
    global H K L ISUP NL cutslabfoclist nxdir nydir AdvGeom err radio

    puts stdout $cutslabfoclist
    # first we must check a variables
    check_var {{H int} {K int} {L int} {ISUP int} {NL posint} {nxdir posint}\
	    {nydir posint}} $cutslabfoclist
    # check_var return $err; if $err --> mistake -> RETURN
    if $err { return }

    # if we come so far --> everything is OK
    set n [xcAdvGeomState new]
    set AdvGeom($n,slab) "$H $K $L \n$ISUP $NL\n"
    destroy $t

    set radio(cellmode) "prim"

    # t.k
    global XCState
    xcDebug "state is: $XCState(state)"
    
    # this  is for UNDO/REDO
    GenCommUndoRedo "Cut a SLAB"
    
    CalStru 
}


proc Cancel {t} {
    destroy $t
    uplevel return
}


proc SelPlane {t} {
    global H K L ISUP NL cutslabfoclist species nxdir nydir SelPlane
    
    set grab 0
    if { [grab current] == "$t" } {
	set grab 1
	catch { grab release $t }
    }
    puts stdout "GRAB IS> [grab status $t]"
    PreSel .sel_plane .mesa "Select a Plane" \
            "For Plane Click on three atoms" PlaneSel 3
    tkwait variable SelPlane(done)
    puts stdout "Plane> $SelPlane(plane)"
    set H [lindex $SelPlane(plane) 0]
    set K [lindex $SelPlane(plane) 1]
    set L [lindex $SelPlane(plane) 2]
    puts stdout "HKL> $H $K $L"

    if { $grab } { 
	catch { grab $t }
    }
}    


# -----------------------------------------------------------------------------
# PART#2>>
#          CUT A CLUSTER
#          ^^^^^^^^^^^^^
proc CutCluster {} {
    global cutCL cutclusfoclist species XCState anal coorn brdmod 
    
    # cutCL .... all variables needed for crystal's CLUSTER data are herein
    # cutclusfoclist .... list of Entries; numbers in Entries must be checked
    #                     when OK button is pressed
    
    #if XCState is not c95, we can not cut a cluster
    if { [xcIsActive c95] == 0 || [xcIsActive render] == 0 }  { return }
    
    # species must be a periodic one
    if { $species == "molecule" } {
	tk_dialog .cutslabwarning "Warning" "WARNING:\"Cut a Cluster\" \
		option can't be applied for MOLECULES !!!" warning 0 OK
	return
    }
    
    if [winfo exists .cutslab] { return }
    set t [toplevel .cutclus]
    xcPlace . .cutclus 530 25
    wm title $t "Cut a Cluster"
    wm iconname $t "Cut a Cluster"
    
    set f1 [frame $t.f1 -relief groove -bd 2]
    set f2 [frame $t.f2 -relief groove -bd 2]
    set f3 [frame $t.f3 -relief groove -bd 2]
    set f4 [frame $t.f4 -relief groove -bd 2]
    set f5 [frame $t.f5]
    pack $f1 -side top -padx 7 -pady 7 \
	    -ipadx 7 -ipady 7 -fill x -expand 1
    pack $f2 -side top -padx 7 \
	    -ipadx 7 -ipady 7 -fill x -expand 1
    pack $f3 -side top -padx 7 -pady 7 \
	    -ipadx 7 -ipady 7 -fill x -expand 1
    pack $f4 -side top -padx 7 \
	    -ipadx 7 -ipady 7 -fill x -expand 1
    pack $f5 -side top -padx 7 -pady 7 \
	    -ipadx 7 -ipady 7 -fill x -expand 1
    # ---f1-f1-f1-f1-f1-f1-f1-f1---------------------------------
    set fr11 [frame $f1.1]
    set fr12 [frame $f1.2]
    set fr13 [frame $f1.3]
    pack $fr11 $fr12 $fr13 -side top -pady 5 -fill x
    
    set l1 [label $fr11.l1 -text "Coordinates of the centre of the cluster:"]
    pack $l1 -side top -expand 1
    
    Entries $fr11 "X: Y: Z:" "cutCL(X) cutCL(Y) cutCL(Z)" 7
    set cutclusfoclist " $fr11.frame.entry1 $fr11.frame.entry2 \
	    $fr11.frame.entry3"
    focus $fr11.frame.entry1
    
    set b [button $fr11.b -text "Select the Coordinates" \
	    -command [list SelClusCoor $t]]
    pack $b -expand 1
    
    set l1 [label $fr12.l1 -text "Maximum number of stars:"]
    set e1 [entry $fr12.e1 -relief sunken -width 3 -textvariable cutCL(NST)]
    set l2 [label $fr13.l2 -text "Radius of sphere centered at XYZ:"]
    set e2 [entry $fr13.e2 -relief sunken -width 5 -textvariable cutCL(RMAX)]
    pack $l1 $e1 $l2 $e2 -side left -expand 1
    append cutclusfoclist " $e1 $e2"
    
    # ---f2-f2-f2-f2-f2-------------------------------------------
    set cutCL(NNA) 0
    set anal No
    RadioButCmd $f2 "Print nearest neighbour\nanalysis of cluster atoms:" \
	    anal NeighAnal left left 0 1 2 "Yes" "No"
    # ---f3-f3-f3-f3-f3-------------------------------------------
    set cutCL(NCN) 0
    set coorn No
    RadioButCmd $f3 "Define user-defined\ncoordination numbers:" \
    	    coorn CoorNum left left 0 1 2 "Yes" "No"
    # ---f4-f4-f4-f4-f4-------------------------------------------
    set cutCL(NMO) 0
    set brdmod No
    RadioButCmd $f4 "Modify border atoms:" \
	    brdmod SetNMO left left 0 1 2 "Yes" "No"
    # ---f5-f5-f5-f5-f5-------------------------------------------
    set okf [frame $f5.okf -relief sunken -bd 2]
    set ok  [button $okf.ok -text "OK" \
	    -command [list CutClusOK $t]]
    set can [button $f5.can -text "Cancel" -command [list Cancel $t]]
    pack $okf -side left -expand 1 -padx 5m -pady 2m 
    pack $can -side left -expand 1 -padx 5m -pady 2m -ipadx 1m -ipady 1m
    pack $ok -padx 1m -pady 1m -ipadx 1m -ipady 1m
}


proc SelClusCoor {w} {
    global done

    if { [winfo exists .how2sel] } { return }
    set t [xcToplevel .how2sel "Choose" "Choose" . 581 74]
    set l1 [label $t.l1 -text "How to select centre of cluster:"]
    set b1 [button $t.b1 -text "Centre on Atom" \
	    -command [list CentreAtomCL $t]]
    set b2 [button $t.b2 -text "Centre on Hole" \
	    -command [list CentreHoleCL $t]]
    pack $l1 -side top -expand 1 -padx 10 -pady 10
    pack $b1 $b2 -side top -padx 10 -pady 5 -ipadx 5 -ipady 5 \
	    -expand 1 -fill x    

    tkwait visibility $t

    tkwait variable done
}


proc CentreAtomCL {w} {
    global done cutCL SelAtomCL

    destroy $w
    PreSel .catom .mesa "Coordinates of the Centre of the Cluster" \
	    "To select the centre of the cluster click on one atom" \
	    SelCentreAtomCL 1
    tkwait variable SelAtomCL(done)
    set cutCL(X) [lindex $SelAtomCL(centre) 0]
    set cutCL(Y) [lindex $SelAtomCL(centre) 1] 
    set cutCL(Z) [lindex $SelAtomCL(centre) 2]
    set done 1
}


proc CentreHoleCL {w} {
    global done cutCL SelHoleCL

    destroy $w
    PreSel .catom .mesa "Coordinates of the Centre of the Cluster" \
	    "A \"hole\" is geometrical centre of several atoms!!!\n\
	    To select the hole click on desired number of atoms" \
	    SelCentreHoleCL 15; # 15 is maximum allowed number of 
                                # selected atoms
    tkwait variable SelHoleCL(done)
    set cutCL(X) [lindex $SelHoleCL(centre) 0]
    set cutCL(Y) [lindex $SelHoleCL(centre) 1] 
    set cutCL(Z) [lindex $SelHoleCL(centre) 2]
    set done 1
}


proc CutClusOK {t} {
    global cutCL cutclusfoclist AdvGeom err system radio
    # t...          toplevel
    # brdmod...     flag for border modification
    
    # first we must check a variables
    check_var {{cutCL(X) real} {cutCL(Y) real} {cutCL(Z) real} \
	    {cutCL(NST) posint} {cutCL(RMAX) real}} $cutclusfoclist
    # check_var return $err; if $err --> mistake -> RETURN
    if $err { return }

    # if we come so far --> everything is OK
    # assign AdvGeom
    set n [xcAdvGeomState new]
    CutCLInput1

    destroy $t

    #
    # set cellmode to PRIMITIVE
    #
    set radio(cellmode) "prim"

    puts stdout "BORDER MODIFICATION: $cutCL(NMO)"
    puts stdout "CLUSTER"
    puts stdout $AdvGeom($n,cluster)
    
    # do we have border modification???
    if { $cutCL(NMO) == 0 } {
	# everything is ok --> cut a slab
	append AdvGeom($n,cluster) "$cutCL(NMO)\n"
	xcRenderCluster
    } else {
	# BORDER MODIFICATION::
	# first prepare for border modification --> do first run to get
        # sequential cluster number
	# for frist run set cutCL(NMO) to 0
	set cutCL(NMO) 0
	append AdvGeom($n,cluster) "$cutCL(NMO)\n"; #THIS IS TEMPORARY
	set firstRun [MakeInput]
	if { ! [RunC95 $system(c95_integrals) {} $firstRun] } {
	    #tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
	    return
	}
     
	puts stdout "CRYSTAL INPUT:\n $firstRun"
	puts stdout "CRYSTAL-95/98/03: OUTPUT"
	set output [ReadFile $system(SCRDIR)/xc_output.$system(PID)]
	# make a toplevel for DISPLAYING CRYSTAL OUTPUT
	set display [DispFirstRun $output]
	# make another toplevel for specifying BORDER MODIFICATION
	How2BordMod
    }
}


proc CutCLInput1 {} { 
    global AdvGeom cutCL
    # if we come so far --> everything is OK

    set n [xcAdvGeomState current]
    set AdvGeom($n,cluster) \
	    "$cutCL(X) $cutCL(Y) $cutCL(Z) $cutCL(NST) $cutCL(RMAX)\n"
    append AdvGeom($n,cluster) "$cutCL(NNA) $cutCL(NCN)\n"
    if { $cutCL(NNA) != 0 } {
	append AdvGeom($n,cluster) "$cutCL(RNNA)\n"
    }
    if { $cutCL(NCN) > 0 } {
	for {set i 1} {$i <= $cutCL(NCN)} {incr i} {
	    append AdvGeom($n,cluster) "$cutCL(L,$i) $cutCL(MCONN,$i)\n"
	}
    }
    return $n
}


proc NeighAnal {var} {
    global cutCL done
    
    # NeighAnal is procedure for neighbour nearest analysis for CLUSTER option
    if { $var == "No" } {
	# set cutCL(NNA) to zero
	set cutCL(NNA) 0
	return
    } elseif { $var == "Yes" } {
	set cutCL(NNA) 1

	if { [winfo exists .neigh] }  { return }
	set top [toplevel .neigh]
	xcPlace . .neigh 720 200
	wm title $top "Radius of Sphere"
	wm iconname $top "Radius"
    
	# some decorations
	set f1 [frame $top.f1 -relief raised -bd 2]
	set f2 [frame $top.f2 -relief raised -bd 2]
	pack $f1 $f2 -side top -padx 0 -pady 0 -ipadx 7 -ipady 7 -fill both
	
	set m [message $f1.msg -justify center -aspect 400 -text \
		"Radius of sphere in which to search for neighbours:"]
	pack $m -expand 1 -padx 5 -pady 5
	set e1 [entry $f1.ent -relief sunken -width 5 \
		-textvariable cutCL(RNNA)]
	pack $e1 -side top -expand 1
	focus $e1
	set ok [button $f2.ok -text "OK" -command [list NeighAnalOK $top $e1]]
	set can [button $f2.can -text "Cancel" \
		-command [list NeighAnalCan $top]]
	pack $ok $can -side left -padx 5 -pady 5 -ipadx 3 -ipady 3 -expand 1
	
	tkwait visibility $top
	catch { grab $top }
    
	tkwait variable done
	catch { grab release $top }

	destroy $top	
    }
}

	
proc NeighAnalOK {t e} {
    global cutCL err done
    
    #check if cutCL(RNNA) is real
    check_var {{cutCL(RNNA) real}} $e
    if $err {return}
    set done 1
}


proc NeighAnalCan {t} {
    global cutCL done anal
    
    set done 1
    set cutCL(NNA) 0 
    set anal "No"
}
	

proc CoorNum {var} {
    global cutCL nfoclist
    # NeighAnal is procedure for neighbour nearest analysis for CLUSTER option
    if { $var == "No" } {
	# set cutCL(NNA) to zero
	set cutCL(NCN) 0
	return
    } elseif { $var == "Yes" } {

	if [winfo exists .coorn] { return }
	set top [toplevel .coorn]
	xcPlace . .coorn 720 250
	wm title $top "Number"
	wm iconname $top "Number"
	
	# some decorations
	set f1 [frame $top.f1 -relief raised -bd 2]
	set f2 [frame $top.f2 -relief raised -bd 2]
	pack $f1 $f2 -side top -padx 0 -pady 0 -ipadx 7 -ipady 7 -fill both
	
	set m [message $f1.msg -justify center -aspect 400 -text \
		"Number of user-defined coordination numbers:"]
	pack $m -expand 1 -padx 5 -pady 5
	if { $cutCL(NCN) == 0 } {set cutCL(NCN) ""}
	set e1 [entry $f1.ent -relief sunken -width 5 \
		-textvariable cutCL(NCN)]
	focus $e1
	pack $e1 -side top -expand 1
	append nfoclist " $e1 "
	
	set ok [button $f2.ok -text "OK" -command [list NCoorOK $top $e1]]
	pack $ok -side top -padx 5 -pady 5 -ipadx 3 -ipady 3 -expand 1
	
	tkwait visibility $top
    }
}


proc NCoorOK {t e} {
    global nfoclist err cutCL
    
    # check if culCL(NCN) is positive integer
    check_var {{cutCL(NCN) point}} $nfoclist
    if $err {return}
    destroy $t
    GetCoorNum
}


proc GetCoorNum {} {
    global cutCL varlist foclist done
    
    if [winfo exists .getcoorn] { return }
    toplevel .getcoorn
    wm title .getcoorn "Custom Coordination Numbers"
    xcPlace . .getcoorn 580 190
    catch { grab .getcoorn }
    
    # bottom frame where OK button will be
    set fb [frame .getcoorn.fb]
    pack $fb -side bottom -expand true -fill both 
    # and one frame where canvas&scrollbar will be!!
    set ft [frame .getcoorn.ft -relief sunken -bd 2]
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
    for {set i 1} {$i <= $cutCL(NCN)} {incr i 1} {	    
	frame $f.fr$i -relief groove -bd 2
	pack $f.fr$i -padx 5 -pady 5
	label $f.fr${i}.label$i -text "Atom N.: $i" 
	pack $f.fr${i}.label$i -anchor w -padx 7 -pady 7
	frame $f.fr${i}.frm$i 
	pack $f.fr${i}.frm$i -side top -anchor w
	Entries $f.fr${i}.frm$i \
		{{Atomic Number:} {Coordination Number:}} \
		[list cutCL(L,$i) cutCL(MCONN,$i)] 4
	# make a varlist & foclist for PROC CHECK_VAR
	lappend varlist [list cutCL(L,$i) nat] [list cutCL(MCONN,$i) posint]
	lappend foclist \
		$f.fr${i}.frm$i.frame.entry1 $f.fr${i}.frm$i.frame.entry2
    }
	
    puts stdout "FOCLIST: $foclist\n\n"
    puts stdout "VARLIST: $varlist"
    set child [lindex [pack slaves $f] 0]
    
    # set the focus to first entry that upper FOR-LOOP create
    focus $f.fr1.frm1.frame.entry1
    
    tkwait visibility $child
    set width [winfo width $f]
    set height [winfo height $f]
    if { $cutCL(NCN) < 3 } {
	$c config -width $width -height $height 
    } else {
	$c config -width $width -height \
		[expr $height / $cutCL(NCN) * 3] \
		-scrollregion "0 0 $width $height"
    }
    
    button $fb.butok -text OK -command \
	    [list CoorNumOK .getcoorn]
    
    button $fb.butcan -text Cancel -command \
	    [list CoorNumCan .getcoorn]
    
    pack $fb.butok $fb.butcan -side left \
	    -expand 1 -padx 10 -pady 10
    
    tkwait variable done

    catch { grab release .getcoorn }
    destroy .getcoorn
}


proc CoorNumOK {t} {
    global cutCL done err varlist foclist
    # check if variables were correctly assigned
    check_var $varlist $foclist
    if $err {return}
    set done 1
}


proc CoorNumCan {t} {
    global cutCL done coorn
    
    set done 1
    set cutCL(NCN) 0 
    set coorn "No"
}


proc SetNMO {var} {
    global cutCL
    if { $var == "No" } {
	# set cutCL(NNA) to zero
	set cutCL(NMO) 0
    } else {
	set cutCL(NMO) 1
    }
}
	
	
proc DispFirstRun {output} {
    # procedure displays first crystal95 run --> to get cluster atoms
    # sequential numbers

    if [winfo exist .disp1run] {return}
    set t [toplevel .disp1run]
    wm title $t "Crystal95: Output"
    wm iconname $t "Crystal95: Output"
    xcPlace . .disp1run -50 70

    # in upper frame crystal output will be displayed
    set text [DispText $t.f1 $output 80 35]
    set f2 [frame $t.f2 -height 10]
    pack $f2 -side bottom -before $t.f1 -fill x
    set close [button $f2.cl -text "Close" -command [list destroy $t]]
    set view [button $f2.vi -text "View Cluster" -command ViewCluster]
    set num [button $f2.fi -text "Cluster Information" \
	    -command [list ClusInfo $text]]
    pack $close $view $num \
	    -side left -expand 1 -ipadx 2 -ipady 2 -pady 10
    return $t
}


proc ClusInfo {tw} {
    #tw .... textwidget 
    
    set index [$tw search -exact "CLUSTER CALCULATION" 1.0]
    $tw yview $index
}


proc ViewCluster {} {
    global cutCL AdvGeom
    CutCLInput1
    # temporary put zero for NMO
    set n [xcAdvGeomState current]
    append AdvGeom($n,cluster) "0\n"
    xcRenderCluster
}


proc How2BordMod {} {
    set t [xcToplevel .how2 "Modification of Border" "Border-Atoms" \
	    . 573 470]
    
    set l1 [label $t.l1 -text "Specify modification of border atoms:"]
    set b1 [button $t.b1 -text "separately to each atom" \
	    -command [list BordMod 1]]
    set b2 [button $t.b2 -text "to a group of atoms" \
	    -command [list BordMod 2]]
    pack $l1 -side top -expand 1 -padx 10 -pady 10
    pack $b1 $b2 -side top -padx 10 -pady 5 -ipadx 5 -ipady 5 \
	    -expand 1 -fill x
}


proc BordModSeq {} {
    global cutCL varlist foclist

    set varlist ""
    set foclist ""

    set t [xcToplevel .bordmod "Border Atoms" "Border Atoms" . 450 135]

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

    for {set i 1} {$i <= $cutCL(NMO)} {incr i 1} {	    
	frame $f.fr$i -relief groove -bd 2
	pack $f.fr$i -padx 5 -pady 5 -side top
	label $f.fr${i}.label$i -text "Atom N.: $i" 
	pack $f.fr${i}.label$i -anchor w -padx 7 -pady 7
	frame $f.fr${i}.frm$i 
	pack $f.fr${i}.frm$i -side top -anchor w
	Entries $f.fr${i}.frm$i \
		{"Label of atom:" "Number of stars:"} \
		[list cutCL(IPAD,$i) cutCL(NVIC,$i)] 4
	# make a varlist & foclist for PROC CHECK_VAR
	lappend varlist [list cutCL(IPAD,$i) posint] \
		[list cutCL(NVIC,$i) posint]
	lappend foclist \
		$f.fr${i}.frm$i.frame.entry1 $f.fr${i}.frm$i.frame.entry2 

	frame $f.fr$i.fram
	pack $f.fr${i}.fram -side top -anchor w
	if ![info exists cutCL(IPAR,$i)] [list set cutCL(IPAR,$i) 0]
	set chk [checkbutton $f.fr$i.fram.chk -text "Hydrogen Saturation" \
		-variable cutCL(IPAR,$i) ]
	set l21 [label $f.fr$i.fram.l1 -text "Hydrogen bond length:"]
	set e21 [entry $f.fr$i.fram.e1 -relief sunken -width 7 \
		-textvariable cutCL(BOND,$i)]
	pack $chk $l21 $e21 -side left -padx 5 -pady 5 -expand 1    
	# update a varlist & foclist for PROC CHECK_VAR
	lappend varlist [list cutCL(BOND,$i) real]
	lappend foclist $f.fr$i.fram.e1
    }
	
    puts stdout "FOCLIST: $foclist\n\n"
    puts stdout "VARLIST: $varlist"
    set child [lindex [pack slaves $f] 0]
    
    # set the focus to first entry that upper FOR-LOOP create
    focus $f.fr1.frm1.frame.entry1
    
    tkwait visibility $child
    set width [winfo width $f]
    set height [winfo height $f]
    if { $cutCL(NMO) < 3 } {
	$c config -width $width -height $height 
    } else {
	$c config -width $width -height \
		[expr $height / $cutCL(NMO) * 3] \
		-scrollregion "0 0 $width $height"
    }
    return $t
}


proc BordMod {var} {
    global cutCL ok varlist foclist

    destroy .how2
        if { $var == 1 } {
	# SEPARATELY TO EACH ATOM
	# HOW MANY ATOMS TO MODIFY
	set t [xcToplevel .brdatm "Modification of Border" \
	    "Border Atoms" . 576 482]
	set f1 [frame $t.f1 -relief raised -bd 2]
	set f2 [frame $t.f2 -relief raised -bd 2]
	pack $f1 $f2 -side top -padx 0 -pady 0 \
		-ipadx 10 -ipady 10 -fill both -expand 1
	set l11 [label $f1.l1 -text "Number of border atoms to be modified:"]
	set cutCL(NMO) ""
	set e11 [entry $f1.e1 -relief sunken -width 3 -textvariable cutCL(NMO)]
	focus $e11
	set b21 [button $f2.ok -text "OK" -command [list NBordM $e11]]
	pack $l11 $e11 -side top -expand 1
	pack $b21 -side top -ipadx 5 -ipady 5 -expand 1
	# wait unitl ok is set, then destroy $t
	tkwait variable ok
	destroy $t
	# specify MODIFICATIONS: for each atom separately
	set t [BordModSeq]; # BordModSeq returns a name for toplevel	
    } else {
	set t [xcToplevel .brdatm "Modification of Border" \
		"Border Atoms" . 340 380]
	set cutCL(NMO) -1
	set f1 [frame $t.f1 -relief groove -bd 2]
	pack $f1 -padx 5 -pady 5 -side top
	Entries $f1 \
		{"Label of first atom:" "Label of last atom:" \
		"Number of stars:"} \
		[list cutCL(IMIN) cutCL(IMAX) cutCL(NVIC,1)] 4
	focus $f1.frame.entry1
	# make a varlist & foclist for PROC CHECK_VAR
	set varlist ""
	lappend varlist [list cutCL(IMIN) posint] \
		[list cutCL(IMAX) posint] [list cutCL(NVIC,1) posint]
	set foclist "$f1.frame.entry1 \
		$f1.frame.entry2 $f1.frame.entry3"
	
	puts stdout "VARLIST: $varlist"
	puts stdout "FOCLIST: $foclist"
	frame $f1.fr
	pack $f1.fr -side top 
	set chk [checkbutton $f1.fr.chk -text "Hydrogen Saturation" \
		-variable cutCL(IPAR,1) ]
	set l21 [label $f1.fr.l1 -text "Hydrogen bond length:"]
	set e21 [entry $f1.fr.e1 -relief sunken -width 7 \
		-textvariable cutCL(BOND,1)]
	pack $chk $l21 $e21 -side left -padx 5 -pady 5 -expand 1    
	# update a varlist & foclist for PROC CHECK_VAR
	lappend varlist [list cutCL(BOND,1) real]
	lappend foclist $f1.fr.e1
    }

    set f3 [frame $t.f3]
    pack $f3 -side top -fill x -expand 1
    set okf [frame $f3.okf -relief sunken -bd 2]
    set ok  [button $okf.ok -text "OK" \
	    -command [list BordModOK $t]]
    set can [button $f3.can -text "Cancel" -command [list Cancel $t]]
    pack $okf -side left -expand 1 -padx 5m -pady 2m 
    pack $can -side left -expand 1 -padx 5m -pady 2m -ipadx 1m -ipady 1m
    pack $ok -padx 1m -pady 1m -ipadx 1m -ipady 1m
}
	
	
proc NBordM {foclist} {
    global cutCL err ok
    
    check_var {{cutCL(NMO) posint}} $foclist
    if $err {return}
    set ok 1
}


proc BordModOK {t} {
    global AdvGeom cutCL foclist varlist err

    # first we must check a variables
    if { $cutCL(NMO) == -1 } {
	# border modification of a group of atoms
	puts stdout "VARLIST: $varlist"
	puts stdout "FOCLIST: $foclist"
	check_var $varlist $foclist
	# check_var return $err; if $err --> mistake -> RETURN
	if $err { return }
	destroy $t
	# first part of CLUSTER INPUT
	set n [xcAdvGeomState current]
	CutCLInput1
	append AdvGeom($n,cluster) "$cutCL(NMO)\n"
	append AdvGeom($n,cluster) "$cutCL(IMIN) $cutCL(IMAX) $cutCL(NVIC,1) \
		$cutCL(IPAR,1) $cutCL(BOND,1)\n"
	puts stdout "\n\nCRYSTAL INPUT::"
	puts stdout $AdvGeom($n,cluster)
    } elseif { $cutCL(NMO) > 0 } {
	check_var $varlist $foclist
	if $err {return}
	destroy $t
	CutCLInput1
	append AdvGeom($n,cluster) "$cutCL(NMO)\n"	
	for {set i 1} {$i <= $cutCL(NMO)} {incr i} {
	    append AdvGeom($n,cluster) "$cutCL(IPAD,$i) $cutCL(NVIC,$i) \
		    $cutCL(IPAR,$i) $cutCL(BOND,$i)\n"
	}
    }
    # calculate CLUSTER
    xcRenderCluster
}


proc xcRenderCluster {} {
    global system n_groupsel err
    
    # this  is for UNDO/REDO
    GenCommUndoRedo "Cut a Cluster"

    CalStru
    
    # update structure
    GenGeomDisplay 1
    
    #set input [MakeInput]
    #puts stdout "CRYSTAL INPUT:\n $input"	
    #if { ! [RunC95 $system(c95_integrals) {} $input] } {
    #	#tk_dialog .err "ERROR" "ERROR: $err" error 0 OK
    #	return
    #}
    #
    ## this  is for UNDO/REDO
    #GenCommUndoRedo "Cut a Cluster"
    #
    ## update structure
    #GenGeomDisplay 1
}
