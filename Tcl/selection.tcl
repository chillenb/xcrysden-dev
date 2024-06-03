#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/selection.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# warning:: some outer procs use select global variable as well, so 
# You should never unset the select variable after done/selected/cancel 
# button was closed
#
proc SelCheck {com oper} {
    global maxsel nsel
    #com .... command to exexute if OK
    #oper ... name of command/operation to execute if OK
    
    if { $maxsel != $nsel } {
	dialog .selcheck Warning "WARNING: You must select $maxsel \
		atoms for $oper" warning 0 OK
	return
    } else {
	eval $com
    }
}


proc Done {w com {tplw {}}} {
    global maxsel nsel select SelPlane SelAtomCL SelHoleCL SelLabel SelLine
    
    # w...name of canvas
    # tlpw...name of selection toplevel

    if { $nsel == $maxsel } {
	for {set i 1} {$i <= $maxsel} {incr i} {
	    selectRecordToFields_ $i
	}
	# this 'ifs' are for PLANE SELECTION -- AdvOption

	puts stderr "Selection Done: com == $com"

	if { [lindex $com 0] != "PlaneSel" && \
		[lindex $com 0] != "PlaneSelSel" && \
		[lindex $com 0] != "SelCentreAtomCLSel" && \
		[lindex $com 0] != "LabelSelSel" && \
		[lindex $com 0] != "LineSelSel" && \
		[lindex $com 0] != "ParalleSel" && \
		[lindex $com 0] != "AtomSel" } {	    
	    eval $com
	    return
	} elseif { [lindex $com 0] == "PlaneSel" } {
	    Plane
	    return
	} elseif { [lindex $com 0] == "SelCentreAtomCL" } {
	    SelCentreAtomCL
	    return
	} elseif { [lindex $com 0] == "PlaneSelSel" } {
	    set SelPlane(plane) [Plane]
	    DeSel $tplw $w
	    set SelPlane(done) 1
	    return $SelPlane(plane)
	} elseif { [lindex $com 0] == "SelCentreAtomCLSel" } {
	    set SelAtomCL(centre) [SelCentreAtomCL]
	    DeSel $tplw $w
	    set SelAtomCL(done) 1
	    return $SelAtomCL(centre)
	} elseif { [lindex $com 0] == "LabelSelSel" } {
	    set SelLabel(label) [LabelSel]
	    set SelLabel(done) 1
	    DeSel $tplw $w
	    return $SelLabel(label)
	} elseif { [lindex $com 0] == "LineSelSel" } {
	    set SelLine(coor) [LineSel]
	    set SelLine(done) 1
	    DeSel $tplw $w
	    return $SelLine(coor)
	} elseif { [lindex $com 0] == "ParalleSel" || \
		[lindex $com 0] == "AtomSel"} {
	    DeSel $tplw $w
	    set select(done) 1
	    return
	}
    } else {
	# for SelCentreHoleCLSel its not needed that $nsel == $maxsel
	for {set i 1} {$i <= $nsel} {incr i} {
	    selectRecordToFields_ $i
	    #set select(Nat$i) [lindex $select(obj,$i) 2]
	    #set select(X$i) [lindex $select(obj,$i) 3]
	    #set select(Y$i) [lindex $select(obj,$i) 4]
	    #set select(Z$i) [lindex $select(obj,$i) 5]
	}
	if { [lindex $com 0] == "SelCentreHoleCL" } {
	    SelCentreHoleCL
	    return 
	} elseif { [lindex $com 0] == "SelCentreHoleCLSel" } {
	    set SelHoleCL(done) 1
	    set SelHoleCL(centre) [SelCentreHoleCL]
	    DeSel $tplw $w
	    return $SelHoleCL(centre)
	}
    }
}


proc NextSel {com} {
    global select maxsel nsel

    xc_select .mesa clean
    xc_select .mesa begin
    set select(initialized) 1 

    # delete all text
    $select(textWid) config -state normal
    $select(textWid) delete 0.0 end 

    # update "iron text"
    set nsel 0
    selectUpdateText $com
    #VoidTextSel $select(textWid) $com $maxsel
}


#-----------------------------------------------------------------------------
# this is for Line-Adding type of ATOMINSE command
proc LineSel {{w {}}} {
    global select maxsel SelLine err system
    
    ########################################
    # CD to $system(SCRDIR)
    cd $system(SCRDIR)
    ########################################

    # is fraction of line-length specified correctly
    check_var {{SelLine(fract) fract}} $SelLine(entry)
    if $err { return }
    
    # R = R1 + t*(R2-R1); t...fraction of R2-R1
    set t $SelLine(fract)
    set x [expr $select(X1) + $t*($select(X2) - $select(X1))] 
    set y [expr $select(Y1) + $t*($select(Y2) - $select(Y1))] 
    set z [expr $select(Z1) + $t*($select(Z2) - $select(Z1))] 

    # translate point to first cell
    set fileID [open "$system(SCRDIR)/xc_tmp.$system(PID)" w]
    puts $fileID "0 $x $y $z"
    flush $fileID
    close $fileID
    
    $select(textWid) config -state normal
    
    # t.k.: atomlab return coordinates in Angs units

    if { [catch {set coor [exec $system(FORDIR)/atomlab 2 \
			       $system(SCRDIR)/xc_struc.$system(PID) \
			       $system(SCRDIR)/xc_tmp.$system(PID)]} errmsg] } {
	$select(textWid) delete [expr $maxsel + 4].0 [expr $maxsel + 5].end
	$select(textWid) insert [expr $maxsel + 5].0 $errmsg
	$select(textWid) config -state disabled
    } else {
	set xx [lindex $coor 0]
	set yy [lindex $coor 1]
	set zz [lindex $coor 2]

	set coor [coorToUnit $select(unit)  $x $y $z]
	set tr_coor [coorToUnit $select(unit)  $xx $yy $zz]

	$select(textWid) insert [expr $maxsel + 4].0 "Coordinates of point are: \
		[eval {format {%6.10f  %6.10f  %6.10f}} $coor]\n"

	if { abs($x-$xx) > 1e-5 || abs($y-$yy) > 1e-5 || abs($z-$zz) > 1e-5 } {
	    $select(textWid) insert [expr $maxsel + 5].0 "Coordinates translated to: \
		[eval {format {%6.10f  %6.10f  %6.10f}} $tr_coor]\n"	    
	}

	$select(textWid) config -state disabled
	return "$xx $yy $zz"
    }
}


# -------------------------------------------------------------------------
# this is for selecting a centre of CLUSTER (AdvGeom option)
# WE NEED TO SELECT ONE ATOM
proc SelCentreAtomCL {{w {}}} {
    global select maxsel

    set coor [coorToUnit $select(unit)  $select(X1) $select(Y1) $select(Z1)]

    $select(textWid) config -state normal
    $select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end
    $select(textWid) insert [expr $maxsel + 5].0 "Coordinates of selected atom are: $coor\n"
    $select(textWid) config -state disabled

    return "$select(X1) $select(Y1) $select(Z1)"
}


# -------------------------------------------------------------------------
# this is for selecting a centre of CLUSTER (AdvGeom option)
# WE NEED TO SELECT A HOLE (geometrical centre of several atoms)
proc SelCentreHoleCL {{w {}}} {
    global select nsel maxsel SelHoleCL system
    #nsel ... number of selected atoms

    ########################################
    # CD to $system(SCRDIR)
    cd $system(SCRDIR)
    ########################################
    
    set xc 0.0
    set yc 0.0
    set zc 0.0
    for {set i 1} {$i <= $nsel} {incr i} {
	set xc [expr $xc + $select(X$i)]
	set yc [expr $yc + $select(Y$i)]
	set zc [expr $zc + $select(Z$i)]
    }
    set xc [expr $xc / $nsel]
    set yc [expr $yc / $nsel]
    set zc [expr $zc / $nsel]
 
    $select(textWid) config -state normal

    $select(textWid) delete [expr $nsel + 3].0 [expr $maxsel + 5].end
    $select(textWid) insert [expr $nsel + 3].0 \
	"\n----------------------------------------------------------------------------------\n\n"

    # if SetHoleCL(transl) exists --> translate hole to basic cell

    if { [info exists SelHoleCL(transl)] } {

	WriteFile "$system(SCRDIR)/xc_tmp.$system(PID)"  "0 $xc $yc $zc"
	
	if { [catch {set label [exec $system(FORDIR)/atomlab 2 \
		$system(SCRDIR)/xc_struc.$system(PID) \
		$system(SCRDIR)/xc_tmp.$system(PID)]} errmsg] } {
	    $select(textWid) delete [expr $maxsel + 4].0 [expr $maxsel + 5].end
	    $select(textWid) insert [expr $maxsel + 5].0 $errmsg

	    $select(textWid) config -state disabled
	} else {
	    # todo: if the point wasn't translated, don't make the
	    # message: "coordinates translated to ..."
	    set x [lindex $label 0]
	    set y [lindex $label 1]
	    set z [lindex $label 2]


	    set coor [coorToUnit $select(unit)  $xc $yc $zc]
	    set tr_coor [coorToUnit $select(unit)  $x $y $z]

	    $select(textWid) insert [expr $maxsel + 4].0 "Coordinates of \"hole\" are: \
		    [eval {format {%6.10f  %6.10f  %6.10f}} $coor]\n"

	    if { abs($x-$xc) > 1e-5 || abs($y-$yc) > 1e-5 || abs($z-$zc) > 1e-5 } {
		$select(textWid) insert [expr $maxsel + 5].0 "Coordinates translated to: \
		    [eval {format {%6.10f  %6.10f  %6.10f}} $tr_coor]\n"
	    }
	    
	    $select(textWid) config -state disabled
	    return "$x $y $z"
	}
    } else {
	
	set coor [coorToUnit $select(unit)  $xc $yc $zc]

	$select(textWid) insert [expr $nsel + 5].0 \
	    "Coordinates of \"hole\" are: [eval {format {%6.10f  %6.10f  %6.10f}} $coor]\n"

	$select(textWid) config -state disabled
	return "$xc $yc $zc"
    }
}
    
    
proc Plane {{w {}}} {
    global select maxsel Plane system


    ########################################
    # CD to $system(SCRDIR)
    cd $system(SCRDIR)
    ########################################

    #----------------------------------------------
    # I THINK THAT $w is NOT NEEDED -CHECK THIS
    # if I gonna remove $w, than it should be also removed from the
    # Angle Dihedral Distance procs & from all calls
    #-----------------------------------------------

    # first check for linear dependence of atomic positions
    set x21 [expr $select(X1) - $select(X2)]
    set y21 [expr $select(Y1) - $select(Y2)]
    set z21 [expr $select(Z1) - $select(Z2)]
    set d21 [expr sqrt( $x21 * $x21 + $y21 * $y21 + $z21 * $z21)] 

    set x32 [expr $select(X2) - $select(X3)]
    set y32 [expr $select(Y2) - $select(Y3)]
    set z32 [expr $select(Z2) - $select(Z3)]
    set d32 [expr sqrt( $x32 * $x32 + $y32 * $y32 + $z32 * $z32)] 

    set ab [expr $x21 * $x32 + $y21 * $y32 + $z21 * $z32]
    set dif [ expr abs($ab - $d21 * $d32) ]
    
    $select(textWid) config -state normal

    if { $dif < 0.01 } {
	$select(textWid) delete [expr $maxsel + 4].0 [expr $maxsel + 5].end
	$select(textWid) insert [expr $maxsel + 5].0 \
		"WARNING: linear dependence of atomic positions !!!"
	return
    }
     
    # puts the coor of 3 sel. points into $system(SCRDIR)/xc_tmp.$pid
    set plncrd [open "$system(SCRDIR)/xc_tmp.$system(PID)" w]
    puts $plncrd "$select(X1) $select(Y1) $select(Z1)"
    puts $plncrd "$select(X2) $select(Y2) $select(Z2)"
    puts $plncrd "$select(X3) $select(Y3) $select(Z3)"
    flush $plncrd
    close $plncrd
    puts stdout "$select(X1) $select(Y1) $select(Z1)"
    puts stdout "$select(X2) $select(Y2) $select(Z2)"
    puts stdout "$select(X3) $select(Y3) $select(Z3)\n##############"
    #eval [list exec cat $system(SCRDIR)/xc_tmp.$system(PID)]
    ReadFile $system(SCRDIR)/xc_tmp.$system(PID)

    # now calculate a plane
    xcDebug "$system(FORDIR)/calplane \
	    $system(SCRDIR)/xc_struc.$system(PID) \
	    $system(SCRDIR)/xc_tmp.$system(PID)"

    set fileID [open "|$system(FORDIR)/calplane \
	    $system(SCRDIR)/xc_struc.$system(PID) \
	    $system(SCRDIR)/xc_tmp.$system(PID)"]
    set Plane [read $fileID]
    set Plane1 [lrange $Plane 0 2]
    set Plane2 [lrange $Plane 3 5]
    close $fileID

    $select(textWid) delete [expr $maxsel + 4].0 [expr $maxsel + 5].end
    $select(textWid) insert [expr $maxsel + 4].0 \
	    "Unrounded Miller indexes: $Plane1\n"
    $select(textWid) insert [expr $maxsel + 5].0 \
	    "Miller indexes rounded to: $Plane2"
    
    $select(textWid) config -state disabled

#    $select(textWid) insert [expr $maxsel + 3].0 "$Plane"
    return $Plane2
}    


proc Dihedral {w} {
    global select maxsel

    set pi 3.14159265358979323844
    # we've three vectors: 21,32,34
    # vector 32 is normal vector
    set x21 [expr $select(X1) - $select(X2)]
    set y21 [expr $select(Y1) - $select(Y2)]
    set z21 [expr $select(Z1) - $select(Z2)]

    # normal vector:
    set x32 [expr $select(X2) - $select(X3)]
    set y32 [expr $select(Y2) - $select(Y3)]
    set z32 [expr $select(Z2) - $select(Z3)]

    set x34 [expr $select(X4) - $select(X3)]
    set y34 [expr $select(Y4) - $select(Y3)]
    set z34 [expr $select(Z4) - $select(Z3)]

    # now vectror ar...projection of 21 on surface R (n is norm.v. of R)
    #             br...projection of 34 on surface R
    # la -lambda a; factor
    # lb -lambda b; factor
    set sd32 [expr $x32 * $x32 + $y32 * $y32 + $z32 * $z32]
    set la [expr ( $x21 * $x32 + $y21 * $y32 + $z21 * $z32 ) / $sd32 ]
    set lb [expr ( $x34 * $x32 + $y34 * $y32 + $z34 * $z32 ) / $sd32 ]

    set xar [expr $x21 - $la * $x32]
    set yar [expr $y21 - $la * $y32]
    set zar [expr $z21 - $la * $z32]
    set ar [expr sqrt($xar * $xar + $yar * $yar + $zar * $zar)]

    set xbr [expr $x34 - $lb * $x32]
    set ybr [expr $y34 - $lb * $y32]
    set zbr [expr $z34 - $lb * $z32]
    set br [expr sqrt($xbr * $xbr + $ybr * $ybr + $zbr * $zbr)]

    set arbr [expr $ar * $br]

    $select(textWid) config -state normal

    if { $arbr < 0.01 } {
	$select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end
	$select(textWid) insert [expr $maxsel + 5].0 \
		"WARNING: Linear dependance of atomic positions !!!"
	return
    }
    set var [expr ( $xar * $xbr + $yar * $ybr + $zar * $zbr ) / \
	    ( $arbr )]
    # $var colud be greater than one, due to roundoff error
    if { $var > 1.0 } {set var 1.0}
    if { $var < -1.0 } {set var -1.0}
    puts stdout "acos:: $var"
    set dihedr [expr acos($var) * 180.0 / $pi]
    $select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end 
    $select(textWid) insert [expr $maxsel + 5].0 [format "Dihedral Angle is %.${select(angl_precision)}f deg" $dihedr]

    $select(textWid) config -state disabled
    return $dihedr
}
    
    
proc Angle {w} {
    global select maxsel nsel

    set pi 3.14159265358979323844
    # this is just for testing
    if { $nsel == $maxsel } {	
	set x21 [expr $select(X1) - $select(X2)]
	set y21 [expr $select(Y1) - $select(Y2)]
	set z21 [expr $select(Z1) - $select(Z2)]
	set dist21 [expr sqrt($x21 * $x21 + $y21 * $y21 + $z21 * $z21)] 
	set x23 [expr $select(X3) - $select(X2)]
	set y23 [expr $select(Y3) - $select(Y2)]
	set z23 [expr $select(Z3) - $select(Z2)]
	set dist23 [expr sqrt($x23 * $x23 + $y23 * $y23 + $z23 * $z23)]
	set var [expr ( $x21 * $x23 + $y21 * $y23 +  $z21 * $z23 ) / \
		( $dist21 * $dist23 )]
	# var could be greater than 1.0, due to roundoff error
	if { $var > 1.0 } {set var 1.0}
	if { $var < -1.0 } {set var -1.0}
	set angl [expr acos($var) * 180.0 / $pi]
	
	$select(textWid) config -state normal

	$select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end 
	$select(textWid) insert [expr $maxsel + 5].0 [format "Angle is %.${select(angl_precision)}f deg" $angl]

	$select(textWid) config -state disabled   
	return $angl
    }
}


proc Distance {w} {
    global select maxsel nsel Const

    # this is just for testing
    if { $nsel == $maxsel } {
	puts stdout "DISTANCE BETWEEN ATOMS: \n\
		$select(X1) $select(Y1) $select(Z1) \n\
		$select(X2) $select(Y2) $select(Z2)"
	set dx [expr $select(X1) - $select(X2)]
	set dy [expr $select(Y1) - $select(Y2)]
	set dz [expr $select(Z1) - $select(Z2)]
	set dist [expr sqrt($dx * $dx + $dy * $dy + $dz * $dz)] 
	set distBohr [expr $dist / $Const(bohr)]
	# first delete, then insert

	$select(textWid) config -state normal

	$select(textWid) delete [expr $maxsel + 5].0 [expr $maxsel + 5].end 
	$select(textWid) insert [expr $maxsel + 5].0 [format "Distance is %.${select(dist_precision)}f ANGSTROMS = %.${select(dist_precision)}f Bohrs"    $dist $distBohr]

	$select(textWid) config -state disabled

	return $dist
    }
}
    

proc AtomInfo {{w {}}} {
    # this is dummy procedure
    return 
}


proc SelectAtom { com x y } {
    global select nsel maxsel Xor Yor df zoom

    xcDebug -debug "SelectAtom"

    # if we click the atom first time --> SELECTION
    # if we click the atom second time --> DESELECTION

    set sel 1
    set obj [xc_select .mesa sqn $x $y]
    # if $obj == 0 -> no atom was selected
    if { $obj == 0 } { 
	xcDebug -debug "AtomSelect End"
	return 
    }

    $select(textWid) config -state normal

    for {set i 1} {$i <= $nsel} {incr i} {

	if { $obj == $select(atomID,$i) } {

	    # DESELECT ATOM

	    xc_deselect .mesa atom $obj
	    
	    set sel 0
	    for {set j $i} {$j < $nsel} {incr j} {
		set jj [expr $j + 1]
		set select(atomID,$j) $select(atomID,$jj)
		set select(obj,$j) $select(obj,$jj)
	    }
	    
	    incr nsel -1		   		
	    selectUpdateText $com

	    # this is for selection bond-lines rearranging

	    if { $nsel >= 1 } { xc_deselect .mesa line }
	    for {set j 1} {$j < $nsel} {incr j} {
		set jj [expr $j + 1]		
		set select(lineID,$j) [xc_select .mesa line \
					   $select(atomID,$j) $select(atomID,$jj)]
	    }
	    # finish changes
	    xc_select .mesa finish
	    return
	}
    }
    	
    if { $nsel < $maxsel && $sel == 1} {
	incr nsel
	set select(atomID,$nsel) $obj	
	set select(obj,$nsel) [xc_select .mesa atom $obj]
	selectUpdateText $com
    }

    # draw a line between select atom & previous selected atom

    if { $nsel <= $maxsel && $nsel > 1 } {
	set nsel1 [expr $nsel - 1]
	set select(lineID,$nsel1) [xc_select .mesa line \
		$select(atomID,$nsel1) $select(atomID,$nsel)]
    }

    $select(textWid) config -state disabled

    # now flush to display all changes
    xc_select .mesa finish
    xcDebug -debug "AtomSelect End"
    update
}

proc PreSel {w c topname title com m} {
    global maxsel nsel SelLine select check periodic

    # PreSel ...PrepareSelection
    # w ...     name of topwidget
    # c ...     name of canvas
    # topname . name of toplevel
    # title ... title 
    # com ..... command to exexute
    # m ....... set maxsel $m

    if { $check(perspective) } {
	set select(revert_to_perspective) 1
	set check(perspective) 0
	Perspective
    } else {
	set select(revert_to_perspective) 0
    }

    set select(selection_mode) 1

    # some initializations
    
    $c config -cursor cross
    bind $c <Button-1> [list SelectAtom $com %x %y]

    # due to configuring of $c the display of $c vanish, 
    # to fix that call "update"
    update

    xc_select .mesa begin
    set select(initialized)    1
    set nsel 0
    set maxsel $m
    
    # make selection toplevel widget

    if { [winfo exists $w] } { return }
    xcToplevel $w $topname "Selection" . 0 0 1
    AlwaysOnTopON . $w
    bind $w <Destroy> [list DeSel $w $c]
    focus $w

    set topfrm [frame $w.topfrm -class StressText]
    set titl [label $topfrm.lbl -text $title -relief groove -bd 2]
    pack $topfrm -side top -expand 1 -fill x 
    pack $titl -expand 1 -fill x -ipady 3 -ipadx 10 -pady 10 -padx 10
    #set frm [frame $w.frm -bd 1]

    if { $maxsel <= 4 } {
	set height [expr $maxsel + 5]
	set select(textWid) [text $w.text -relief sunken \
		-width 82 -height $height -font TkFixedFont]
    } else {
	# for SelCentreHoleCL maxsel could be up to 15
	set height 12
	set fts [frame $w.f1]
	pack $fts -side top -expand 1
	set select(textWid) [text $fts.text -relief sunken -width 82 -height \
		$height -yscrollcommand "$fts.sy set" -font TkFixedFont]
	scrollbar $fts.sy -orient vert -command "$select(textWid) yview" 
	pack $fts.sy -side right -fill y
    }
    
    # lengt units radiobuttons

    if { ! [info exists select(unit)] } {
	set select(unit) angs
    }
    set midfrm [labelframe $w.midfrm -text "Display coordinates in units:"]
    foreach unit {
	angs bohr conv prim alat
    } txt {
	Angstrom Bohr Crystal-Convetional Crystal-Primitive Alat
    } {
	radiobutton $midfrm.$unit -text $txt -variable select(unit) -value $unit -anchor w
	pack $midfrm.$unit -side left -padx 2 -pady 0 -fill x -expand 1
    }
    
    if { $periodic(dim) == 0 } {
	# disable conv prim alat radiobuttons
	foreach unit {conv prim alat} {
	    $midfrm.$unit configure -state disabled
	}
    }

    _select_unit $com select unit ""
    trace add variable select(unit) write [list _select_unit $com]

    # buttons "Done", "Next", "Close"

    set botfrm [frame $w.botfrm]

    # if $comm == "LineSel" -> there must also be entry for specifying the 
    # line-length fraction
    if { $com == "LineSel" } {
	set llf [frame $botfrm.llf -relief groove -bd 2]
	pack $llf -side left -expand 1 -padx 20 -pady 10
	Entries $llf {{Line-Length Fraction:}} SelLine(fract) 10
	# remember which entry to focus in case of error
	set SelLine(entry) $llf.frame.entry1
    }
    
    if { $com != "AtomInfo" && $com != "AtomSel" && $com != "ParalleSel" } {
	set but1 [button $botfrm.done -text "Done" -command \
		[list Done $c "$com $c"]]
	pack $but1 -side left -padx 5 -pady 3 -expand 1
    }
    
    set but2 [button $botfrm.next -text "Next" \
	    -command [list NextSel $com]]
    set but3 [button $botfrm.can -text "Close" -command \
	    [list DeSel $w $c]]

    pack $select(textWid) -side top -padx 10 -pady 3 -fill x -expand 1    
    pack $midfrm -side top -padx 10 -fill x -expand 1
    pack $botfrm -side top -fill both -expand 1
    pack $but2 $but3 -side left -padx 5 -pady 3 -expand 1

    # extra button when we are selecting a Plane or "Centre of
    # Cluster" for AdvGeom
    if { $com == "PlaneSel" || $com == "SelCentreAtomCL" || \
	    $com == "SelCentreHoleCL" || $com == "LabelSel" || \
	    $com == "LineSel" } {
	set but4 [button $botfrm.sel -text "Selected" -command \
		[list Done $c "${com}Sel $c" $w]]
	pack $but4 -side left -padx 5 -pady 3 -expand 1
    }
    
    # extra button when we are selecting the IsoSpace/IsoPlane
    if { $com == "AtomSel" || $com == "ParalleSel" } {
	set but4 [button $botfrm.sel -text "Selected" -command \
		[list Done $c ${com} $w]]
	pack $but4 -side left -padx 5 -pady 3 -expand 1
    }
    
    selectUpdateText $com
    #VoidTextSel $select(textWid) $com $maxsel 
    return $w
}


#proc VoidTextSel {textWid com maxsel} {
#    $textWid config -state normal
#    $textWid insert 1.0 \
#	"                     ID Sym Atm.Num  X/ANGSTROMS    Y/ANGSTROMS    Z/ANGSTROMS\n"
#    $textWid insert 2.0 \
#	"--------------------------------------------------------------------------------\n"
#    for {set i 1} {$i <= $maxsel} {incr i} {
#	if { $com != "SelCentreHoleCL" } {
#	    $textWid insert [expr $i + 2].0 "Selected Atom No.$i: \n"
#	}
#    }
#    # if $com == "SelCentreHoleCL", we do not know how many atoms will be 
#    # selected, but for all other $com we do know
#    if { $com != "SelCentreHoleCL" } {
#	$textWid insert [expr $maxsel + 3].0 \
#	    "--------------------------------------------------------------------------------\n\n"
#    }
#    
#    $textWid config -state disabled
#}

proc DeSel {w c} {
    global xcCursor select check
    
    xcDebug "DeSel"
    if { $select(initialized) } {
	# in order to clean, initialization should be performed first !!!
	xc_select .mesa clean
    }
    set select(selection_mode) 0
    set select(initialized)    0
    
    $c config -cursor $xcCursor(default)
    bind $c <Button-1> {}

    if { [winfo exists $w] } { 
	bind $w <Destroy> {}
	destroy $w 
    }
    
    if { $select(revert_to_perspective) } {
	set check(perspective) 1
	Perspective
    }
}


proc selectRecordToFields_ {i} {
    global select

    set select(Sqn$i) [lindex $select(obj,$i) 0]
    set select(Sym$i) [lindex $select(obj,$i) 1]
    set select(Nat$i) [lindex $select(obj,$i) 2]
    set select(X$i)   [lindex $select(obj,$i) 3]
    set select(Y$i)   [lindex $select(obj,$i) 4]
    set select(Z$i)   [lindex $select(obj,$i) 5]
}


proc selectDisplayCoor_ {i} {
    global select

    selectRecordToFields_ $i
    
    set coor [coorToUnit $select(unit)  $select(X$i) $select(Y$i) $select(Z$i)]
    
    set select(disp_X$i) [lindex $coor 0]
    set select(disp_Y$i) [lindex $coor 1] 
    set select(disp_Z$i) [lindex $coor 2]
}


proc selectUpdateText {{com ""}} {
    global select maxsel nsel

    switch -- $select(unit) {
	bohr { set unit Bohr }
	prim { set unit CrystalPrim }
	conv { set unit CrystalConv }
	alat { set unit Alat }
	angs - default { set unit Angstrom }
    }

    for {set i 1} {$i <= $nsel} {incr i} {
	selectDisplayCoor_ $i
    }

    $select(textWid) config -state normal
    $select(textWid) delete 0.0 end 
    $select(textWid) insert 1.0 \
	[format "%21s ID Sym Atm.Num  X/%-12s Y/%-12s Z/%-9s\n" { } $unit $unit $unit]
    $select(textWid) insert 2.0 \
	"----------------------------------------------------------------------------------\n"

    for {set i 1} {$i <= $nsel} {incr i} {
	set data [format " %-3d %-3s %-3d  %+14.9f %+14.9f %+14.9f" \
		      $select(Sqn$i) $select(Sym$i) $select(Nat$i) $select(disp_X$i) $select(disp_Y$i) $select(disp_Z$i)]  
	$select(textWid) insert [expr $i + 2].0 "Selected Atom No.$i: $data\n"
    }
    for {set i [expr $nsel + 1]} {$i <= $maxsel} {incr i} {
	if { $com != "SelCentreHoleCL" } {
	    $select(textWid) insert [expr $i + 2].0 "Selected Atom No.$i: \n"
	}
    }

    # if $com == "SelCentreHoleCL", we do not know how many atoms will be 
    # selected, but for all other $com we do know
    if { $com != "SelCentreHoleCL" } {
	$select(textWid) insert [expr $maxsel + 3].0 \
	    "----------------------------------------------------------------------------------\n\n"
    }
    
    $select(textWid) config -state disabled
}


proc _select_unit {com name1 name2 ops} {
    global select
    
    if { $name1 != "select" && $name2 != "unit" } {
	return
    }
    selectUpdateText $com
}
