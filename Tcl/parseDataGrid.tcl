#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/parseDataGrid.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc parseDGInfo info {
    global sInfo DG system

    #
    # parse datagrid INFO record
    #
    set DG(n_block)  [lindex $info 0]

    for {set i 1} {$i <= $DG(n_block)} {incr i} {
	set im  [expr $i - 1]
	set block  [lindex $info $i]
	set DG(type,$im)  [lindex $block 0]
	set DG(ident,$im) [lindex $block 1]
	set DG(n_subblock,$im)  [lindex $block 2]
	set subblock  [lindex $block 3]
	xcDebug "   Block Number: $i"
	xcDebug "                type: $DG(type,$im)"
	xcDebug "               ident: $DG(ident,$im)"
	xcDebug "      # of subblocks: $DG(n_subblock,$im)"
	xcDebug "            subblock: $subblock"
	for {set j 0} {$j < $DG(n_subblock,$im)} {incr j} {
	    set DG(subident,$im,$j) [lindex $subblock $j]
	}
    }
    
    set t [xcToplevel .dg "DataGrid: data presentation" \
	    "DataGrid" . 50 0 1]
    set f1 [frame $t.1]
    set f2 [frame $t.2]
    set c  [canvas $f1.c -yscrollcommand [list $f1.sy set] \
		-xscrollcommand [list $f2.sx set] \
		-width 550 -height 300 \
		-bg "#ffffff"]
    set sy [scrollbar $f1.sy -orient vertical \
	    -command [list $f1.c yview]]
    set sx [scrollbar $f2.sx -orient horizontal \
	    -command [list $f1.c xview]]
    set pad [expr [$sy cget -width] + \
	    2 * ([$sy cget -bd] + [$sy cget -highlightthickness])]
    set fpad [frame $f2.pad -width $pad -height $pad]

    pack $f1 -side top -fill both -expand 1
    pack $f2 -side top -fill x
    pack $c -side left -fill both -expand true
    pack $sy -side left -fill y
    pack $fpad -side right   
    pack $sx -side bottom -fill x

    button .__b
    set font [.__b cget -font]
    set datagrid_font [ModifyFont $font .__b \
	    -family courier -size 20 -weight bold -default 1]
    destroy .__b

    # image dimension: 109x39
    image create photo datagrid -format gif \
    	    -file $system(BMPDIR)/datagrid.gif
	
    set x0  30; # x position of main vertical line
    set y0  30; # y center position of datagrid image 
    set x1  60; # x center position of datagrdi image
    set xvl 70; # x position of subblock vertical line

    $c create line $x0 $y0 $x1 $y0 \
	    -width 4 -tags mainLine0 -arrow last \
	    -tags mainLine0 -fill #000
    $c create image $x1 $y0 -image datagrid -anchor center -tags image

    if { ![info exists DG(radio)] } {
	set DG(radio) 0
	set DG(cb0,0) 1
    }
    set Ys $y0
    for {set i 0} {$i < $DG(n_block)} {incr i} {
	DG_CreateBlock $c $x1 $Ys $i
	if {$i == 0} {
	    $c coords mainLine0 $x0 $y0 $x0 $DG(ystart) $x1 $DG(ystart)
	    $c coords winBlock0 $x1 $DG(ystart)
	} else {
	    $c create line $x0 $Ys_old $x0 $DG(ystart) $x1 $DG(ystart) \
		    -width 4 -tags mainLine$i -arrow last \
		    -tags mainLine$i -fill #000
	}

	for {set j 0} {$j < $DG(n_subblock,$i)} {incr j} {	 
	    set ys [expr $DG(ystart) + $j * ($DG(bh) + $DG(yspace))]
	    DG_CreateSubBlock $c $xvl $ys $i $j
	}
	set Ys_old $DG(ystart)
	set DG(ystart) [expr $ys + 2 * $DG(bh) + $DG(yspace) + $DG(YSpace)]
	set Ys     $DG(ystart)
    }

    set bbox [$c bbox all]
    set x2 [expr $DG(bw) + $x0 + $xvl]
    set y2 [expr [lindex $bbox 3] + $DG(bh)]

    $c config -scrollregion [list 0 0 $x2 $y2]    
    
    DG_RadioBind $DG(r0) 0

    ## OK & Cancel Button
    set f3 [frame $t.f3]
    pack $f3 -side top -fill x
    DefaultButton $f3.ok -text OK -command [list DataGridOK $t]
    button $f3.can -text Cancel -command [list CancelProc $t]
    pack $f3.can $f3.ok -padx 10 -pady 10 -side left -expand 1
    # t.k.
    return $t
}

proc DG_CreateBlock {c xs ys ib} {
    global DG

    set f [frame $c.f$ib -relief raised -bd 1 -class DataGridBlock]
    $c create window $xs $ys -anchor w -window $f -tags winBlock$ib
    
    if ![info exists DG(blockFont)] {
	label $f.l$ib -text "Block #$ib;   dim = $DG(type,$ib)" -anchor w
	set DG(blockFont) [ModifyFont [$f.l$ib cget -font] $f.l$ib \
		-underline 1 -default 1]
	$f.l$ib config -font $DG(blockFont)
    } else {
	label $f.l$ib -text "Block #$ib;   dim = $DG(type,$ib)" \
		-font $DG(blockFont) -anchor w
    }
    
    set DG(r$ib) [radiobutton $f.r$ib -variable DG(radio) \
		      -text "Identifier:    $DG(ident,$ib)" \
		      -command DG_RadiobuttonCmd \
		      -value $ib -anchor w]
    bind $f.r$ib <Button-1> [list DG_RadioBind %W $ib]
    pack $f.l$ib $f.r$ib -side top -fill x

    update
    if ![info exists DG(bh)] {
	set DG(bh) [expr [winfo height $f] + \
		2 * ([$f cget -bd] + [$f cget -highlightthickness])]
	set DG(YSpace) 15; # y-space between block windows
	set DG(yspace) 5; # y-space between subblock windows
	set DG(ystart) [expr $ys + 19 + $DG(YSpace) + $DG(bh)/2]; #y center position for block0
    }
    set bw [expr [winfo width $f] + \
		2 * ([$f cget -bd] + [$f cget -highlightthickness])]
    if ![info exists DG(bw)] { 
	set DG(bw) $bw 
    } elseif { $bw > $DG(bw) } {
	set DG(bw) $bw
    }
}

proc DG_CreateSubBlock {c xs ys ib isb} {
    global DG xcColors
    
    set xsb [expr $xs + 20]
    set sf [frame $c.sf${ib}_$isb -relief raised -bd 1 -class DataGridSubBlock]
    $c create line $xs $ys \
	    $xs  [expr $ys + $DG(bh) + $DG(yspace)] \
	    $xsb [expr $ys + $DG(bh) + $DG(yspace)] \
	    -width 2 -tags horSubLine -arrow last -fill #000 
    $c create window $xsb [expr $ys + $DG(bh) + $DG(yspace)] -anchor w \
	    -window $sf -tags winSubBlock_${ib}_${isb}
    label $sf.l -text "Sub-block #$isb" -font $DG(blockFont) -anchor w
   
    if ![info exists DG(envar$ib,$isb)] {
	set DG(envar$ib,$isb) 1.0
    }
    frame $sf.1
    set DG(c$ib,$isb) [checkbutton $sf.r -variable DG(cb$ib,$isb) \
	    -text "Identifier:    $DG(subident,$ib,$isb)" -anchor w \
	    -command [list DG_CheckCom $ib $isb]]
    set DG(l$ib,$isb) [label $sf.l2 -text "  Multiply factor: " \
	    -fg $xcColors(disabled_fg)]
    set DG(e$ib,$isb) [entry $sf.e -relief sunken -bd 0 \
	    -textvariable DG(envar$ib,$isb) -width 4]

    pack $sf.l $sf.1 -side top -fill x
    pack $sf.r $sf.l2 -side left -fill x -in $sf.1   
    pack $sf.e -side left -fill x -in $sf.1 -padx 1m
    update
    set bw [expr [winfo width $sf] + \
	    2 * ([$sf cget -bd] + [$sf cget -highlightthickness])]
    if { $bw > $DG(bw) } {
	set DG(bw) $bw
    }
}


proc DG_RadioBind {w ib} {
    global DG xcColors

    for {set i 0} {$i < $DG(n_block)} {incr i} {	
	for {set j 0} {$j < $DG(n_subblock,$i)} {incr j} {	 
	    if { $i == $ib } {
		$DG(c$i,$j) config -state normal
		DG_CheckCom $i $j
	    } else {
		$DG(c$i,$j) config -state disabled
		$DG(l$i,$j) config -fg $xcColors(disabled_fg)
		$DG(e$i,$j) config -bd 0
	    }
	}
    }
}


proc DG_CheckCom {ib isb} {
    global DG xcColors

    for {set j 0} {$j < $DG(n_subblock,$ib)} {incr j} {	 
	if { $DG(cb$ib,$isb) } {
	    $DG(l$ib,$isb) config -fg $xcColors(enabled_fg)
	    $DG(e$ib,$isb) config -bd 1
	} else {
	    $DG(l$ib,$isb) config -fg $xcColors(disabled_fg)
	    $DG(e$ib,$isb) config -bd 0
	}
    }
}


proc DataGridOK {{t {}}} {
    global DG DataGrid
    
    set ib $DG(radio)
    set datarec $ib
    for {set i 0} {$i < $DG(n_subblock,$ib)} {incr i} {
	if $DG(cb$ib,$i) {
	    if { ! [number DG(envar$ib,$i) real] } { return }
	    append datarec " $i $DG(envar$ib,$i)"
	}
    }

    if { $t != "" } {
	CancelProc $t
    }

    eval xc_isodatagrid $datarec

    if { $DG(type,$ib) == "2D" } {
	
	set DataGrid(launch_command) IsoControl2D
	set DataGrid(dim) "2D"
	DataGrid2Isosurf
	IsoControl2D
	return "2D"
    } elseif { $DG(type,$ib) == "3D" } {
	set DataGrid(launch_command) IsoControl
	set DataGrid(dim) "3D"
	DataGrid2Isosurf
	IsoControl
	return "3D"
    }
}


proc DataGrid2Isosurf {} {
    global isosurf prop isosurf_struct isosign isodata
    
    set isosurf(minvalue)   [xc_iso minvalue]
    set isosurf(maxvalue)   [xc_iso maxvalue]    
    set isosurf(rangevalue) [expr $isosurf(maxvalue) - $isosurf(minvalue)]
    set prop(type_of_run)   RHF

    # set isosurf_struct(isosign) ""
    # set isosurf_struct(isodata) ""
    # set isosurf_struct(spin)    ""
    
    set isosurf_struct(3Dinterpl_degree) 1
    trace variable isosurf(3Dinterpl_degree) w xcTrace

    set isosign       ""
    set isodata       ""
    SetIsoSurfArray
    ConvertTwoSideVar; # to get correct value of isosurf(twoside_lighting)
    Set_UpdateIsosurf_Struct
}


proc DataGrid {} {
    global sInfo DG DataGrid
    
    if { ! [info exists DataGrid(first_time)] } {
	set DataGrid(first_time) exists
    } else {
	# DataGrid(first_time) already exists, simple launch the last datagrid
	eval $DataGrid(launch_command)
    }
	
    # .dg is the toplevel for DataGrid
    if { [winfo exists .dg] || [winfo exists .iso] || \
	    [winfo exists .iso2D] } { return }

    set info [xc_isodatagrid info]
    update

    if { [array exists DG] } { unset DG }
    # t.k.
    return [parseDGInfo $info]
}


proc DG_RadiobuttonCmd {} {
    global DG
    
    set DG(cb${DG(radio)},0) 1
}
