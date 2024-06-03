#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wigner.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc SetWignerSeitzInit {} {
    global periodic check ws

    if { $periodic(dim) < 3 } {
	return
    }
    if ![info exists periodic(igroup)] {
	return
    }

    if { $check(wigner) == 0 } {
	if ![info exists ws(not_config_yet)] {	
	    WignerSeitzInit
	} else {
	    if $ws(not_config_yet) {
		WignerSeitzInit
	    }
	}
    }

    set igroup $periodic(igroup)
    if { $igroup == 9 } {
	set igroup 8
    }	
    SetWignerSeitz $igroup
}


proc SetWignerSeitz nc {
    global ws wsc wsp ws_lfpos ws_npos

    set ws(c_type) $nc
    set ws(normal_pos_color)   "#ef0000"
    set ws(selected_pos_color) "#33ff66"
    #
    # cage definitions: BEGIN
    #
    set CON13 0.3333333333333333
    set CON23 0.6666666666666667

    set pc   1
    set ac   2
    set bc   3
    set cc   4
    set fc   5
    set ic   6
    set rc   7
    set hc   8
    set rcc  9
    set hcc  10
    set rpcc 11

    set ws_npos($pc)   8
    set ws_npos($ac)   10
    set ws_npos($bc)   10
    set ws_npos($cc)   10
    set ws_npos($ic)   9
    set ws_npos($fc)   14
    set ws_npos($hc)   12
    set ws_npos($rc)   10
    set ws_npos($hcc)  8
    set ws_npos($rcc)  8
    set ws_npos($rpcc) 8

    set lpos($pc) {
     -0.341936 -0.270062 -0.748429
     -0.614030 -0.589584  0.159244
     -0.341936  0.673770 -0.418063
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($ac) {
     -0.341936 -0.270062 -0.748429
     -0.481125  0.042093 -0.129410
     -0.614030 -0.589584  0.159244
     -0.341936  0.673770 -0.418063
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.481125 -0.042093  0.129410
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($bc) {
     -0.341936 -0.270062 -0.748429
      0.003142 -0.471916 -0.165183
     -0.614030 -0.589584  0.159244
     -0.341936  0.673770 -0.418063
     -0.003142  0.471916  0.165183
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($cc) {
     -0.341936 -0.270062 -0.748429
      0.136047  0.159761 -0.453837
     -0.614030 -0.589584  0.159244
     -0.136047 -0.159761  0.453837
     -0.341936  0.673770 -0.418063
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($fc) {
     -0.341936 -0.270062 -0.748429
      0.136047  0.159761 -0.453837
     -0.481125  0.042093 -0.129410
      0.003142 -0.471916 -0.165183
     -0.614030 -0.589584  0.159244
     -0.136047 -0.159761  0.453837
     -0.341936  0.673770 -0.418063
     -0.003142  0.471916  0.165183
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.481125 -0.042093  0.129410
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($ic) {
     -0.341936 -0.270062 -0.748429
      0.000000  0.000000  0.000000
     -0.614030 -0.589584  0.159244
     -0.341936  0.673770 -0.418063
     -0.614030  0.354248  0.489610
      0.614030 -0.354248 -0.489610
      0.341936 -0.673770  0.418063
      0.614030  0.589584 -0.159244
      0.341936  0.270062  0.748429
    }
    set lpos($hc) {
	-0.500000 -0.615661 -0.788011
	0.500000 -0.080352 -0.571732
	-0.500000  0.454958 -0.355452
	-0.500000 -0.990268  0.139173
	0.500000 -0.454958  0.355452
	-0.500000  0.080352  0.571732
	-1.500000  0.990268 -0.139173
	-1.500000  0.615661  0.788011
	1.500000 -0.615661 -0.788011
	1.500000 -0.990268  0.139173
	0.500000  0.990268 -0.139173
	0.500000  0.615661  0.788011
    }
    set lpos($hcc) {
	-0.500000 -0.615661 -0.788011
	-0.500000 -0.990268  0.139173
	-1.500000  0.990268 -0.139173
	-1.500000  0.615661  0.788011
	1.500000 -0.615661 -0.788011
	1.500000 -0.990268  0.139173
	0.500000  0.990268 -0.139173
	0.500000  0.615661  0.788011
    }
    set lpos($rc) {
	-0.500000 -0.615661 -0.788011
	0.500000 -0.205220 -0.262670
	-0.500000  0.205220  0.262670
	-0.500000 -0.990268  0.139173
	-1.500000  0.990268 -0.139173
	-1.500000  0.615661  0.788011
	1.500000 -0.615661 -0.788011
	1.500000 -0.990268  0.139173
	0.500000  0.990268 -0.139173
	0.500000  0.615661  0.788011
    }
    set lpos($rcc) {
	-0.500000 -0.615661 -0.788011
	-0.500000 -0.990268  0.139173
	-1.500000  0.990268 -0.139173
	-1.500000  0.615661  0.788011
	1.500000 -0.615661 -0.788011
	1.500000 -0.990268  0.139173
	0.500000  0.990268 -0.139173
	0.500000  0.615661  0.788011
    }
    set lpos($rpcc) {
	0.212675 -0.368363 -2.412276
	0.063295 -1.264331 -0.630444
	-0.921512  0.441406 -0.977740
	-1.070892 -0.454562  0.804092
	1.070892  0.454562 -0.804092
	0.921512 -0.441406  0.977740
	-0.063295  1.264331  0.630444
	-0.212675  0.368363  2.412276
    }
    
    set ws_lfpos($pc) {
      0.000000  0.000000  0.000000
      0.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($ac) {
      0.000000  0.000000  0.000000
      0.000000  0.500000  0.500000
      0.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.500000  0.500000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($bc) {
      0.000000  0.000000  0.000000
      0.500000  0.000000  0.500000
      0.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      0.500000  1.000000  0.500000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($cc) {
      0.000000  0.000000  0.000000
      0.500000  0.500000  0.000000
      0.000000  0.000000  1.000000
      0.500000  0.500000  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($fc) {
      0.000000  0.000000  0.000000
      0.500000  0.500000  0.000000
      0.000000  0.500000  0.500000
      0.500000  0.000000  0.500000
      0.000000  0.000000  1.000000
      0.500000  0.500000  1.000000
      0.000000  1.000000  0.000000
      0.500000  1.000000  0.500000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.500000  0.500000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($ic) {
      0.000000  0.000000  0.000000
      0.500000  0.500000  0.500000
      0.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($hc) {
      0.000000  0.000000  0.000000
      0.666667  0.333333  0.000000
      0.333333  0.666667  0.000000
      0.000000  0.000000  1.000000
      0.666667  0.333333  1.000000
      0.333333  0.666667  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    set ws_lfpos($rc) {
      0.000000  0.000000  0.000000
      0.666667  0.333333  0.333333
      0.333333  0.666667  0.666667
      0.000000  0.000000  1.000000
      0.000000  1.000000  0.000000
      0.000000  1.000000  1.000000
      1.000000  0.000000  0.000000
      1.000000  0.000000  1.000000
      1.000000  1.000000  0.000000
      1.000000  1.000000  1.000000
    }
    #
    # cage definitions: END
    #
    
    #
    # Widgets building block
    #
    if [winfo exists .wgnset] { return }
    set t [xcToplevel .wgnset \
	    "Wigner-Seitz Cell" "Wigner-Seitz Sell" . 0 0 1]

    set f0 [frame $t.f0]
    set f1 [frame $t.f1]
    set f2 [frame $t.f2 -relief raised -bd 2]

    pack $f0 $f1 $f2 -side top -fill both

    proc SetWignerSeitz_Show {type f f1 f2 b1 b2} {
	global ws wsp wsc
	
	set ws(type) $type

	pack forget $f2
	pack $f1 -in $f -side top -fill both
	$b1 configure -bd 3
	$b2 configure -bd 1
    }

    set ws(prim_b) [button $f0.prim \
	    -text "Wigner-Seitz settings\nfor primitive cell mode" \
	    -command [list SetWignerSeitz_Show \
	    prim $f1 $t.f1p $t.f1c $f0.prim $f0.conv] \
	    -highlightthickness 0 \
	    -bd 3]
    set ws(conv_b) [button $f0.conv \
	    -text "Wigner-Seitz settings\nfor conventional cell mode" \
	    -command [list SetWignerSeitz_Show \
	    conv $f1 $t.f1c $t.f1p $f0.conv $f0.prim] \
	    -highlightthickness 0 \
	    -bd 1]
    pack $ws(prim_b) $ws(conv_b) -side left -fill both

    set f1p [frame $t.f1p -relief raised -bd 2]
    set f1c [frame $t.f1c -relief raised -bd 2]
    pack $f1p -in $f1 -side top -fill both

    set can [button $f2.can -text "Cancel" \
	    -command [list SetWignerSeitz_Cancel $t]]
    set tst [button $f2.tst -text "Test It" \
	    -command [list SetWignerSeitz_OK test $t]]
    set ok  [button $f2.ok  -text "OK" \
	    -command [list SetWignerSeitz_OK OK $t]]
    pack $can $tst $ok -side left -expand 1 -padx 5 -pady 5

    #
    # make widgets for primitive & conventional modes
    #
    if { ![info exists wsp(mode)] } {
	set wsp(mode) "every"
    }
    if { ![info exists wsc(mode)] } {
	set wsc(mode) "every"
    }
    foreach f [list [list $f1c wsc] [list $f1p wsp]] { 
	set w       [lindex $f 0]
	set varname [lindex $f 1]
	if { $varname == "wsp" } {
	    if { $nc == $hc } {
		set nc $hcc
	    } elseif { $nc == $rc } {
		set nc $rpcc
	    } else {
		set nc $pc
	    }
	}
	
	set pp $pc
	if { $nc == $hc || $nc == $hcc } {
	    set pp $hcc
	} elseif { $nc == $rc } {
	    set pp $rcc
	} elseif { $nc == $rpcc } {
	    set pp $rpcc
	}

	set off 2.0
	if { $nc == $hc || $nc == $rc || $nc == $hcc || $nc == $rpcc } {
	    set off 4.0
	}
	set off2 [expr $off / 2.0]

	puts stdout "$varname $nc"
	upvar #0 $varname var
	set f1l [frame $w.l -relief sunken -bd 2]
	set f1r [frame $w.r]
	pack $f1l $f1r -side left -fill both -padx 2m -pady 2m
    
	#
	# Canvas
	#
	set w 170
	set h 170
	set var(can) [canvas $f1l.c -width $w -height $h -bg "#ffffff"]
	pack $var(can)
	
	set rfw $w
	set rfh $h
		
	set Zmin 0.0
	set Zmax 0.0
	for {set i 0} {$i < $ws_npos($nc)} {incr i} {
	    set pos($i,0) [lindex $lpos($nc) [expr $i * 3 + 0]]
	    set pos($i,1) [lindex $lpos($nc) [expr $i * 3 + 1]]
	    set pos($i,2) [lindex $lpos($nc) [expr $i * 3 + 2]]
	    if { $pos($i,2) < $Zmin } { set Zmin $pos($i,2) }
	    if { $pos($i,2) > $Zmax } { set Zmax $pos($i,2) }
	}
	set dZ [expr $Zmax - $Zmin]
	
	for {set i 0} {$i < $ws_npos($pp)} {incr i} {
	    set p($i,0) [lindex $lpos($pp) [expr $i * 3 + 0]]
	    set p($i,1) [lindex $lpos($pp) [expr $i * 3 + 1]]
	    set p($i,2) [lindex $lpos($pp) [expr $i * 3 + 2]]
	    
	    set p($i,0) [expr ($p($i,0) + $off2) * $rfw / $off]
	    set p($i,1) [expr ($off2 - $p($i,1)) * $rfh / $off]
	    
	    #puts stdout "$i $p($i,0) $p($i,1)"
	}
	
	$var(can) create polygon \
		$p(1,0) $p(1,1) $p(3,0) $p(3,1) \
		$p(7,0) $p(7,1) $p(5,0) $p(5,1) \
		-fill {} -outline $ws(normal_pos_color) -width 3 -tags cagepoly
	$var(can) create polygon \
		$p(3,0) $p(3,1) $p(2,0) $p(2,1) \
		$p(6,0) $p(6,1) $p(7,0) $p(7,1) \
		-fill {} -outline $ws(normal_pos_color) -width 3 -tags cagepoly
	$var(can) create polygon \
		$p(4,0) $p(4,1) $p(5,0) $p(5,1) \
		$p(7,0) $p(7,1) $p(6,0) $p(6,1) \
		-fill {} -outline $ws(normal_pos_color) -width 3 -tags cagepoly
	
	$var(can) create line $p(0,0) $p(0,1) \
		$p(1,0) $p(1,1) -fill $ws(normal_pos_color) \
		-width 1 -arrow last -arrowshape {15 15 4} -tags cageline
	$var(can) create line $p(0,0) $p(0,1) \
		$p(2,0) $p(2,1) -fill $ws(normal_pos_color) \
		-width 1 -arrow last -arrowshape {15 15 4} -tags cageline
	$var(can) create line $p(0,0) $p(0,1) \
		$p(4,0) $p(4,1) -fill $ws(normal_pos_color) \
		-width 1 -arrow last -arrowshape {15 15 4} -tags cageline
	
	$var(can) create text [expr $p(4,0) + 8] $p(4,1)  \
		-text a -anchor w -tags text
	$var(can) create text $p(2,0) [expr $p(2,1) - 8]  \
		-text b -anchor s -tags text
	$var(can) create text [expr $p(1,0) - 11] $p(1,1) \
		-text c -anchor e -tags text
		
	for {set i 0} {$i < $ws_npos($nc)} {incr i} {
	    set pos($i,0) [expr ($pos($i,0) + $off2) * $rfw / $off]
	    set pos($i,1) [expr ($off2 - $pos($i,1)) * $rfh / $off]
	    set r [expr 4.0 + ($pos($i,2) - $Zmin) / $dZ * 3.0]
	    $var(can) create oval \
		    [expr $pos($i,0) - $r] [expr $pos($i,1) - $r] \
		    [expr $pos($i,0) + $r] [expr $pos($i,1) + $r] \
		    -fill $ws(normal_pos_color) \
		    -outline $ws(normal_pos_color) -tags "cage pos p$i"
	    #puts stdout "$i $pos($i,0) $pos($i,1)"
	}
	
	#
	# maybe Wigner-Seitz cell option was not configured yet;
	# and this is consistent with the default settings 
	# the DEFAULT setting is: the 0 0 0 position selected
	#
	if { $ws(not_config_yet) } {
	    set id [$var(can) find withtag p0]
	    set var(sel,$id) 1
	    $var(can) itemconfigure $id -fill $ws(selected_pos_color) \
		    -outline $ws(selected_pos_color)
	} else {
	    foreach id [$var(can) find withtag pos] {
		if [info exists var(sel,$id)] {
		    if $var(sel,$id) {
			$var(can) itemconfigure $id \
				-fill $ws(selected_pos_color) \
				-outline $ws(selected_pos_color)
		    }
		}
	    }
	}

	$var(can) bind pos <1> [list SetWignerSeitz_BindCan %W %x %y] 
	#
	# left frame:: radiobuttons, checkbuton & button
	#
	set f1rt [frame $f1r.t -relief groove -bd 2]
	pack $f1rt -side top -padx 2m -pady 0

	radiobutton $f1r.rb0 \
		-text "Display Wigner-Seitz cell on every node" \
		-value "every" \
		-variable ${varname}(mode) \
		-command [list SetWignerSeitz_Radio $var(can)] \
		-anchor w
	radiobutton $f1r.rb1 \
		-text "Display Wigner-Seitz on selected node" \
		-value "selected" \
		-variable ${varname}(mode) \
		-command [list SetWignerSeitz_Radio $var(can)] \
		-anchor w
	
	set f1rf [frame $f1r.f]
	checkbutton $f1rf.ck \
		-text "Transparent Wigner-Seitz cell" -anchor w \
		-variable ws(transparent)
	proc SetWignerSeitz_ColorOK {type t} {
	    global ws mody_col
	    
	    if { $type == "OK" } {
		set cID [xcModifyColorGetID]
		set ws(color) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	    }
	    destroy $t
	}
	proc SetWignerSeitz_Color t {
	    global ws

	    set t [xcToplevel [WidgetName] \
		    "Set Wigner-Seitz Cell's Color" \
		    "Set Wigner-Seitz Cell's Color" $t 300 90 1]
	    xcModifyColor $t "Set Wigner-Seitz Cell's Color:" \
		    $ws(color) \
		    groove left left 100 100 70 5 20

	    set ok  [DefaultButton [WidgetName $t] -text "OK" \
		    -command [list SetWignerSeitz_ColorOK OK $t]]
	    set can [button [WidgetName $t] -text "Cancel" \
		    -command [list SetWignerSeitz_ColorOK Cancel $t]]
	    pack $ok $can -padx 10 -pady 10 -expand 1
	}

	button $f1rf.b -text "Color" -command [list SetWignerSeitz_Color $t]
	pack $f1r.rb0 $f1r.rb1 -in $f1rt -side top -fill x -pady 10
	pack $f1rf -side top -fill x -pady 10
	pack $f1rf.ck $f1rf.b -side left -expand 1
    }

    set ws(type) prim
    SetWignerSeitz_Radio $wsp(can)
    set ws(type) conv
    SetWignerSeitz_Radio $wsc(can)

    SetWignerSeitz_Show \
	    prim $f1 $t.f1p $t.f1c $f0.prim $f0.conv
}

proc SetWignerSeitz_Radio can {
    global ws wsp wsc
    
    if { $ws(type) == "prim" } {
	set mode $wsp(mode)
	set varname wsp
    } else {
	set mode $wsc(mode)
	set varname wsc
    }
    upvar #0 $varname var
    if { $mode == "every" } {
	# cageline
	foreach item [$can find withtag cageline] {
	    set ${varname}(color,$item) [$can itemcget $item -fill]
	}
	foreach item [$can find withtag cagepoly] {
	    set ${varname}(color,$item) [$can itemcget $item -outline]
	}
	foreach item [$can find withtag text] {
	    set ${varname}(color,$item) [$can itemcget $item -fill]
	}
	foreach item [$can find withtag pos] {
	    set ${varname}(color,$item) [$can itemcget $item -fill]
	}
	
	$can itemconfigure cageline -fill $ws(can_disable_color)
	$can itemconfigure cagepoly -outline $ws(can_disable_color)
	$can itemconfigure text -fill $ws(can_disable_color)
	$can itemconfigure pos -fill $ws(can_disable_color) \
		-outline $ws(can_disable_color)	
    } else {
	foreach item [$can find withtag cageline] {
	    if [info exists var(color,$item)] {
		$can itemconfigure $item -fill $var(color,$item)
	    }
	}
	foreach item [$can find withtag cagepoly] {
	    if [info exists var(color,$item)] {	    
		$can itemconfigure $item -outline $var(color,$item)
	    }
	}
	foreach item [$can find withtag text] {
	    if [info exists var(color,$item)] {
		$can itemconfigure $item -fill $var(color,$item)
	    }
	}
	foreach item [$can find withtag pos] {
	    if [info exists var(color,$item)] {
		$can itemconfigure $item -fill $var(color,$item) \
		    -outline $var(color,$item)
	    }
	}
    }
}

proc SetWignerSeitz_BindCan {can x y} {
    global ws wsp wsc

    if { $ws(type) == "prim" } {
	set mode $wsp(mode)
	set varname wsp
    } else {
	set mode $wsc(mode)
	set varname wsc
    }
    upvar #0 $varname var

    if { $mode == "every" } {
	return
    }
    
    set id [$can find closest $x $y]
    
    # did "find closest find the "pos" item or some other item ???
    set find 0
    foreach item [$can find withtag pos] {
	if { $id == $item } {
	    set find 1
	}
    }
    if !$find { return }

    if ![info exists var(sel,$id)] {
	set var(sel,$id) 0
    }

    #
    # if position was unselected -> select it
    # if position was selected   -> unselect it
    #
    if $var(sel,$id) {
	# UNSELECT
	$can itemconfigure $id -fill $ws(normal_pos_color) \
		-outline $ws(normal_pos_color)
	set var(sel,$id) 0
    } else {
	# SELECT
	$can itemconfigure $id -fill $ws(selected_pos_color) \
		-outline $ws(selected_pos_color)
	set var(sel,$id) 1
    }
}


proc SetWignerSeitz_OK {mode {t {}}} {
    global ws wsp wsc ws_lfpos ws_npos system radio check

    #
    # now Wigner_Seitz settings has been configured
    #
    set ws(not_config_yet) 0
    proc GetWignerSeitz_NSel varname {
	upvar #0 $varname var
	set n 0
	if ![winfo exists $var(can)] { return 0 }
	foreach item [$var(can) find withtag pos] {
	    if [info exists var(sel,$item)] {
		if $var(sel,$item) {
		    incr n
		}
	    }
	}
	return $n
    }	
    #
    # settings for PRIMITIVE cell mode
    #
    set render after
    if { $radio(cellmode) == "prim" } {
	set render now
    }
    if { [info exists wsp(mode)] && $wsp(mode) == "selected" } {
	if { [winfo exists $wsp(can)] } {
	    set n     [GetWignerSeitz_NSel wsp]
	    set outID [open $system(SCRDIR)/xc_ndsfp.$system(PID) w]
	    puts $outID "$n"
	    for {set i 0} {$i < $ws_npos($ws(pc))} {incr i} {
		set id [$wsp(can) find withtag p$i]
		if $wsp(sel,$id) {
		    puts $outID "\
			[lindex $ws_lfpos($ws(pc)) [expr $i * 3 + 0]] \
			[lindex $ws_lfpos($ws(pc)) [expr $i * 3 + 1]] \
			[lindex $ws_lfpos($ws(pc)) [expr $i * 3 + 2]]"
		}
	    }
	    flush $outID
	    close $outID
	}
	xc_wigner .mesa prim \
	    -nodesfile    $system(SCRDIR)/xc_ndsfp.$system(PID) \
	    -transparency $ws(transparent) \
	    -color        [rgb_h2f $ws(color)] \
	    -render       $render
    } else {
	xc_wigner .mesa prim \
	    -transparency $ws(transparent) \
	    -color        [rgb_h2f $ws(color)] \
	    -render       $render
    }

    #
    # settings for CONVENTIONAL cell mode
    #    
    set render after
    if { $radio(cellmode) == "conv" } {
	set render now
    }
    if { [info exists wsc(mode)] && $wsc(mode) == "selected" } {
	if { [winfo exists $wsp(can)] } {
	    set n     [GetWignerSeitz_NSel wsc]
	    set outID [open $system(SCRDIR)/xc_ndsfc.$system(PID) w]
	    puts $outID "$n"
	    for {set i 0} {$i < $ws_npos($ws(c_type))} {incr i} {
		set id [$wsc(can) find withtag p$i]
		if $wsc(sel,$id) {
		    puts $outID "\
			[lindex $ws_lfpos($ws(c_type)) [expr $i * 3 + 0]] \
			[lindex $ws_lfpos($ws(c_type)) [expr $i * 3 + 1]] \
			[lindex $ws_lfpos($ws(c_type)) [expr $i * 3 + 2]]"
		}
	    }
	    flush $outID
	    close $outID
	}
	xc_wigner .mesa conv \
	    -nodesfile    $system(SCRDIR)/xc_ndsfc.$system(PID) \
	    -transparency $ws(transparent) \
	    -color        [rgb_h2f $ws(color)] \
	    -render       $render
    } else {
	xc_wigner .mesa conv \
	    -transparency $ws(transparent) \
	    -color        [rgb_h2f $ws(color)] \
	    -render       $render
    }

    if { $mode == "OK" } {
	if { [winfo exists $t] } { destroy $t }
	set check(wigner) 1	
    }
}


proc SetWignerSeitz_Cancel {{t {}}} {
    global check wsp ws

    xc_wigner .mesa clear
    set check(wigner) 0
    # this is temporary solution
    # set ws(not_config_yet) 1
    if { $t != {} } {
	CancelProc $t
    }
}


proc WignerSeitz {} {
    global ws wsp wsc periodic check system radio
    
    #
    # so far just for CRYSTALS
    #
    if { $periodic(dim) < 3 } {
	set check(wigner) 0
	return
    }

    #
    # maybe we are in "Wigner-Seitz Settings" mode
    #
    if { [winfo exists .wgnset] } { 
	if { $check(wigner) == 0 } {
	    set check(wigner) 1
	} else {
	    set check(wigner) 0
	}
	return
    }
    if { ![info exists ws(not_config_yet)] } {	
	set ws(not_config_yet) 1
    }

    if { $ws(not_config_yet) } {
	#WignerSeitzInit
    
	if $check(wigner) {
	    #
	    # INITIAL settings for PRIMITIVE cell mode
	    #    
	    set wsp(mode) "selected"
	    set outID [open $system(SCRDIR)/xc_ndsfp.$system(PID) w]
	    puts $outID "1\n0.0 0.0 0.0"
	    flush $outID
	    close $outID
	    
	    set render after
	    if { $radio(cellmode) == "prim" } {
		set render now
	    }
	    xc_wigner .mesa prim \
		    -nodesfile    $system(SCRDIR)/xc_ndsfp.$system(PID) \
		    -transparency $ws(transparent) \
		    -color        [rgb_h2f $ws(color)] \
		    -render       $render
	    
	    #
	    # INITIAL settings for CONVENTIONAL cell mode
	    #
	    set wsc(mode) "selected"
	    set outID [open $system(SCRDIR)/xc_ndsfc.$system(PID) w]
	    puts $outID "1\n0.0 0.0 0.0"
	    flush $outID
	    close $outID
	    
	    set render after
	    if { $radio(cellmode) == "conv" } {
		set render now
	    }
	    xc_wigner .mesa conv \
		    -nodesfile    $system(SCRDIR)/xc_ndsfc.$system(PID) \
		    -transparency $ws(transparent) \
		    -color        [rgb_h2f $ws(color)] \
		    -render $render
	}
    } elseif $check(wigner) {
	SetWignerSeitz_OK test	
    }
    
    if { !$check(wigner) } {
	SetWignerSeitz_Cancel
    }
}


proc WignerSeitzInit {} {
    global ws wsp wsc

    xcDebug "In WignerSeitzInit"

    set ws(pc)   1
    set ws(ac)   2
    set ws(bc)   3
    set ws(cc)   4
    set ws(fc)   5
    set ws(ic)   6
    set ws(rc)   7
    set ws(hc)   8
    set ws(rcc)  9
    set ws(hcc)  10
    set ws(rpcc) 11

    set ws(normal_pos_color)   "#ef0000"
    set ws(selected_pos_color) "#33ff66"    
    set ws(can_disable_color) #eeeeee
    set ws(color) #55eeff
    set ws(transparent) 0
    
    #
    # initialize wsx(sel,$i); 50 is more then enough 
    # (there is never more then 50 atoms)
    for {set i 0} {$i < 50} {incr i} {
	set wsp(sel,$i) 0
	set wsc(sel,$i) 0
    }
    set ws(not_config_yet) 1
}
