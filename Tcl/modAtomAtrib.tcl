#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/modAtomAtrib.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#=============================================================================#
#             ATOMIC COLORS --- ATOMIC COLORS --- ATOMIC COLORS               #
#=============================================================================#
proc ModAtomCol {} {
    global atcol Alist mody
    # Alist - AtomNameLists

    set top .atomcol
    # just in case in window already exists
    if { [winfo exists $top] } { return } 
    toplevel .atomcol
    wm title $top "Atomic Colors"
    wm iconname $top "Atomic Colors"
    #grab $top

    # place $top according to "."
    xcPlace . $top 100 50

    # there will be three frames left, right and bottom
    set l [frame $top.l -relief raised -bd 2]
    set r [frame $top.r -relief raised -bd 2]
    set b [frame $top.b -relief raised -bd 2]

    pack $b -side bottom -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1
    pack $l $r -side left -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1

    # in LEFT-FRAME there will be Entrie & scrolllistbox
    Entries $l "Atom:" atcol(atom) 6
    bind $l.frame.entry1 <Return> { AtomColBind }

    set lbox [ScrolledListbox2 $l.lb -width 12 -height 20 -setgrid true]
    bind $lbox <ButtonRelease-1> [list AtomColSelect %W %y]

    # we must rebuilt Alist so that after each atom we add "\n"
    set n 0
    foreach elem $Alist {
	if { $n < 10 } {
	    lappend Atoms "  $n:         $elem"
	} elseif { $n < 100 } {
	    lappend Atoms  " $n:         $elem"
	} else {
	    lappend Atoms   "$n:         $elem"
	}	
	incr n
    }
    eval {$lbox insert 0} $Atoms
    
    # RIGHT-FRAME --- RIGHT-FRAME
    # before we select any atom from above Listbox, a default will be H
    set color        [xc_getvalue $mody(D_ATCOL_ONE) 0]
    set atcol(red)   [lindex $color 0]
    set atcol(green) [lindex $color 1]
    set atcol(blue)  [lindex $color 2]     
    set atcol(atom)  "X"
    set atcol(nat)   0

    set atcol(hxred)   [ d2h [expr int($atcol(red) * 255)  ] ] 
    set atcol(hxgreen) [ d2h [expr int($atcol(green) * 255)] ] 
    set atcol(hxblue)  [ d2h [expr int($atcol(blue) * 255) ] ] 
    
    set fr  [frame $r.1 -relief sunken -bd 2]
    set col [frame $fr.col  -bd 0 \
	    -bg "#$atcol(hxred)$atcol(hxgreen)$atcol(hxblue)" \
	    -width 150 -height 150]
    
    scale $r.red -from 0 -to 1 -length 150 -variable atcol(red) \
	    -orient horizontal -label "Red:" -tickinterval 1.0 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command Color
    scale $r.green -from 0 -to 1 -length 150 \
	    -variable atcol(green) \
	    -orient horizontal -label "Green:" -tickinterval 0.5 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command Color
    scale $r.blue -from 0 -to 1 -length 150 \
	    -variable atcol(blue) \
	    -orient horizontal -label "Blue:" -tickinterval 0.5 \
	    -digits 3 -resolution 0.01 -showvalue true -width 10 \
	    -command Color

    pack $fr -side top -fill both -expand 1 -padx 30 -pady 15 -ipadx 0 -ipady 0
    pack $col -side top -fill both -expand 1 -padx 0 -pady 0 
    pack $r.red $r.green $r.blue -side top -fill both -expand 1

    # BOTTOM FRAME --- BOTTOM FRAME
    # OK, Default & Cancel Button
    set col [button $b.upd -text "Update Color" \
	    -command AtomColUpD]
    set def [button $b.def -text "Default Color" \
	    -command [list AtomColLoad default]]
    set can [button $b.can -text "Reset All" \
	    -command AtomColReset]
    set clo [button $b.clo -text "Close" -command [list AtomColOK $top]]
    
    pack $col $def $can $clo -side left -expand 1 -pady 10 -padx 5
            
    return
}


proc AtomColBind {} {
    global atcol

    set atcol(nat) [Aname2Nat $atcol(atom)]
    if ![regexp {[0-9]+} $atcol(nat)] {
	tk_dialog .atradup "ERROR" $atcol(nat) error 0 OK
	return
    }
    
    AtomColLoad
    return
}


proc Color { {trash {}} } {
    global atcol

    set atcol(hxred) [ d2h [expr int($atcol(red) * 255)] ]
    set atcol(hxgreen) [ d2h [expr int($atcol(green) * 255)] ]
    set atcol(hxblue) [ d2h [expr int($atcol(blue) * 255)] ] 
    
    .atomcol.r.1.col configure \
	    -bg "#$atcol(hxred)$atcol(hxgreen)$atcol(hxblue)"

} 


proc AtomColSelect { w y } {
    global atcol mody

    # $w is the name of listbox widget who contains groups
    # $y is the vertical position of "selection"
    puts stdout "$w select anchor [$w nearest $y]"
    $w select anchor [$w nearest $y]
    set nline [$w curselection]
    set atom [$w get $nline]

    # now we must purify atom variable, because it's sometning like "  8:  O"
    # we will take a number instead of name (communication with xcrys)
    regexp {[0-9]+} $atom atcol(nat)
    
    # to update entry we also need atom name
    regexp {[A-Za-z]+} $atom atcol(atom)

    AtomColLoad

    return
}


proc AtomColLoad { {type {}} } {
    global mody atcol
    
    if { $type == "default"} {
	set color [xc_getdefault $mody(D_ATCOL_ONE) $atcol(nat)]	
    } else {
	puts stdout "GETVALUE::"
	set color [xc_getvalue $mody(D_ATCOL_ONE) $atcol(nat)]
    }

    set atcol(red) [lindex $color 0]
    set atcol(green) [lindex $color 1]
    set atcol(blue) [lindex $color 2] 

    Color
    puts stdout "ATOM=\"$atcol(nat)\"; NAME=$atcol(atom)"
    flush stdout
}


proc AtomColOK {top} {
    global atcol

    if { [winfo exists $top] } { 
	#grab release $top
	destroy $top
    }
    return
}


proc AtomColUpD {} {
    global mody atcol

    puts stdout "xc_newvalue .mesa $mody(L_ATCOL_ONE) $atcol(nat) $atcol(red) $atcol(green) $atcol(blue)"
    flush stdout
    xc_newvalue .mesa $mody(L_ATCOL_ONE) $atcol(nat) \
	    $atcol(red) $atcol(green) $atcol(blue)

    return
}


proc AtomColReset {} {
    global mody atcol
    
    # reset atomic colors
    xc_resetvar .mesa $mody(R_ATCOL)
    # update display
    AtomColLoad

    return
}


#=============================================================================#
#               ATOMIC RADII --- ATOMIC RADII --- ATOMIC RADII                #
#=============================================================================#
proc ModAtomRad {} {
    global atrad Alist mody
    # Alist - AtomNameLists

    set top .atomrad
    # just in case in window already exists
    if { [winfo exists $top] } { return } 
    toplevel $top
    wm title $top "Atomic Radii"
    wm iconname $top "Atomic Radii"
    #grab $top

    # place $top according to "."
    xcPlace . $top 100 50

    # there will be three frames left, right and bottom
    set l [frame $top.l -relief raised -bd 2]
    set r [frame $top.r -relief raised -bd 2]
    set b [frame $top.b -relief raised -bd 2]

    pack $b -side bottom -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1
    pack $l $r -side left -padx 0 -pady 0 -ipadx 0 -ipady 0 \
	    -fill both -expand 1
    
    # in LEFT-FRAME there will be a scrolllistbox
    set lbox [ScrolledListbox2 $l.lb -width 12 -height 15 -setgrid true]
    bind $lbox <ButtonRelease-1> [list AtomRadSelect %W %y]

    # we must rebuilt Alist so that after each atom we add "\n"
    set n 0
    foreach elem $Alist {
	if { $n < 10 } {
	    lappend Atoms "  $n:         $elem"
	} elseif { $n < 100 } {
	    lappend Atoms  " $n:         $elem"
	} else {
	    lappend Atoms   "$n:         $elem"
	}	
	incr n
    }
    eval {$lbox insert 0} $Atoms
    
    # RIGHT-FRAME --- RIGHT-FRAME
    # before we select any atom from above Listbox, a default will be X
    set atrad(covf)   [xc_getvalue $mody(D_COVF)]
    set atrad(scale)  [xc_getvalue $mody(D_ATRAD_SCALE)]
    set atrad(atom)   "X"
    set atrad(rad)    [xc_getvalue $mody(D_ATRAD_ONE) 0]
    set atrad(covrad) [xc_getvalue $mody(D_RCOV_ONE)  0]
    set atrad(nat)    0

    set atrad(covf_old)  $atrad(covf)
    set atrad(scale_old) $atrad(scale)

    set f11 [frame $r.1 -relief groove -bd 2]
    Entries $f11 {"Chemical connectivity factor:"} atrad(covf) 8 0
    set b12 [button $f11.frame.b12 -text "Clear" \
	    -command [list AtRadClr atrad(covf) $f11.frame.entry1]]
    set b11 [button $f11.frame.b11 -text "Default" \
	    -command [list AtRadDef atrad(covf) D_COVF]]
    set b13 [button $f11.frame.b13 -text "Update" \
		 -command AtRadUpdateConF]
    pack $f11 -side top -padx 5 -pady 5 -fill x -ipady 10
    pack $b12 $b11 $b13 -side left -side left -padx 5 -pady 5 
    focus $f11.frame.entry1

    set f21 [frame $r.2 -relief groove -bd 2]
    Entries $f21 {"SpaceFill/Ball scale factor:"} atrad(scale) 8 0
    set b22 [button $f21.frame.b22 -text "Clear" \
	    -command [list AtRadClr atrad(scale) $f21.frame.entry1]]
    set b21 [button $f21.frame.b21 -text "Default" \
	    -command [list AtRadDef atrad(scale) D_ATRAD_SCALE]]
    set b23 [button $f21.frame.b23 -text "Update" \
	    -command AtRadUpdateScale]
    pack $f21 -side top -padx 5 -pady 5 -expand 1 -fill x -ipady 10
    pack $b22 $b21 $b23 -side left -side left -padx 5 -pady 5
    
    # t.k.: start here ...
    set f31 [frame $r.3 -relief groove -bd 2]
    set lbl [label $f31.lbl -text "Change radius:" -relief flat]
    pack $lbl -side top
    #set f311 [frame $f31.1]
    #set f312 [frame $f31.2]
    #pack $f311 $f312 -side top -fill x
    Entries $f31 {"Atom:" "Display radius:" "Covalent radius:"} {
	{atrad(atom)} {atrad(rad)} {atrad(covrad)}
    } 8 0

    set b3 [button $f31.frame.b3 -text "Clear" -command [list AtRadClr atrad(rad) $f31.frame.entry2]]
    pack $b3 -side left -side left -padx 5 -pady 5

    set b31 [button $f31.b31 -text "Update Radius"  -command AtRadUpD]
    set res [button $f31.res -text "Reset All"      -command AtRadReset]
    set b32 [button $f31.b32 -text "Default Radius" -command [list AtRadLoad default]]
    pack $f31 -side top -padx 5 -pady 5 -fill x -ipady 10
    pack $b31 $res $b32  -side left -side left -padx 5 -pady 5 -expand 1
    # /t.k.
    
    # BOTTOM FRAME --- BOTTOM FRAME
    # OK button
    set clo [button $b.clo -text "OK" -command [list AtRadOK $top]]
    
    pack $clo -side left -expand 1 -pady 10 -padx 5
            
    return
}


proc AtRadUpdateConF {} {
    global atrad mody
    if { $atrad(covf_old) != $atrad(covf) } {
	xc_newvalue .mesa $mody(L_COV_SCALE) $atrad(covf)
    }
    set atrad(covf_old)  $atrad(covf)
}


proc AtRadUpdateScale {} {    
    global atrad mody
    if { $atrad(scale_old) != $atrad(scale) } {
	xc_newvalue .mesa $mody(L_ATRAD_SCALE) $atrad(scale)
    }   
    set atrad(scale_old) $atrad(scale)
}


proc AtRadUpD {} {
    global mody atrad

    # from atrad(atom) --> atrad(nat)
    set atrad(nat) [Aname2Nat $atrad(atom)]
    if ![regexp {[0-9]+} $atrad(nat)] {
	tk_dialog .atradup "ERROR" $atrad(nat) error 0 OK
	return
    }
    
    xc_newvalue .mesa $mody(L_ATRAD_ONE) $atrad(nat) $atrad(rad)
    xc_newvalue .mesa $mody(L_RCOV_ONE)  $atrad(nat) $atrad(covrad)
    
    return
}


proc AtRadLoad { {type {}} } {
    global mody atrad
    
    # from atrad(atom) -> atrad(nat)
    puts stdout "AtRadLoad:: $atrad(nat)"
    flush stdout
    if { $type == "default"} {
	set atrad(rad)    [xc_getdefault $mody(D_ATRAD_ONE) $atrad(nat)]
	set atrad(covrad) [xc_getdefault $mody(D_RCOV_ONE) $atrad(nat)]
    } else {
	set atrad(rad)    [xc_getvalue $mody(D_ATRAD_ONE) $atrad(nat)]
	set atrad(covrad) [xc_getvalue $mody(D_RCOV_ONE) $atrad(nat)]
    }
    return
}


proc AtRadReset {} {
    global mody

    # reset atomic radius
    xc_resetvar .mesa $mody(R_ATRAD)
    # t.k.: insert the RCOV code here ...
    xc_resetvar .mesa $mody(R_RCOV)
    # update a display
    AtRadLoad

    return
}    


proc AtRadOK {top} {
    global atrad mody
    
    puts stdout "AtRadOK:: mody(L_COV_SCALE) = $atrad(covf)"
    puts stdout "AtRadOK:: mody(L_ATRAD_SCALE) $atrad(scale)"
    flush stdout
    # update all parameters
    # first Chemical connectivity factor
    if { $atrad(covf_old) != $atrad(covf) } {
	xc_newvalue .mesa $mody(L_COV_SCALE) $atrad(covf)
    }
    # spacefill/ball scale factor
    if { $atrad(scale_old) != $atrad(scale) } {
	xc_newvalue .mesa $mody(L_ATRAD_SCALE) $atrad(scale)
    }

    # is atomname ok
    #set atrad(nat) [Aname2Nat $atrad(atom)]
    #if ![regexp {[0-9]+} $atrad(nat)] {
    #	 tk_dialog .atradup "ERROR" $atrad(nat) error 0 OK
    #	 return
    #}
    #xc_newvalue .mesa $mody(L_ATRAD_ONE) $atrad(nat) $atrad(rad)
    
    # now we can exit
    if { [winfo exists $top] } { 
	#grab release $top
	destroy $top
    }
}


proc AtomRadSelect { w y } {
    global atrad mody

    # $w is the name of listbox widget who contains groups
    # $y is the vertical position of "selection"
    puts stdout "$w select anchor [$w nearest $y]"
    $w select anchor [$w nearest $y]
    set nline [$w curselection]
    set atom [$w get $nline]

    # now we must purify atom variable, because it's sometning like "  8:  O"
    # we will take a number instead of name (communication with xcrys)
    regexp {[0-9]+} $atom atrad(nat)
    
    # to update entry we also need atom name
    regexp {[A-Za-z]+} $atom atrad(atom)

    AtRadLoad

    return
}


proc AtRadDef {var const} {
    upvar #0 $var varn
    global mody 

    set varn [xc_getdefault $mody($const)]

    return
}

proc AtRadClr {var w} {
    upvar #0 $var varn
    
    set varn ""
    focus $w
    return
}
