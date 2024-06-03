#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/bz.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc BAND_Init {} {

    if { ! [Bz_MakeToplevel] } {
	return 0
    }
    return 1
}


proc Bz_MakeToplevel {{what {}}} {
    global Bz BzOK gengeom periodic


    #---
    # For k-labels lookup: try to determine the Bravais-lattice type
    # from periodic(igroup)
    #
    load_kLabels
    set Bz(lattice_type) [igroup2BravaisLattice $periodic(igroup)]
    #xcDebug -stderr "tk: LATTICE-TYPE: $Bz(lattice_type)"
    #kLabels_Note
    #---

    set Bz(what)   $what
    set BzOK(done) 0

    if { ! [info exists BzOK(wien_kpath)] } {
	set Bz(tplw) [xcToplevel [WidgetName] "Band Path Selection" \
			  "Path Selection" . 100 100 1]
    } else {
	# k-path selection for WIENXX program
	update
	set dir [file tail $BzOK(wien_dir)]
	set Bz(tplw) [xcToplevel [WidgetName] \
			  "*** XCrySDen *** K-path selection: $dir" \
			  "Path Selection" . 290 000 1]
    }
    catch { grab $Bz(tplw) }
    wm resizable $Bz(tplw) 0 0
    
    set f1 [frame $Bz(tplw).f1 -highlightthickness 0]
    set f2 [frame $Bz(tplw).f2 -highlightthickness 0]
    pack $f1 $f2 -side top -fill x

    ############
    # FRAME #1 #
    ############
    set Bz(primBZbutton) [button $f1.b1 \
	    -text "Primitive Brillouin Zone" \
	    -bd 3 \
	    -highlightthickness 0 \
	    -command [list Bz_ShowBZ prim]]
    set Bz(convBZbutton) [button $f1.b2 \
	    -text "Conventional Brillouin Zone" \
	    -bd 1 \
	    -highlightthickness 0 \
	    -command [list Bz_ShowBZ conv]]
    pack $Bz(primBZbutton) $Bz(convBZbutton) -side left -fill y

    if { $periodic(dim) < 3 } {
	$Bz(convBZbutton) configure -state disabled
    }

    ############
    # FRAME #2 #
    ############
    set f21 [frame $f2.1]
    set f22 [frame $f2.2]
    set Bz(primBZframe) $f21
    set Bz(convBZframe) $f22
    #
    # check if primitive or conventional BZ really exists
    #
    foreach type {prim conv} {
	if { [xc_bz exists $type] } {
	    Bz_RenderBZ $type
	} else {
	    $Bz(${type}BZbutton) config -state disabled
	}
    }
    
    if { [xc_bz exists prim] } {
	# primitive BZ is default
	pack $f21 -fill both
	set Bz(rendered) prim
    } elseif { [xc_bz exists conv] && $periodic(dim) == 3 } {
	# if there is no primitive BZ, show conventional
	Bz_ShowBZ conv
    } else {
	ErrorDialog "ERROR" "Something has gone wrong. I can't find any Brillouin Zone. It is either bug in the program or you've done something wrong"
	destroy $Bz(tplw)
	return 0
    }

    tkwait variable BzOK(done)
    return $BzOK(done)
}


proc Bz_ShowBZ {type} {
    global Bz

    # maybe already "active" button was pressed
    if { $Bz(rendered) == $type } {
	return
    }
    if { $type == "prim" } {
	pack forget $Bz(convBZframe)
	pack $Bz(primBZframe)
	set Bz(rendered) prim
	$Bz(primBZbutton) configure -bd 3
	$Bz(convBZbutton) configure -bd 1
    } else {
	pack forget $Bz(primBZframe)
	pack $Bz(convBZframe)
	set Bz(rendered) conv
	$Bz(primBZbutton) configure -bd 1
	$Bz(convBZbutton) configure -bd 3
    }
}


proc Bz_RenderBZ {type} {
    global Bz BzOK

    set fl [frame $Bz(${type}BZframe).l]
    set fr [frame $Bz(${type}BZframe).r]
    pack $fl $fr -side left -fill y -expand 1

    ##############
    # frame: $fl #
    ##############
    frame $fl.df -relief raised -bd 3
    pack $fl.df -ipadx 3 -ipady 3
    set can [xc_bz init $type $fl.can \
	    {-width 15c -height 15c -bg "#ffffff"}]
    pack $can -side top -fill both -expand 1 -in $fl.df

    # clear previous k-point coordinates
    
    if { [info exists  Bz($can,coor1)] } {
        array unset Bz $can,coor*
    }
    

    ### some buttons for rotations
    set f [frame $can.f]
    $can create window 1 1 \
	    -anchor nw \
	    -window $f

    set text "Primitive Brillouin Zone"
    if { $type == "conv" } {
	set text "Conventional Brillouin Zone (use only for curiosity)"
    }
    $can create text 7.5c 0.2c \
	    -text $text \
	    -anchor n

    set Bz($can,rot) 5
    set rxp [button $f.rxp -image rotXplus\
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can +x]]
    set rxm [button $f.rxm -image rotXmin \
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can -x]]
    set ryp [button $f.ryp -image rotYplus \
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can +y]]
    set rym [button $f.rym -image rotYmin \
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can -y]]
    set rzp [button $f.rzp -image rotZplus \
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can +z]]
    set rzm [button $f.rzm -image rotZmin \
	    -highlightthickness 0 \
	    -command [list Bz_ManualRotate $can -z]]
    set zup [button $f.zup -image zoomUp \
	    -highlightthickness 0 \
	    -command [list Bz_ManualZoom $can  0.05]]
    set zdn [button $f.zdn -image zoomDown \
	    -highlightthickness 0 \
	    -command [list Bz_ManualZoom $can -0.05]]
    set prn [button $f.prn -image printer \
	    -highlightthickness 0 \
	    -command [list xcPrintCanvas $can]]
    pack $rxp $rxm $ryp $rym $rzp $rzm $zup $zdn $prn -side top

    set Bz(${type}can)         $can
    set Bz($can,npoint)        [xc_bz get $can npoint]
    set Bz($can,npoly)         [xc_bz get $can npoly]
    set Bz($can,nselected)     [xc_bz get $can nselected]
    set Bz($can,last_selected) -1
    set Bz($can,B1down)        0
    set Bz($can,motionB1down)  0
    set Bz($can,state_points)  1
    
    xc_bz viewport $can
    xc_bz render $can
    
    bind $can <Configure>        [list Bz_ViewPort %W]
    bind $can <ButtonPress>      [list Bz_ButtonPressed %W %x %y]
    bind $can <Motion>           [list Bz_Motion %W %x %y]
    bind $can <B1-ButtonRelease> [list Bz_ButtonReleased %W]
    bind $can <Destroy>          [list Bz_BindDestroy %W]
    bind $can <Button-4>         [list Bz_MouseWheel %W +]
    bind $can <Button-5>         [list Bz_MouseWheel %W -]

    #$can bind point   <1> [list Bz_PointButtonPressed %W %x %y]
    #$can bind polygon <1> [list Bz_PolyButtonPressed %W %x %y]

    set bf [frame $fl.fb]
    pack $bf -side bottom -fill x
    
    checkbutton $bf.r1 \
	    -relief raised \
	    -bd 2 \
	    -variable Bz($can,state_points) \
	    -anchor w \
	    -text "Display Special Points" \
	    -command [list Bz_SetState $can points]

    checkbutton $bf.r2 \
	    -relief raised \
	    -bd 2 \
	    -variable Bz($can,state_vectors) \
	    -anchor w \
	    -text "Display Reciprocal Vectors" \
	    -command [list Bz_SetState $can vectors]
    
    #checkbutton $bf.r3 \
	    #	 -text "Display Symbols" \
	    #	 -relief raised \
	    #	 -bd 2 \
	    #	 -state disabled
    pack $bf.r1 $bf.r2 -side left -fill x -expand 1
    
    ##############
    # frame: $fr #
    ##############
    frame $fr.f
    set h 0
    set b1 [button $fr.f.b1 \
	    -highlightthickness 1 \
	    -text "Delete Last\nSelected Point" \
	    -height 2 \
	    -command [list Bz_DeletePoints last $can] \
	    -width 11]
    set bhlt [$b1 cget -highlightthickness]
    incr h [winfo reqheight $b1]
    xcDebug "tcl: button-width: [winfo reqheight $b1]"
    set b2 [button $fr.f.b2 \
	    -highlightthickness 1 \
	    -text "Delete All\nSelected Points" \
	    -height 2 \
	    -command [list Bz_DeletePoints all $can] \
	    -width 11]    
    # make the font size 10
    set Bz(font) [$b2 cget -font]
    set Bz(font) [ModifyFont $Bz(font) $b1 -size 10]
    #$b1 config -font $Bz(font)
    #$b2 config -font $Bz(font)
    set f1 [frame $fr.f1 -relief raised -bd 2 \
	    -highlightthickness 1]
    set l1 [label $f1.l1 \
	    -width 18 \
	    -anchor w \
	    -height 1 \
	    -text "Rotation Step:"]   
    incr h [winfo reqheight $l1]
    set e1 [entry $f1.e1 \
	    -textvariable Bz($can,rot) \
	    -bd 1 \
	    -width 3]

    set f2 [frame $fr.f2 -relief raised -bd 2 \
	    -highlightthickness 1]
    set l2 [label $f2.l2 \
	    -width 18 \
	    -anchor w \
	    -height 1 \
	    -text "# of Selected Points:"]   
    incr h [winfo reqheight $l2]
    set e2 [entry $f2.e2 \
	    -textvariable Bz($can,nselected) \
	    -bd 1 \
	    -width 3]

    #set f2a [frame $fr.f2a -relief raised -bd 2\
    #	     -highlightthickness 1]
    #set l2a [label $f2a.l2 \
    #	     -width 18 \
    #	     -anchor w \
    #	     -height 2 \
    #	     -justify left \
    #	     -text "Total # of k-points\nalong the path:"]   
    #incr h [winfo reqheight $l2a]
    #set e2a [entry $f2a.e2 \
    #	     -textvariable Bz($can,nK) \
    #	     -bd 1 \
    #	     -width 3]
    #set Bz($can,nK_entry) $e2a

    pack $fr.f $f1 $f2 -side top -fill x
    pack $b1 $b2 -side left -fill x -expand 1
    pack $l1 $l2 -side left
    pack $e1 $e2 -side left -pady 2 -padx 2 -expand 1 -fill x
    incr h 4
    set efont [$e1 cget -font]
    set Bz(efont) [ModifyFont $efont $b1 -size 10]
    
    #
    # SCROLLED-CANVAS
    #
    set f3 [frame $fr.f3 -relief raised -bd 2 \
	    -highlightthickness 1]
    pack $f3 -fill both
    set Bz($can,scrollcan) [canvas $f3.can \
	    -width 200 \
	    -yscrollcommand  [list $f3.yscroll set]]
    set scb [scrollbar $f3.yscroll \
	    -orient vertical -command [list $Bz($can,scrollcan) yview]]
    pack $scb -side right -fill y
    pack $Bz($can,scrollcan) -side left -fill y -expand true -padx 2

    # create FRAME to hold everything
    set f [frame $Bz($can,scrollcan).f -bd 0]
    set Bz($can,frame) $f
    $Bz($can,scrollcan) create window 1 1 -anchor nw -window $f

    # now create the first line and compute its width & height
    set cf1 [frame $f.f0]
    set cl1 [label $cf1.l1 -text "#" \
	    -width 2 \
	    -height 1 \
	    -anchor w \
	    -relief sunken -bd 1 \
	    -font $Bz(font)]
    set bg [$cl1 cget -bg]
    # format 6-2-6-2-6 -> width==22
    set Bz(text1) "  reciprocal coordinates"
    set ce1 [entry $cf1.e1 -textvariable Bz(text1) \
	    -bg $bg \
	    -width 22 \
	    -relief sunken -bd 1 \
	    -font $Bz(efont)]
    set Bz(text2) "label"
    set ce2 [entry $cf1.e2 -textvariable Bz(text2) \
	    -bg $bg \
	    -width 7 \
	    -relief sunken -bd 1 \
	    -font $Bz(efont)]
    pack $cf1 -side top -padx 2
    pack $cl1 $ce1 $ce2 -side left
    $ce1 config -state disabled
    $ce2 config -state disabled
    set lhlt [$cl1 cget -highlightthickness]
    set ehlt [$ce1 cget -highlightthickness]
    
    # compute the desired canvas width & height
    set w [expr [winfo reqwidth  $ce1] + [winfo reqwidth  $ce2] + \
	    [winfo reqwidth $cl1] + 2 * $lhlt + 2 * $ehlt]
    set h [expr [winfo fpixels . 15c] - (1.5 * $h)]
    set Bz($can,width) $w
    $Bz($can,scrollcan) config -width $w -height $h
    # compute how many window items can be rendered
    set item_height [expr [winfo reqheight $cl1] + 2]
    set Bz($can,item_height) $item_height
    set Bz($can,n_item) [expr int($h / $item_height)]
    xcDebug "tcl: n_item: $Bz($can,n_item)"
    for {set i 1} {$i <= $Bz($can,n_item) - 1} {incr i} {
	set cf [frame $Bz($can,frame).f$i]
	set cl [label $cf.l -text $i \
		-width 2 \
		-height 1 \
		-anchor w \
		-relief sunken -bd 1 \
		-font $Bz(font)]
	# format 6-2-6-2-6 -> width==22
	set ce1 [entry $cf.e1 -textvariable Bz($can,coor$i) \
		-width 22 \
		-relief sunken -bd 1 \
		-font $Bz(efont) \
		-state disabled]
	set Bz($can,coor_entry$i) $ce1 
	set ce2 [entry $cf.e2 -textvariable Bz($can,label$i) \
		-width 7 \
		-relief sunken -bd 1 \
		-font $Bz(efont) \
		-state normal]
	pack $cf -side top -padx 2
	pack $cl $ce1 $ce2 -side left
    }

    set f4 [frame $fr.f4 -relief raised -bd 2 -highlightthickness 1]
    pack $f4 -side top -fill both -expand 1
    set ok [DefaultButton $f4.ok -text "OK" -command [list Bz_OK $can $type]]
    set can [button $f4.can -text "Cancel" \
	    -command [list CancelProc $Bz(tplw) BzOK(done)]]
    pack $ok $can -side left -expand 1
}


proc Bz_DeletePoints {what can} {
    global Bz
    
    xcDebug "COMMAND:: Bz_DeletePoints"
    if { $what == "last" } {
	set Bz($can,coor$Bz($can,nselected)) {}	
	set Bz($can,label$Bz($can,nselected)) {}	
	set Bz($can,nselected) [xc_bz deselect $can]
	set Bz($can,last_selected) [xc_bz get $can last_selected]
    } elseif { $what == "all" } {
	for {set i 1} {$i <= $Bz($can,nselected)} {incr i} {
	    set Bz($can,coor$i) {}
	    set Bz($can,label$i) {}
	}
	set Bz($can,nselected) [xc_bz deselectall $can]
	set Bz($can,last_selected) [xc_bz get $can last_selected]	
    }
}


#proc Bz_PolyButtonPressed {can x y} {
#
#    set id [$can find closest $x $y]
#    xcDebug "Id: $id"
#    for {set i 0} {$i < 12} {incr i} {
#	set pid [$can find withtag p$i]
#	if { $id == $pid } {
#	    xcDebug "LINE:: polygon #: $i"
#	}
#    }
#}


#proc Bz_PointButtonPressed {can x y} {
#    global Bz
#    
#    set id [$can find closest $x $y]
#    xcDebug "Id: $id"
#    for {set i 0} {$i < $Bz($can,npoint)} {incr i} {
#	set pid [$can find withtag pt$i]
#	if { $id == $pid } {
#	    xcDebug "tcl: Selected point #: $i"
#	    xc_bz select $can $i
#	}
#    }
#}


proc Bz_ViewPort can {

    xcDebug "COMMAND:: Bz_ViewPort"
    xc_bz viewport $can
    xc_bz render   $can
}


proc Bz_ButtonPressed {can x y} {
    global Bz

    xcDebug "COMMAND:: Bz_ButtonPressed"

    # debuging
    xcDebug "LINE:: polyID [$can find overlapping \
    	    [expr $x - 3] [expr $y - 3]  [expr $x + 3] [expr $y + 3]]"
    xcDebug "LINE:: STACKING ORDER == [$can find withtag polygon]\n"

    # maybe a point was selected
    set ps  [$can find enclosed \
    	    [expr $x - 10] [expr $y - 10]  [expr $x + 10] [expr $y + 10]]    
    set nps [expr [llength $ps] - 1]
    
    # upper items are at the end of ps list, so:
    for {set i $nps} {$i > -1} {incr i -1} {
    	 set id [lindex $ps $i]
    	 for {set j 0} {$j < $Bz($can,npoint)} {incr j} {
    	     set pid [$can find withtag pt$j]
    	     if { $id == $pid } {
    		 if { $j == $Bz($can,last_selected) } {
    		     # deselect
    		     xcDebug "tcl: Deselected point #: $j"
    		     set Bz($can,coor$Bz($can,nselected)) {}
		     set Bz($can,label$Bz($can,nselected)) {}
    		     set Bz($can,nselected) [xc_bz deselect $can]
    		     set Bz($can,last_selected) [xc_bz get $can last_selected]
    		     xcDebug "tcl:      last_selected: $Bz($can,last_selected)"
    		 } else {
    		     # select
    		     xcDebug "tcl:   Selected point #: $j"
    		     set list [xc_bz select $can $j]
    		     set Bz($can,nselected) [lindex $list 0]

		     #
		     # assigning the coordinate of selected k-point
		     #
                     set Bz($can,coor$Bz($can,nselected)) ""
    		     foreach item [lrange $list 1 end] {
    			 append Bz($can,coor$Bz($can,nselected)) \
				 [format {% 1.5f} $item]
			 append Bz($can,coor$Bz($can,nselected)) { }
    		     }
		     #
		     # For k-labels lookup: try to set the K-points label
		     #
		     set Bz($can,label$Bz($can,nselected)) [eval {getKLabel $Bz(lattice_type)} $Bz($can,coor$Bz($can,nselected))]

		     # number of selected k-points
    		     set Bz($can,last_selected) $j

		     # take care if nselected is greater than number of entries
		     if { $Bz($can,nselected) > $Bz($can,n_item) - 1 } {
			 # make new entries that will hold coords & label
			 set cf [frame $Bz($can,frame).f$Bz($can,nselected)]
			 set cl [label $cf.l -text $Bz($can,nselected) \
				 -width 2 \
				 -height 1 \
				 -anchor w \
				 -relief sunken -bd 1 \
				 -font $Bz(font)]
			 # format 6-2-6-2-6 -> width==22
			 set ce1 [entry $cf.e1 \
				 -textvariable \
				 Bz($can,coor$Bz($can,nselected)) \
				 -width 22 \
				 -relief sunken -bd 1 \
				 -font $Bz(efont) \
				 -state normal]
			 set ce2 [entry $cf.e2 \
				 -textvariable \
				 Bz($can,label$Bz($can,nselected)) \
				 -width 4 \
				 -relief sunken -bd 1 \
				 -font $Bz(efont) \
				 -state normal]
			 pack $cf -side top -padx 2
			 pack $cl $ce1 $ce2 -side left
			 incr Bz($can,n_item)
			 set h [expr ($Bz($can,nselected) + 1) * \
				 $Bz($can,item_height) + 5]
			 $Bz($can,scrollcan) config \
				 -scrollregion "0 0  $Bz($can,width) $h"
		     }
    		 }
    		 set i -1
    	     }
    	 }
    }
    
    set Bz($can,B1down) 1
}


proc Bz_ButtonReleased can {
    global Bz

    xcDebug "COMMAND:: Bz_ButtonReleassed"
    set Bz($can,B1down) 0
    set Bz($can,motionB1down) 0
}
 
   
proc Bz_Motion {can x y} {
    global Bz

    xcDebug "COMMAND:: Bz_Motion"

    if $Bz($can,B1down) {
	if !$Bz($can,motionB1down) {
	    set Bz($can,motionB1down) 1
	} else {
	    xc_bz rotate $can \
		    [expr $x - $Bz($can,oldX)] [expr $y - $Bz($can,oldY)] 0
	}
	
	set Bz($can,oldX) $x
	set Bz($can,oldY) $y
    }
}


proc Bz_BindDestroy can {
    global Bz

    set ren $Bz(rendered)
    unset Bz
    set Bz(rendered) $ren; #this is used for --wien_kpath
}


proc Bz_SetState {can what} {
    global Bz

    xcDebug "COMMAND:: Bz_SetState $what $Bz($can,state_$what)"
    if { $what == "points" } {
	xc_bz state $can points $Bz($can,state_points)
    } elseif { $what == "vectors" } {
	xc_bz state $can vectors $Bz($can,state_vectors)
    }
}


proc Bz_ManualRotate {can type} {
    global Bz

    xcDebug "COMMAND:: Bz_ManualRotate"
    switch -exact -- $type {
        "+x" { xc_bz degrotate $can -$Bz($can,rot) 0 0 }
	"-x" { xc_bz degrotate $can $Bz($can,rot) 0 0 }
	"+y" { xc_bz degrotate $can 0 -$Bz($can,rot) 0 }
	"-y" { xc_bz degrotate $can 0 $Bz($can,rot) 0 }
	"+z" { xc_bz degrotate $can 0 0 -$Bz($can,rot) }
	"-z" { xc_bz degrotate $can 0 0 $Bz($can,rot) }
    }
}


proc Bz_OK {can type} {
    global Bz prop BzOK fillEntries

    set input {}
    set title {}
    if { $type == "conv" } {
	append input "CONVCELL\n"
    }
    append input "BAND\n"

    # at least two points must have been selected; check that
    if { $Bz($can,nselected) < 2 } {
	tk_dialog [WidgetName] ERROR "ERROR !\nAt least two points must be selected for \"band structure\". Please do so !" error 0 OK
	focus $Bz($can,coor_entry1)
	return
    }
        
    set BzOK(iss) [xc_bz iss $can]
    
    # t.k.: Thu Feb 26 22:17:52 CET 2004
    global wnKP
    set wnKP(npoi) $Bz($can,nselected)
    set wnKP(M)    $BzOK(iss)
    set wnKP(type) $type
    #/

    # COPY significant data from Bz() array before it will destroyed by
    # <Destroy> binding command
    set nselected $Bz($can,nselected)
    for {set i 1} {$i <= $Bz($can,nselected)} {incr i} {
	#
	foreach c $Bz($can,coor$i) {
	    append coorlabel "[format {% 4d} [expr round($c * $BzOK(iss))]] "
	}
	append coorlabel "          $Bz($can,label$i)\n"	

	# t.k:
	set ipol 1
	foreach c $Bz($can,coor$i) {	
	    set wnKP(poi,$i,$ipol) $c
	    incr ipol
	}
	#/
    }
    
    catch { grab release $Bz(tplw) }
    if { $Bz(tplw) != "." } {
	destroy $Bz(tplw)
    }
    
    # maybe user want's to modify something
    set t [xcToplevel [WidgetName] "Band Structure Script" \
	    "BAND" . 100 100 1]
    catch { grab $t }

    set f1 [frame $t.f1 -relief raised -bd 2]    
    set f2 [frame $t.f2 -relief raised -bd 2]
    pack $f1 $f2 -side top -fill both -expand 1

    set text {
	NOTES:
	
	1. Please check if k-point coordinates looks OK. These k-points have integer coordinates, calculated as (kx*M, ky*M, kz*M), where M is an integer multiplier. Actual k-point coordinates are then obtained as: 1/M * (kx*M, ky*M, kz*M) !!!

	2. You can also specify some new k-points manually.
	
	3. The k-path is obtained by connecting the k-points together.
    }

    message $f1.msg -aspect 300 -relief ridge -bd 2 -text $text
    pack $f1.msg -side top -padx 5 -pady 5
    
    FillEntries $f1 {
	"M multiplier was set to:"
	"Total number of k-points along the path:"
    } [list BzOK(iss) BzOK(nK)] 37 8 top left
    
    # t.k: are the any side-effect for setting prop(n_band) to non-number ???
    if { ! [info exists prop(n_band)] } {
	set prop(n_band) XX
	set prop(firstband) 1
	set prop(lastband) 10
    }

    set BzOK(foclist) $fillEntries
    if { ! [info exists BzOK(wien_kpath)] } {
	set fl3 [SelBandIntv $f1]
	set BzOK(varlist) {
	    {BzOK(iss)       posint} 
	    {BzOK(nK)        posint} 
	    {prop(firstband) posint} 
	    {prop(lastband)  posint}
	}
    } else {
	# specify energy range
#	set l1 [label $f1.l1 -text \
#		"Please enter energy interval for BAND structure calculation !!!" \
#		-relief flat -justify left -anchor w]
#
#	pack $l1 -side top -pady 5 -fill x -expand 1
#	FillEntries $f1 {
#	    "Minimum Energy:" 
#	    "Maximum Energy:"
#	} {BzOK(Emin) BzOK(Emax)} 37 8 top left
#	set BzOK(foclist) [concat $BzOK(foclist) $fillEntries]
	set BzOK(varlist) {
	    {BzOK(iss)  posint} 
	    {BzOK(nK)   posint} 
	}
#	    {BzOK(Emin) real} 
#	    {BzOK(Emax) real}

	# this is needed to complete c95_BAND_script
	set prop(firstband) 1
	set prop(lastband)  10
    }

    set f2l [frame $f2.l]
    set f2r [frame $f2.r]
    set f2ll [frame $f2l.l]
    set f2lr [frame $f2l.r]
    set f2llu [frame $f2ll.u]
    set f2llb [frame $f2ll.b]

    pack $f2l -side left -fill y -padx 5 -pady 5
    pack $f2r -side left -fill both -expand 1
    pack $f2ll $f2lr   -side left -fill y
    pack $f2llu $f2llb -side top


    set l1 [text $f2ll.l1 -width 30 -height 1]
    $l1 insert 1.0 [format "%7s %4s %4s %4s  %7s" Format: M*kx M*ky M*kz label]
    $l1 config -state disabled
    $l1 config -bg [$f2l cget -bg]

    set t1_h [expr $nselected + 1] 
    if { $t1_h < 10 } {
	set t1_h 10
    }
    set t1 [text $f2ll.t1 \
		-setgrid true \
		-wrap none \
		-width 30 \
		-height $t1_h \
		-yscrollcommand "$f2lr.sy set"]
    set sy [scrollbar $f2lr.sy -orient vert \
		-command "$f2ll.t1 yview"]
	    
    pack $l1 $t1 -side top -padx 2 -pady 2 -fill y
    pack $sy -side left -padx 2 -pady 2 -fill y

    # this ...:
    $t1  insert 1.0 $coorlabel
    
    set ok [DefaultButton $f2r.ok -text "OK" \
		-command [list Bz_OK_OK $t $type $t1]] 
    set can [button $f2r.can -text "Cancel" \
		 -command [list CancelProc $t BzOK(done)]]
    pack $ok $can -side top -expand 1 -padx 5 -pady 5
}


proc Bz_OK_OK {tplw type coor_label} {
    global BzOK prop properties kpath wnKP

    set coorlabeltext  [$coor_label get 1.0 end]
    xcDebug -stderr "Coorlabel text:\n $coorlabeltext"
    
    # first check the entry-variables
    if { ! [check_var $BzOK(varlist) $BzOK(foclist)] } {
	return
    }

    #
    # assign default labels
    #

    # since we have automatic k-point labeling such a naming might be
    # misleading
    #set labels {A B C D E F G H I J K L M N O P R S T U V Z X Y W a b c d e f g h i j k l m n o p r s t u v z x y w}
    for {set i 1} {$i < 100} {incr i} {
	append labels " K.$i"
    }
    
    set errorText {ERROR: you have mistype the k-point coordinates. Please try again !}

    #
    # parse coorlabeltext and assign coordinates and labels
    # 
    set ith   0
    foreach line [split $coorlabeltext \n] {
	#if { [llength $line] < 3 || [llength $line] > 4 } {
	#    tk_dialog [WidgetName] ERROR  $errorText error 0 OK
	#    return
	#}
	if { [llength $line] == 3 || [llength $line] == 4 } {
	    for {set k 0} {$k < 3} {incr k} {
		set poi($ith,$k) [lindex $line $k]
		if { ! [string is integer -strict $poi($ith,$k)] } {
		    ErrorDialog $errorText
		    return
		}
	    }
	    set lab($ith) [lindex $line 3]
	    if { $lab($ith) == {} } {
		set lab($ith) [lindex $labels $ith]
	    }	    
	    # t.k.
	    set wnKP(label,[expr $ith + 1]) $lab($ith)
	    #/
	    incr ith
	}
    }

    set kpath(point_labels) {}
    for {set i 0} {$i < $ith} {incr i} {
	set ii [expr $i + 1]
	set properties(TICK$ii) $lab($i)
	append title $lab($i)-
	append kpath(point_labels) \
	    [format "%3d  %3d  %3d     %s\n" \
		 $poi($i,0) $poi($i,1) $poi($i,2) $lab($i)]
	#
    }
    set title [string trimright $title -]

    if { ! [info exists prop(newk_script)] } {
	set prop(newk_script) ""
    }
    set prop(c95_BAND_script) $prop(newk_script)
    if { $type == "conv" } {
	set kpath(basis) "CONVENTIONAL"
	append prop(c95_BAND_script) "CONVCELL\n"
    } else {
	set kpath(basis) "PRIMITIVE"
    }

    set prop(NLINE) [expr $ith - 1]
    append prop(c95_BAND_script) "BAND\n$title\n"
    append prop(c95_BAND_script) \
	    "$prop(NLINE) $BzOK(iss) $BzOK(nK) $prop(firstband) $prop(lastband) 1 0\n"
    for {set i 0} {$i < $prop(NLINE)} {incr i} {
	set ii [expr $i + 1]
	append prop(c95_BAND_script) \
		"$poi($i,0) $poi($i,1) $poi($i,2)    $poi($ii,0) $poi($ii,1) $poi($ii,2)\n"
    }
    
    xcDebug -stderr "CRYSTALxx Band Structure Script:"
    xcDebug -stderr "--------------------------------"
    xcDebug -stderr $prop(c95_BAND_script)

    catch { grab release $tplw }

    destroy $tplw
    set BzOK(done) 1
}


proc Bz_ManualZoom {can initstep} {
    global Bz
    xc_bz zoom $can [expr $initstep * $Bz($can,rot)]
}

proc Bz_MouseWheel {can dir} {
    set zoom 0.075
    xc_bz zoom $can ${dir}$zoom
}

#proc Bz_CheckScript can {
#    global Bz
#
#    # CHECK THE FOLLOWING::
#    # 1.) at least two points must have been selected
#    # 2.) check if the coordinates are specified correctly. 
#    
#    # 1.)
#    #if { $Bz($can,nK) == "" } {
#    #	 tk_dialog [WidgetName] ERROR "ERROR !\nYou forget to specify \
#    #		 the \"Bz($can,nK)\" variable. Please do so !" error 0 OK
#    #	 focus $Bz($can,nK_entry)
#    #	 return 0
#    #}
#    #if ![check_var [list [list Bz($can,nK) posint]] $Bz($can,nK_entry)] {
#    #	 return 0
#    #}
#	 
#    # 1.)
#    if { $Bz($can,nselected) < 2 } {
#	 tk_dialog [WidgetName] ERROR "ERROR !\nAt least two points must be selected for \"band structure\". Please do so !" error 0 OK
#	 focus $Bz($can,coor_entry1)
#	 return 0
#    }
#    
#    ## 2.)
#    #for {set i 1} {$i <= $Bz($can,nselected)} {incr i} {
#    #	 if { $Bz($can,coor$i) == "" } {
#    #	     tk_dialog [WidgetName] ERROR "ERROR !\nYou forget to specify \
#    #		     coordinates for point #$i. Please do so !" error 0 OK
#    #	     focus $Bz($can,coor_entry$i)
#    #	     return 0
#    #	 }	
#    #	 set n 0
#    #	 foreach num $Bz($can,coor$i) {
#    #	     if [catch {expr abs($num)}] {
#    #		 dialog .number1 ERROR "ERROR !\nYou have specified a character instead of numbers for coordinates of point #$i. Try again !" error 0 OK
#    #		 focus $Bz($can,coor_entry$i)
#    #		 return 0
#    #	     }
#    #	     incr n
#    #	 }
#    #	 if { $n < 3 } {
#    #	     dialog .number1 ERROR "ERROR !\nYou have badly specified the coordinates of point #$i. Intead of 3 numbers, You specified just $n numbers. Try again !" error 0 OK
#    #	     focus $Bz($can,coor_entry$i)
#    #	     return 0
#    #	 }	
#    #}
#
#    return 1
#}
