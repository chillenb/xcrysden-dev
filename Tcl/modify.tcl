#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/modify.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#===========================
# OPTION: NumCellDrawn
#===========================
proc NumCellDrawn {} {
    global species nxdir nydir nzdir nxdirold nydirold nzdirold \
	     NumCellDrawn periodic XCState myParam

    if { $periodic(dim) == 0 } { return }

    #just in any case
    if { [info exists NumCellDrawn(done)] } { unset NumCellDrawn(done) }

    # this option is only for $species != molecule
    #if { $periodic(dim) == 0 } {
    #	tk_dialog .nmdwarning "Warning" "WARNING:\"Number of Units Drawn\" \
    #		option can not be applied for MOLECULES/CLUSTERS !!!" \
    #		warning 0 OK
    #	return
    #}

    # if $species is not defined (whan we are using viewer just for viewing)
    if { [info exist species] == 0 && [xcIsActive properties] == 0 } { return }

    if [winfo exists .ncd] { return }
    set t [xcToplevel .ncd "Modify Number of Units Drawn" "Modify" \
	    . 100 50 1]
    
    # some decorations
    set f1 [frame $t.f1 -relief groove -bd 2]
    pack $f1 -side top -padx 7 -pady 7 -ipadx 7 -ipady 7 -fill x
    
    set m [message $f1.msg -justify center -width 200 -text \
	    "Specify number of cells to be drawn for each dimension:"]
    pack $m -expand 1 -padx 5 -pady 5

    if [info exists nxdir] {set nxdirold $nxdir}
    if [info exists nydir] {set nydirold $nydir}
    if [info exists nzdir] {set nzdirold $nzdir}

    puts stdout "nxdir,nydir,nzdir> $nxdir $nydir $nzdir"
    # display entries for specifying nxdir,nydir,nzdir
    if { $periodic(dim) == 3 } {
	if { [info exists myParam(CRYSTAL_MAXCELL)] } {
	    set max $myParam(CRYSTAL_MAXCELL)
	} else {
	    set max  10
	}
	set dirs {x y z}
    } elseif { $periodic(dim) == 2 } {
	if { [info exists myParam(SLAB_MAXCELL)] } {
	    set max $myParam(SLAB_MAXCELL)
	} else {
	    set max  20
	}
	set dirs {x y}
    } elseif { $periodic(dim) == 1 } {
	if { [info exists myParam(POLYMER_MAXCELL)] } {
	    set max $myParam(POLYMER_MAXCELL)
	} else {
	    set max  50
	}
	set dirs x
    }
    set res  [expr {int($max/10)}]
    foreach i $dirs {
	set sc [scale $f1.sc$i -from 1 -to $max \
		    -length       200 \
		    -variable     n${i}dirold \
		    -orient       horizontal \
		    -label        "In [string toupper $i]-dir:" \
		    -tickinterval $res \
		    -digits       1 \
		    -resolution   1 \
		    -showvalue    true \
		    -width        10]
	pack $sc -side top -fill x -expand 1 -padx 5 -pady 5
    }

    set okf [frame $t.okf -relief sunken -bd 2]
    set ok  [button $okf.ok -text OK -command [list NumCellDrawnOK $t]]
    set upd [button $t.upd  -text Update \
	    -command [list NumCellDrawnOK $t update]]
    set can [button $t.can -text Cancel -command [list NumCellDrawnCan $t]]
    pack $okf -side left -expand 1 -padx 5m -pady 2m 
    pack $upd -side left -expand 1 -padx 5m -pady 2m -ipadx 1m -ipady 1m
    pack $can -side left -expand 1 -padx 5m -pady 2m -ipadx 1m -ipady 1m
    pack $ok -padx 1m -pady 1m -ipadx 1m -ipady 1m
}

 
proc NumCellDrawnOK {t {update 0}} {
    global nxdir nydir nzdir nxdirold nydirold nzdirold \
	    species mode  NumCellDrawn err system periodic
    
    # for $species == molecule it is inpossible to come here !!!
    # ndirx must be specified for POLYMERS, SLABS, CRYSTALS
    check_var {{nxdirold posint}} .ncd.f1.frame.entry1
    # check_var return $err; if $err --> mistake -> RETURN
    if $err { return }

    if { $periodic(dim) > 1 } {
	check_var {{nydirold posint}} .ncd.f1.frame.entry2
    }
    if $err { return }

    if { $periodic(dim) == 3 } {
	check_var {{nzdirold posint}} .ncd.f1.frame.entry3
    }
    if $err { return }

    # if we come so far everything is good, update a structure
    if { $update == 0 } {
	destroy $t
    }
    # check for AdvGeom options; 
    # in future this should be replaced with $periodic(dim) !!!
    if { [info exists AdvGeom(slab)] } { set nzdirold 0 }
    # if OK/UPDATE button was pressed, than set nxdir/nydir/nzdir values
    if { [info exists nxdirold] } { set nxdir $nxdirold }
    if { [info exists nydirold] } { set nydir $nydirold }
    if { [info exists nzdirold] } { set nzdir $nzdirold }
    
    GenGeomDisplay 1
}


proc NumCellDrawnCan {t} {
    destroy $t
    return
}
    



