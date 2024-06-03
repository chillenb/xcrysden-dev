#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnFS.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# EXECUTED BY BUTTON:  "Generate k-mesh"
#
proc wnGenKMesh {} {
    global wn

    #
    # order of questions in KGEN:
    #
    # 1. "NUMBER OF K-POINTS IN WHOLE CELL:"
    # 2. "Shift of k-mesh allowed. Do you want to shift:" 

    ###########
    cd $wn(dir)
    ###########

    set input    "$wn(fs_nkp)\n"
    append input "0\n"; # shift of k-mesh is not allowed

    xcDebug -debug "KGEN INPUT: $input"

    WriteFile xc_kgen.inp $input w
    catch {exec x kgen < xc_kgen.inp > xc_kgen.out} exit_status

    set out [ReadFile xc_kgen.out]
    dialog [WidgetName] "Notification" \
	"K-MESH GENERATED !!!.\nOutput of kgen:\n\n $out\n\nExit status:\n$exit_status" \
	info 0 Done
}


#
# EXECUTED BY BUTTON:   "Render Fermi Surface"
#
proc wnFSGo {outkgen {spin {}}} {
    global wn system

    SetWatchCursor
    update

    xcDebug -stderr "DEBUG -- args: $outkgen $spin"
    ###########
    cd $wn(dir)
    ###########

    set outkgen [file tail $outkgen]

    # try outkgen
    # try $wn(dir)/$wn(filehead).fs_outputkgen
    # try $wn(dir)/$wn(filehead).outputkgen !!! DANGEROUS
    if { ! [file exists $outkgen] } {
	tk_dialog [WidgetName] ERROR \
	    "Please generate the k-mesh first !!!" error 0 OK
	ResetCursor; return
    }
    
    set ok 0
    foreach line [split [ReadFile $outkgen] \n] {
	if [string match "*k-vectors:*" $line] {
	    set ok 1
	    break
	}
    }
    
    if { ! $ok } {
	ErrorDialog "Please generate the k-mesh first !!!"
	ResetCursor
	return	
    }

    #
    # read output1 file by wn_readbands prog !!!
    #
    append deffile "6,'$wn(filehead).outputba$spin','unknown','formatted',0\n"
    append deffile "7,'$wn(filehead).output1$spin' ,'old',    'formatted',0\n"
    append deffile "8,'$wn(filehead).outputbw$spin','unknown','formatted',0\n"
    
    set wn($spin,fs_bandfile)      $wn(filehead).outputba$spin
    set wn($spin,fs_bandwidthfile) $wn(filehead).outputbw$spin
    
    WriteFile band.def $deffile w
    if { [catch {exec $system(BINDIR)/wn_readbands band.def} error_msg] } {
	ErrorDialog "while executing wn_readbands program !!!" $error_msg
	ResetCursor
	return
    }
    
    #
    # execute program wn_readbakgen
    #
    set deffile {}
    append deffile "7, '$wn(filehead).outputba$spin','unknown','formatted',0\n"
    append deffile "8, '$outkgen'                   ,'old'    ,'formatted',0\n"
    append deffile "10,'$wn(filehead).outputfs$spin','unknown','formatted',0\n"

    WriteFile bakgen.def $deffile w
    if { [xcCatchExec $system(BINDIR)/wn_readbakgen bakgen.def] } {
	ResetCursor
	return
    }
    set wn($spin,fs_fsfile) $wn(filehead).outputfs$spin    
    file rename -force -- $wn($spin,fs_fsfile) $system(SCRDIR)/$wn(filehead)
    set wn($spin,fs_fsfile) $system(SCRDIR)/$wn(filehead)/$wn($spin,fs_fsfile)

    ResetCursor
    update

    ##############################################################
    #
    # BAND SELECTION
    #
    ##############################################################
    
    set wn(fs_Efermi) 0.0
    catch {set wn(fs_Efermi) \
	       [exec grep :FER $wn(filehead).output2$spin | \
		    tail -1 | awk "{print \$NF}"]}
    
    OneEntryToplevel [WidgetName] "Fermi Energy" "Ferm Energy" \
	"Specify the Fermi Energy:" 15 wn(fs_Efermi) float 300 20

    set ncol 1
    set nrow 500; # maximum number of bands to display in a single row
		  # for selection checkbuttons
    
    #
    # display the bandwidths in a Text widget
    #
    set text [FS_displayBandWidths $wn($spin,fs_bandwidthfile) $spin]    
    wm geometry $text +0-0
    raise $text 
    update; update idletasks

    SetWatchCursor

    set cw [DisplayUpdateWidget "Please wait" "Please wait: parsing bands data"]

    # read number of bands & band-widhts
    
    set nl 0
    foreach lin [split [ReadFile $wn($spin,fs_bandwidthfile)] \n] {
	if { $nl == 0 } {
	    set wn($spin,nbands) $lin
	} else {
	    if { $nl <= $wn($spin,nbands) } {
		set wn($spin,$nl,minE) [lindex $lin 1]
		set wn($spin,$nl,maxE) [lindex $lin 2]
	    }
	}
	incr nl
    }
    
    if { $wn($spin,nbands) < $nrow } {
	#
	# make a band-width graph !!!
	#
	set xlabel "Band Widths" 
	if { $spin != {} } {
	    append xlabel " (spin type: $spin)"
	}
	GraphInit
	grapher_BARGraph $wn($spin,fs_bandwidthfile) \
	    -Xtitle $xlabel \
	    -Ytitle "E / Ry" \
	    -Yline  $wn(fs_Efermi) \
	    -Yline_text Ef
	set graph [Grapher BARGraph]
	wm geometry $graph -0-0
	raise $graph
    }

    SetWatchCursor
    
    #
    # select bands window
    #

    set t [xcToplevel [WidgetName] "Select bands" "Select Bands" .fs_init 0 0 1]
    wm geometry $t -0+0
    raise $t    
    tkwait visibility $t

    SetWatchCursor
    
    label $t.l \
	-text "Select bands for Fermi Surface drawing:" \
	-relief ridge -bd 2
    pack $t.l -side top -expand 1 -fill x -padx 2m -pady 3m \
	-ipadx 2m -ipady 2m
    
    #
    # we should make a scrolled window
    #
    # CANVAS & SCROLLBAR in CANVAS
    set scroll_frame [frame $t.f -relief sunken -bd 1]
    set xroll_frame   [frame $t.x -relief sunken -bd 1]
    pack $scroll_frame -side top -expand true  -fill y -padx 5    
    pack $xroll_frame  -side top -expand false -fill x
    set c [canvas $scroll_frame.canv \
	       -yscrollcommand [list $scroll_frame.yscroll set] \
	       -xscrollcommand [list $xroll_frame.xscroll set]]
    set scb  [scrollbar $scroll_frame.yscroll \
		 -orient vertical -command [list $c yview]]
    set scx  [scrollbar $xroll_frame.xscroll \
		 -orient horizontal -command [list $c xview]]
    pack $scb -side right -fill y
    pack $c   -side left  -fill both -expand true
    pack $scx -side bottom -fill x
    
    # create FRAME to hold all checkbuttons
    set f [frame $c.f -bd 0]
    $c create window 0 0 -anchor nw -window $f -tags frame

    bind $c <4> { %W yview scroll -5 units }
    bind $c <5> { %W yview scroll +5 units }
    bind $f <4> " $c yview scroll -5 units "
    bind $f <5> " $c yview scroll +5 units "

    # CHECKBUTTONS
    
    if { $wn($spin,nbands) < $nrow } {
	
	for {set i 1} {$i <= $wn($spin,nbands)} {incr i} {
	    set wn($spin,$i,band_selected) 0
	    set cb [checkbutton [WidgetName $f] -text "Band number: $i" \
			-variable wn($spin,$i,band_selected) -relief ridge -bd 2]
	    pack $cb -padx 2m -pady 1 -fill both -expand 1
	}
	#
	set child [lindex [pack slaves $f] 0]
	#
    } else {
	# arrange checkbuttuon in matrix with $nrow rows
	set ncol [expr $wn($spin,nbands) / $nrow]
	if { $ncol*$nrow < $wn($spin,nbands) } {
	    incr ncol
	}

	set nb 1
	for {set ir 0} {$ir < $nrow} {incr ir} {
	    for {set ic 0} {$ic < $ncol} {incr ic} {
		if { $nb <= $wn($spin,nbands) } {
		    set wn($spin,$nb,band_selected) 0
		    set cb [checkbutton [WidgetName $f] -text "Band #: $nb" \
				-variable wn($spin,$nb,band_selected) -relief ridge -bd 2]
		    grid $cb -column $ic -row $ir -padx 1 -pady 1

		    bind $cb <4> " $c yview scroll -5 units "
		    bind $cb <5> " $c yview scroll +5 units "

		    incr nb
		}
	    }
	}
	#
	set child [lindex [grid slaves $f] 0]
	#
    }
    
    # make correct DISPLAY
    
    tkwait visibility $child
    update
    
    set width  [winfo width $f]
    set height [winfo height $f]

    set wsc [expr min([winfo screenwidth  $t],1920) - 100]
    set hsc [expr min([winfo screenheight $t],1080) - 100]

    set w1 [expr $width  / $ncol]
    set h1 [expr $height / min($nrow,$wn($spin,nbands))]

    set nw [expr min(10,$wsc/$w1+1)]
    set nh [expr min(20,$hsc/$h1+1)]

    #puts stderr "wsc=$wsc, w1=$w1, nw=$nw"
    #puts stderr "hsc=$hsc, h1=$h1, nh=$nh"
    #
    #if { [info exists ncol] } {
    #	puts stderr "FS: width  = $width  --> [expr min(10,$ncol)*($width  / $ncol)] ([winfo reqwidth  $f])"
    #}
    #puts stderr "FS: height = $height --> [expr min(20,$nrow)*($height / $wn($spin,nbands))] ([winfo reqheight $f])"
    #puts stderr "FS: nbands = $wn($spin,nbands)"
    
    if { $wn($spin,nbands) < $nh } {
	$c config -width $width -height $height 
    } else {
	if { $wn($spin,nbands) < $nrow } {
	    $c config \
		-width $width -height [expr $nh*$h1] -scrollregion "0 0 $width $height"
	} else {
	    $c config \
		-width  [expr $nw*$w1] \
		-height [expr $nh*$height / $nrow] \
		-scrollregion "0 0 $width $height"	    
	}
    }
    
    #
    # press the "Selected" button when done
    #    
    set b [button [WidgetName $t] -text "Selected" \
	       -command [list wnToFS $t $spin]]
    pack $b -side top -expand 1 -fill x -padx 2m -pady 3m    

    if { $spin != {} } {
	set l [label [WidgetName $t] -text " Spin type: $spin" \
		   -relief ridge -bd 4 -anchor w]
	pack $l -side top -expand 1 -fill x -padx 1m -pady 1m \
	    -ipadx 2m -ipady 2m
    }

    ResetCursor
    destroy $cw
}


proc wnToFS {t {spin {}}} {
    global wn fs xcMisc

    CancelProc $t

    set fs(titlefile)      $wn(filehead)
    set fs($spin,nbands)   $wn($spin,nbands)
    set fs(Efermi)         $wn(fs_Efermi)
    for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	set fs($spin,$i,band_selected) $wn($spin,$i,band_selected)
	set fs($spin,$i,minE)          $wn($spin,$i,minE) 
	set fs($spin,$i,maxE)          $wn($spin,$i,maxE)
	set fs($spin,$i,isolevel)      $fs(Efermi)
    }

    # lets read (band)XSF file; xc_readbandXSF return info structure
    # which look like:
    #
    # 1 { grid_index 3D grid_ident grid_nband 1 {subgrid0_ident}}
    set sinfo [xc_readbandXSF $wn($spin,fs_fsfile)]
    xcDebug -debug "DEBUG: FS-sinfo: $sinfo"
    # parse $sinfo
    set slist [lindex $sinfo 1]
    set fs($spin,grid_index)    [lindex $slist 0]
    set fs($spin,grid_subindex) [expr [lindex $slist 4] - 1]
    if { [lindex $slist 3] != $wn($spin,nbands) } {
	dialog [WidgetName] "ERROR" \
	    "Mismatch occured while reading FermiSurface File: $wn($spin,fs_fsfile)" error 0 Done
    }
    FS_GoFermi $spin
}

