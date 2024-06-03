# 
# FS_bandSelect --
#
# Select a bands for Fermi-surface plotting among the all the bands from
# small toplevel window. Also a abnd-widths are ploted for orientation
#
proc FS_bandSelect {spin} {
    global fs

    set ncol 1
    set nrow 500; # maximum number of bands to display in a single row
		  # for selection checkbuttons

    
    # specify fermi energy

    OneEntryToplevel [WidgetName] "Fermi Energy" "Ferm Energy" \
	"Specify the Fermi Energy:" 15 fs(Efermi) float 300 20

    #
    # display the bandwidths in a Text widget
    #
    set text [FS_displayBandWidths $fs($spin,bandwidthfile)]
    wm geometry $text +0-0
    raise $text 
    update; update idletasks
    
    SetWatchCursor

    set cw [DisplayUpdateWidget "Please wait" "Please wait: parsing bands data"]

    # read number of bands & band-widhts
    
    set nl 0
    foreach lin [split [ReadFile $fs($spin,bandwidthfile)] \n] {
	if { $nl == 0 } {
	    set fs($spin,nbands) $lin
	} else {
	    if { $nl <= $fs($spin,nbands) } {
		set fs($spin,$nl,minE) [lindex $lin 1]
		set fs($spin,$nl,maxE) [lindex $lin 2]
	    }
	}
	incr nl
    }    
    
    if { $fs($spin,nbands) < $nrow } {
	#
	# make a band-width graph !!!
	#
	set xlabel "Band Widths" 
	if { $spin != {} } {
	    append xlabel " (spin type: $spin)"
	}    
	GraphInit
	grapher_BARGraph $fs($spin,bandwidthfile) \
	    -Xtitle     $xlabel \
	    -Ytitle     "E" \
	    -Yline      $fs(Efermi) \
	    -Yline_text Ef
	set graph [Grapher BARGraph]
	wm geometry $graph +0+0
	raise $graph
    }

    SetWatchCursor
    
    #
    # select bands to plot Fermi surface
    #

    set t [xcToplevel [WidgetName] "Select bands" "Select Bands" . 0 0 0]
    raise $t

    SetWatchCursor
    
    label $t.l \
	-text "Select bands for Fermi Surface drawing:" \
	-relief ridge -bd 2
    pack $t.l -side top -expand 1 -fill x -padx 2m -pady 3m \
	-ipadx 2m -ipady 2m

    
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
    
    if { $fs($spin,nbands) < $nrow } {
	
	for {set i 1} {$i <= $fs($spin,nbands)} {incr i} {
	    set fs($spin,$i,band_selected) 0
	    set cb [checkbutton [WidgetName $f] -text "Band number: $i" \
			-variable fs($spin,$i,band_selected) -relief ridge -bd 2]
	    pack $cb -padx 2m -pady 1 -fill both -expand 1
	}
	#
	set child [lindex [pack slaves $f] 0]
	#
    } else {
	# arrange checkbuttuon in matrix with $nrow rows
	set ncol [expr $fs($spin,nbands) / $nrow]
	if { $ncol*$nrow < $fs($spin,nbands) } {
	    incr ncol
	}

	set nb 1
	for {set ir 0} {$ir < $nrow} {incr ir} {
	    for {set ic 0} {$ic < $ncol} {incr ic} {
		if { $nb <= $fs($spin,nbands) } {
		    set fs($spin,$nb,band_selected) 0
		    set cb [checkbutton [WidgetName $f] -text "Band #: $nb" \
				-variable fs($spin,$nb,band_selected) -relief ridge -bd 2]
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
    set h1 [expr $height / min($nrow,$fs($spin,nbands))]

    set nw [expr min(10,$wsc/$w1+1)]
    set nh [expr min(20,$hsc/$h1+1)]

    #puts stderr "wsc=$wsc, w1=$w1, nw=$nw"
    #puts stderr "hsc=$hsc, h1=$h1, nh=$nh"
    #
    #if { [info exists ncol] } {
    #	puts stderr "FS: width  = $width  --> [expr min(10,$ncol)*($width  / $ncol)] ([winfo reqwidth  $f])"
    #}
    #puts stderr "FS: height = $height --> [expr min(20,$nrow)*($height / $fs($spin,nbands))] ([winfo reqheight $f])"
    #puts stderr "FS: nbands = $fs($spin,nbands)"
    
    if { $fs($spin,nbands) < $nh } {
	$c config -width $width -height $height 
    } else {
	if { $fs($spin,nbands) < $nrow } {
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
	       -command [list FS_bandSelect:selected $t $spin]]
    pack $b -side top -expand 1 -fill x -padx 2m -pady 3m    
    
    if { $spin != {} } {
	set l [label [WidgetName $t] -text " Spin type: $spin" \
		   -relief ridge -bd 4 -anchor w]
	pack $l -side top -expand 1 -fill x -padx 1m -pady 1m \
	    -ipadx 2m -ipady 2m
    }

    ResetCursor
    destroy $cw
    
    tkwait variable fs($spin,selected)
}
proc FS_bandSelect:selected {t spin} {
    global fs

    set fs($spin,selected) 1
    CancelProc $t
}
