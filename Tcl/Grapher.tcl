#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/Grapher.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# ------------------------------------------------------------------#
# grafdata ARRAY definition ----------------------------------------#
# ------------------------------------------------------------------#
#

# grafdata(type) ............ what kind of Graph is rendered (DOSS, BAND, BAR)
# grafdata([1..9],[1..9],[1..9]).... values of function (point,segment,graph)
# grafdata(magnify$g)    .... magnify data in graph $g
# grafdata(X[1..9])      .... X point values 
# grafdata(N_graf)       .... number of graphs
# grafdata(N_segment)    .... number of different line segments
# grafdata(N_point)      .... number of points in each segment
# grafdata(Xmin,$g)      .... first X value
# grafdata(Xmax,$g)      .... last X value
# grafdata(RorigX,$g)    .... relative (inside graph-place) origin coord
# grafdata(RsizeX,$g)    .... relative (inside graph-place) X size
# grafdata(Xoffset,$g)   .... offset before and after last point in X 
#                             directions
# grafdata(dX,$g)        .... X increment
# grafdata(Yline[1..9])  .... special horizontal line that occur at Y
# grafdata(Xline[1..9])  .... special vertical line that occur at X
# grafdata(Xline[1..9]_text)  text that goes with Xline[1..9]  
# grafdata(XYline_textX) .... X coor of X/Yline_text relative to top-right
#                             corner 
# grafdata(N_Yline)      .... number of special Y line
# grafdata(N_Xline)      .... number of special X line
# grafdata(Xtitle)       .... title that appear on Xwindow frame
# grafdata(Xicon)        .... name of icon
# grafdata(X_title)      .... title of X axe
# grafdata(Y_title)      .... title of Y axe
# grafdata(X_title_X)    .... X coordinate of X title
# grafdata(N_text)       .... number of arbitrary texts
# grafdata(text[1..9])   .... arbitrary text, specific for some graph;
#                             this is an array: text1, text2, ....
# grafdata(text[1..9]_X) .... X coord of arbitrary text
# grafdata(text[1..9]_font) . font of arbitrary text
# grafdata(text[1..9]_fontsize) . font-size of arbitrary text
# grafdata(text[1..9]_bbox) . 1/0 -> if bbox is drawn around text
# grafdata(N_MXtick)     .... number of MAJOR X ticks
# grafdata(N_mXtick)     .... number of MINOR X ticks between two major ticks
# grafdata(Xtick[1..9])  .... value of MAJOR X ticks
# grafdata(Xtick_draw)  0/1 whether to draw
# grafdata(Xtick[1..9]_text)  text for Xtick[1..9] tick
# grafdata(Xtick_text_draw) 0/1 whether to draw 
# ===========================================================================
#
# ===========================================================================
# grafsize(bbox,$g)      .... bbox for graph $g
# grafsize(axe_title_offset)  offset of axis's X_title & Y_title
# grafsize(X1axe_width)  .... width of X1 & Y1 axis
# grafsize(X2axe_width)  .... width of X2 & Y2 axis
# grafsize(margin_X1)    .... size of bottom X margin around world
# grafsize(canW)         .... width of graph's canvas
# grafsize(conH)         .... heigth of control frame
# grafsize(Yscroll)      .... Y scroll size of canvas; if there is no 
#                             need for Yscroll --> $grafsize(Yscroll) 
#                             should be 0       
# grafsize(Mtick_size)   ... major tick size
# grafsize(mtick_size)   ... minor tick size
# grafsize(tt_offset)    ... thich_title offset
# grafsize(curve_width)  ... width of curve (for all curves, 
#                            if not specified otherwise)
# grafsize(curve_color)  ... color of curve (for all curves, 
#                            if not specified otherwise)
# grafsize(curve_width_$i) . width of i-th curve
# grafsize(curve_color_$1) . color of i-th curve
# grafsize(Xline[1..9]_fill)
# grafsize(Xline[1..9]_width)
# grafsize(Xline[1..9]_stipple)
# ..... the same for Yline .....
# grafsize(bar[1..9]_attrib) ... attrib -> attributes for BAR (BARGraph)
# ==========================================================================
#
# ==========================================================================
# graf(X1) left   X start of graphplace
# graf(Y1) upper  Y start od graphplace
# graf(X2) right  X start of graphplace
# graf(Y1) bottom Y start od graphplace
# graf(active_graf) which graph is currently selected (in case of multygraph)
# ==========================================================================
# 
# ==========================================================================
# grafselection(allow_move) .......... 1/0 is item currently allowed to move
# grafselection(item_allowed_to_move). 1/0 is item allowed to move by 
#                                      definition
# grafselection(active)     .......... when selection os on
# grafselection(id)         .......... ID of selected item
# grafselection(ButtonPressed) ....... 1/0 is button pressed and held down
# grafselection(lastX)      .......... last x-motion coordinate
# grafselection(lastY)      .......... last y-motion coordinate
# grafselection(exist_position_text).. does position text already exists
# grafselection(position_text) ....... glbal variable for position text
# grafselection(label_window) ........ name of label window, that hold 
#                                      position text
# grafselection(item)       .......... decsriptor for selected item; it is the
#                                      tag of item
# grafselection(it_is_curve) ......... curve was selected
# grafselection(sel_curve_num) ....... number of selected curve 
# grafselection(selcurve_fill) ....... fill atribute of selected curve
# grafselection(selcurve_width) ...... width atribute of selected curve
proc GrapherID {} {
    global grapher

    if ![info exists grapher(id)] {
	set grapher(id) 1
    } else {
	incr grapher(id)
    }

    return $grapher(id)
}


proc CurrentGrapherID {} {
    global grapher
    if ![info exists grapher(id)] {
	return 0
    }
    
    return $grapher(id)
}


proc NextGrapherID {} {
    global grapher
    if ![info exists grapher(id)] {
	return 1
    }

    return [expr $grapher(id) + 1]
}


proc Grapher {graphtype} {
    global grafdata grafsize

    set gID [GrapherID]
    xcDebug "In Grapher"
    if ![info exists grafdata($gID,type)] {
	set grafdata($gID,type) uknown
    }
    if ![info exists grafdata($gID,Xtitle)] {
	set grafdata($gID,Xtitle) $graphtype
    }
    if ![info exists grafdata($gID,Xicon)]  {
	set grafdata($gID,Xicon) $graphtype
    }

    set top [xcToplevel .graph$gID $grafdata($gID,Xtitle) \
	    $grafdata($gID,Xicon) . 100 100 0]
    
    # INITIALIZATIONS
    # ---------------
    # can ... canvas
    # con ... control widget
    if ![info exists grafsize($gID,canW)] {
	set grafsize($gID,canW) $grafsize(canW)
    }
    if ![info exists grafsize($gID,canH)] {
	set grafsize($gID,canH) $grafsize(canH)
    }
    set grafsize($gID,conW) $grafsize($gID,canW)
    set grafsize($gID,conH) 70
    
    PlaceGrapher $gID $graphtype

    #tkwait visibility $top
    update

    set grafsize($gID,tplwH) [winfo height $top]
    set grafsize($gID,tplwW) [winfo width  $top]
    set grafsize($gID,Yoffset) \
	    [expr $grafsize($gID,tplwH) - $grafsize($gID,canH) \
	    - $grafsize($gID,conH)]
    set grafsize($gID,Xoffset) \
	    [expr $grafsize($gID,tplwW) - $grafsize($gID,canW)]

    bind $top <Configure> [list ResizeGrapher $gID $graphtype]
    return $top
}


proc ResizeGrapher {gID com} {
    global grafsize grafdata graf

    set h   [winfo height .graph$gID]
    set w   [winfo width .graph$gID]

    set can $graf($gID,can)
    # update only if size of "toplevel .graph" has changed
    if { $grafsize($gID,tplwH) !=  $h || $grafsize($gID,tplwW) != $w } {
	set grafsize($gID,tplwH) $h
	set grafsize($gID,tplwW) $w
	set grafsize($gID,canW) [expr $w - $grafsize($gID,Xoffset)]
	set grafsize($gID,canH) [expr $h - $grafsize($gID,conH) - \
		$grafsize($gID,Yoffset)]
	$can configure \
		-height $grafsize($gID,canH) -width $grafsize($gID,canW)
	$can delete graf

	#
	# do we have radiobuttons
	#
	if { $grafdata($gID,N_graf) > 1 } {
	    $can coords radiobuttons [expr $grafsize($gID,canW) - 10] 10
	}

	$com $gID $can
    }
}

    
proc PlaceGrapher {gID command} {
    global grafsize grafdata graf

    set f [frame .graph$gID.f -bg "#bbb" -relief raised -bd 2]

    if { ![info exists grafsize($gID,Yscroll)] } {
	set grafsize($gID,Yscroll) 0
    }

    if { $grafsize($gID,Yscroll) == 0 } {
	set can [canvas $f.can -width $grafsize($gID,canW) \
		-height $grafsize($gID,canH) -bg "#ffffff"]
    } else {
	set region [list 0 0 $grafsize($gID,canW) $grafsize($gID,Yscroll)]
	set can [canvas $f.can -width $grafsize($gID,canW) \
		-height $grafsize($gID,canH) -bg "#ffffff" \
		-scrollregion $region \
		-yscrollcommand [list $f.yscroll set]]
	set scroll [scrollbar $f.yscroll -orient vertical \
		-command [list $f.can yview]]
    }
    set com [frame .graph$gID.con -relief raised -width $grafsize($gID,conW) \
		 -height $grafsize($gID,conH)]
	
    pack $f -side top -fill both -expand 1
    pack $can -side left -expand 1 -fill both -padx 3 -pady 3

    if { [info exists scroll] } { pack $scroll -side right -fill y }
    pack $com -side top -expand 1 -fill x

    #
    # get the deafult canvas font
    #
    $can create text 0 0 -text a -tags testfont
    set grafdata($gID,deffont) [$can itemcget testfont -font]
    # default font should not be bold
    set grafdata($gID,deffont) \
	    [ModifyFont $grafdata($gID,deffont) $can \
	    -default 1 \
	    -weight  normal]
    $can delete testfont

    #
    # Canvas(item) Bindings
    #
    focus $can
    $can bind graf <1> [list GrapherSelect $gID $can %x %y]
    $can bind graf <Double-1> \
	    [list GrapherSelect $gID $can %x %y $command double $command]
    $can bind graf <2> \
	    [list GrapherSelect $gID $can %x %y $command double $command]
    $can bind graf <Button1-ButtonRelease> [list GrapherSelectRelease $gID]

    bind $can <Motion> [list GrapherSelectMotion $gID $can %x %y]

    bind .graph$gID <Left>    [list GrapherKeyBindings $gID $can left]
    bind .graph$gID <Right>   [list GrapherKeyBindings $gID $can right]
    bind .graph$gID <Up>      [list GrapherKeyBindings $gID $can up]
    bind .graph$gID <Down>    [list GrapherKeyBindings $gID $can down]
    bind .graph$gID <Delete>  [list GrapherKeyBindings $gID $can delete]
    bind .graph$gID <Destroy> [list QuitGrapher $gID]

    ##################################################
    # Command Buttons                                #
    ##################################################
    set ranges [button $com.ranges \
	    -text    "XY Ranges" \
	    -command [list ConfigGrapher $gID $command ranges]]

    set tick [button $com.tics \
	    -text "Set Tics" \
	    -command [list ConfigGrapher $gID $command tics]]

    set grid [button $com.grid \
	    -text "Set Grid" \
	    -command [list ConfigGrapher $gID $command grid]]

    set text [button $com.text \
	    -text "Text" \
	    -command [list ConfigGrapher $gID $command init_text]]

    set marg [button $com.marg \
	    -text "Margins" \
	    -command [list ConfigGrapher $gID $command margin]]
    set empty [frame $com.f -relief raised -bd 2 -highlightthickness 0]
    set dos [button $com.dos \
	    -text "DOS" \
	    -state disabled \
	    -command [list ConfigGrapher $gID $command dos]]
    if { $grafdata($gID,type) == "DOSS" } {
	$dos configure -state normal
    }
    
    set print [button $com.ps \
	    -text "PostScript" \
	    -command [list GraphPrint $gID $can]]
    set quit [button $com.quit \
	    -text "Close" \
	    -command [list QuitGrapher $gID]]
    pack $ranges $tick $grid $text $marg $dos -side left
    pack $empty -side left -fill both -expand 1
    pack $print $quit -side left
    ##################################################

    #
    # if there will be more then one graph ->
    # radiobuttons for selectiong the graph
    #
    set graf($gID,active_graf) 1; # don't remove, please; it won't work
    if { $grafdata($gID,N_graf) > 1 } {	
	set f [frame $can.f]
	$can create window [expr $grafsize($gID,canW) - 10] 10 \
		-anchor ne \
		-window $f \
		-tags "radiobuttons"
	for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {	    
	    set r($g) [radiobutton $f.r$g \
		    -variable graf($gID,active_graf) \
		    -text $g \
		    -value $g \
		    -command [list GrapherActivate $gID $can] \
		    -relief raised \
		    -bd 2]
	    pack $r($g) -side left
	}
    }

    #
    # remember the pathname of canvas
    #
    set graf($gID,can) $can

    #update
    tkwait visibility $can

    ##################
    $command $gID $can
}


proc GrapherActivate {gID can} {
    global graf
    #
    # when some new graph is selected, the selected graph should blink
    #
    set g $graf($gID,active_graf)
    BlinkingRectangle $can 3 \
	    $graf($gID,Xstart,$g) $graf($gID,Ystart,$g) \
	    $graf($gID,Xend,$g) $graf($gID,Yend,$g)
}    


proc UpdateGrapher {gID can graphcommand} {
    $can delete graf
    $graphcommand $gID $can
}


proc QuitGrapher gID {
    global graf grafdata grafsize grafselection graf_setgrid graf_setaxe \
	    graf_settic graf_setticX graf_setticY graf_setmar configGrapher \
	    quitGrapher

    #
    # because we have <Destroy> binding for toplevel, this procedure can
    # be entered many times; prevent that
    #
    if ![info exists quitGrapher($gID)] { 
	set quitGrapher($gID) 0 
    } else { 
	set quitGrapher($gID) 1
    }

    if $quitGrapher($gID) {
	return
    }

    destroy .graph$gID
    set array_to_unset {graf grafdata grafsize graf_setgrid graf_setaxe \
	    graf_setticX graf_setticY graf_setdos graf_setmar \
	    configGrapher}

    foreach array $array_to_unset {
	if [array exists $array] {
	    foreach elem [array names $array $gID,*] {
		upvar $array ar
		unset ar($elem)
	    }
	}
	xcDebug "array: $array; elem: [array names $array *]\n"
    }
}
   

proc GraphPrint {gID can} {
    global grafsize

    # delete selection
    if { $grafsize($gID,Yscroll) != 0 } {
	####################
	# canvas is scrolled
	xcPrintCanvas $can -canvas_scrollheight $grafsize($gID,Yscroll)
    } else {
	########################
	# canvas is not scrolled
	xcPrintCanvas $can
    }
}


proc gr_isInside {gID x y g} {
    global grafdata
    # is_inside_graph_box ???

    if { $x >= $grafdata($gID,Xmin,$g) && $x <= $grafdata($gID,Xmax,$g) && \
	    $y >= $grafdata($gID,Ymin,$g) && $y <= $grafdata($gID,Ymax,$g) } {
	return 1
    }
    return 0
}


proc gr_OIInterp {gID x1 y1 x2 y2 g} {
    global grafdata
    # OI - OutIn interpolation

    set k [expr ($y2 - $y1) / ($x2 - $x1)]
    set n [expr $y1 - $k * $x1]

    set Xl $grafdata($gID,Xmin,$g)
    if { $x1 <= $Xl && $x2 >= $Xl} {
	return "$Xl [expr $k * $Xl + $n]"
    }
    if { $k == 0 } {
	return "$x2 $y2"
    }
    set Yt $grafdata($gID,Ymax,$g)
    if { $y1 >= $Yt && $y2 <= $Yt } {
	return "[expr ($Yt - $n) / $k] $Yt"
    }
    set Yb $grafdata($gID,Ymin,$g)
    if { $y1 <= $Yb && $y2 >= $Yb } {
	return "[expr ($Yb - $n) / $k] $Yb"
    }

    # this can happend the fisrt time we enter gr_OIInterp
    return "$x1 $y1"
}


proc gr_IOInter {gID x1 y1 x2 y2 g} {
    global grafdata
    # IO - InOut interpolation

    set k [expr ($y2 - $y1) / ($x2 - $x1)]
    set n [expr $y1 - $k * $x1]

    set Xr $grafdata($gID,Xmax,$g)
    if { $x1 <= $Xr && $x2 >= $Xr} {
	return "$Xr [expr $k * $Xr + $n]"
    }
    if { $k == 0 } {
	return "$x2 $y2"
    }
    set Yt $grafdata($gID,Ymax,$g)
    if { $y1 <= $Yt && $y2 >= $Yt } {
	return "[expr ($Yt - $n) / $k] $Yt"
    }
    set Yb $grafdata($gID,Ymin,$g)
    if { $y1 >= $Yb && $y2 <= $Yb } {
	return "[expr ($Yb - $n) / $k] $Yb"
    }
    xcDebug "gr_IOInter ERROR: no interpolation found::"
    xcDebug "$x1 $y1 $x2 $y2; $Yt $Yb $Xr"
}


proc XYGraph {gID can} {
    global grafdata grafsize graf
    
    GetGraphWorld $gID $can

    xcDebug "XYGraph"

    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	set Xs  $graf($gID,Xstart,$g)
	set Xr  $graf($gID,Xratio,$g)
	set Xm  $grafdata($gID,Xmin,$g)
	set XM  $grafdata($gID,Xmax,$g)
	set Ys  $graf($gID,Ystart,$g)
	set Yr  $graf($gID,Yratio,$g)
	set Ym  $grafdata($gID,Ymin,$g)
	set YM  $grafdata($gID,Ymax,$g)
	for {set i 1} {$i <= $grafdata($gID,N_segment,$g)} {incr i} {
	    if [array exists line] {unset line}
	    # display curve
	    ################test-setup
	    #set nl        1
	    #set oldinside 1
	    #set inside    1
	    ################
	    set nl        0
	    set oldinside 0
	    set xold $grafdata($gID,X1,$g)
	    set yold $grafdata($gID,1,$i,$g)
	    for {set j 2} {$j <= $grafdata($gID,N_point,$g)} {incr j} {
		set x $grafdata($gID,X$j,$g)
		set y $grafdata($gID,$j,$i,$g)
		set inside [gr_isInside $gID $x $y $g]
		#xcDebug "$oldinside $inside"
		if { $oldinside && $inside } {
		    set X [expr $Xs + ($x - $Xm) * $Xr]
		    set Y [expr $Ys - ($y - $Ym) * $Yr]
		    set xold $x
		    set yold $y
		    append line($nl) "$X $Y "
		    set inside $oldinside
		    continue
		}
		if { !$oldinside && $inside } {
		    set xy [gr_OIInterp $gID $xold $yold $x $y $g]
		    set Xold [expr $Xs + ([lindex $xy 0] - $Xm) * $Xr] 
		    set Yold [expr $Ys - ([lindex $xy 1] - $Ym) * $Yr]
		    set X [expr $Xs + ($x - $Xm) * $Xr]
		    set Y [expr $Ys - ($y - $Ym) * $Yr]
		    set xold $x
		    set yold $y
		    incr nl
		    append line($nl) "$Xold $Yold $X $Y "
		    set oldinside $inside
		    continue
		}
		if { $oldinside && !$inside } {
		    set xy [gr_IOInter $gID $xold $yold $x $y $g]
		    set x [lindex $xy 0]
		    set y [lindex $xy 1]
		    set X [expr $Xs + ($x - $Xm) * $Xr]
		    set Y [expr $Ys - ($y - $Ym) * $Yr]
		    set xold $x
		    set yold $y
		    append line($nl) "$X $Y "
		    set oldinside $inside
		    continue
		}
		if { !$oldinside && !$inside } {
		    set xold $x
		    set yold $y
		    continue
		}
	    }
	    set curve_width $grafsize(curve_width)
	    set curve_color $grafsize(curve_color)
	    if [info exists grafsize($gID,curve_width_$i,$g)] {
		set curve_width $grafsize($gID,curve_width_$i,$g)
	    }
	    if [info exists grafsize($gID,curve_color_$i,$g)] {
		set curve_color $grafsize($gID,curve_color_$i,$g)
	    }	
	    for {set il 1} {$il <= $nl} {incr il} {		
		eval {$can create line} $line($il) {-width $curve_width \
			-tags "graf curve c$g,$i" -smooth 1 \
			-fill $curve_color -joinstyle round}
	    }
	}
    }


    #####################################################################
    # it may easily happen, that some curve goes outside the graph world,
    # get rid off that --> make white rectangles (JUST IN CASE)
    GrapherDummyRect $gID $can
    
    GrapherArbitraryText $gID $can    
    XYAxis $gID $can
    GraphTickLines $gID $can

    # we must raise the radiobuttons
    $can raise radiobuttons dummy_rectangles
}


proc BARGraph {gID can} {    
    global grafdata grafsize graf prop

    GetGraphWorld $gID $can

    set grafdata($gID,type) BAR

    # display BARS
    for {set i 1} {$i <= $grafdata($gID,N_point,1)} {incr i} {
	# bar is (1/N_point) - 2*Xoffset*(1/N_point) tick
	set x1 [expr $grafdata($gID,X$i,1) - \
		(1.0 - 2.0 * $grafdata($gID,Xoffset,1)) * \
		(0.5 / $grafdata($gID,N_point,1))]
	set x2  [expr $grafdata($gID,X$i,1) + \
		(1.0 - 2.0 * $grafdata($gID,Xoffset,1)) * \
		(0.5 / $grafdata($gID,N_point,1))]

	set Y1 [expr $graf($gID,Ystart,1) - \
		($grafdata($gID,$i,1,1) - $grafdata($gID,Ymin,1)) * \
		$graf($gID,Yratio,1)] 
	set Y2 [expr $graf($gID,Ystart,1) - \
		($grafdata($gID,$i,2,1) - $grafdata($gID,Ymin,1)) * \
		$graf($gID,Yratio,1)]

	if { $prop(type_of_run) == "UHF" } {
	    set x4 $x2
	    set x2 $grafdata($gID,X$i,1)	    
	    set X4 [expr $graf($gID,Xstart,1) + \
		    ($x4 - $grafdata($gID,Xmin,1)) * $graf($gID,Xratio,1)]
	    set Y3 [expr $graf($gID,Ystart,1) - \
		    ($grafdata($gID,$i,1,2) - $grafdata($gID,Ymin,1)) * \
		    $graf($gID,Yratio,1)] 
	    set Y4 [expr $graf($gID,Ystart,1) - \
		    ($grafdata($gID,$i,2,2) - $grafdata($gID,Ymin,1)) * \
		    $graf($gID,Yratio,1)]
	}

	set X1 [expr $graf($gID,Xstart,1) + \
		($x1 - $grafdata($gID,Xmin,1)) * $graf($gID,Xratio,1)]
	set X2 [expr $graf($gID,Xstart,1) + \
		($x2 - $grafdata($gID,Xmin,1)) * $graf($gID,Xratio,1)]
	set X3 $X2; # used just for UHF


	#
	# convert "none" to {}
	#
	if { $grafsize($gID,bar${i}_stipple,1) == "none" } {
	    set grafsize($gID,bar${i}_stipple,1) {}
	}
	if { $prop(type_of_run) == "UHF" } {
	    if { $grafsize($gID,bar${i}_stipple,2) == "none" } {
		set grafsize($gID,bar${i}_stipple,2) {}
	    }
	}

	#
	# shadow bar
	#
	set off 0
	if [info exists grafdata($gID,barshadow,1)] {
	    set off 1
	    if $grafdata($gID,barshadow,1) {
		$can create rectangle \
			[expr $X1 + 2] [expr $Y1 + 2] \
			[expr $X2 + 2] [expr $Y2 + 2] \
			-fill #ffffff \
			-outline #ffffff \
			-tags "graf"		
		$can create rectangle \
			[expr $X1 + 2] [expr $Y1 + 2] \
			[expr $X2 + 2] [expr $Y2 + 2] \
			-fill    $grafsize($gID,bar${i}_shadow,1) \
			-outline $grafsize($gID,bar${i}_shadow,1) \
			-stipple $grafsize($gID,bar${i}_stipple,1) \
			-width   $grafsize($gID,bar${i}_width,1) \
			-stipple $grafsize($gID,bar${i}_stipple,1) \
			-tags    "graf sbar sb1,$i"
		if { $prop(type_of_run) == "UHF" } {
		    $can create rectangle \
			    [expr $X3 + 2] [expr $Y3 + 2] \
			    [expr $X4 + 2] [expr $Y4 + 2] \
			    -fill #ffffff \
			    -outline #ffffff \
			    -tags "graf"	
		    $can create rectangle \
			    [expr $X3 + 2] [expr $Y3 + 2] \
			    [expr $X4 + 2] [expr $Y4 + 2] \
			    -fill    $grafsize($gID,bar${i}_shadow,2) \
			    -outline $grafsize($gID,bar${i}_shadow,2) \
			    -stipple $grafsize($gID,bar${i}_stipple,2) \
			    -width   $grafsize($gID,bar${i}_width,2) \
			    -stipple $grafsize($gID,bar${i}_stipple,2) \
			    -tags    "graf sbar sb2,$i"
		}	    
	    }
	}

	$can create rectangle \
		[expr $X1 - $off] $Y1 [expr $X2 - $off] $Y2 \
		-fill #ffffff \
		-outline #ffffff \
		-tags "graf"	
	$can create rectangle \
		[expr $X1 - $off] $Y1 [expr $X2 - $off] $Y2 \
		-fill    $grafsize($gID,bar${i}_fill,1) \
		-outline $grafsize($gID,bar${i}_outline,1) \
		-stipple $grafsize($gID,bar${i}_stipple,1) \
		-width   $grafsize($gID,bar${i}_width,1) \
		-stipple $grafsize($gID,bar${i}_stipple,1) \
		-tags    "graf bar b1,$i"

	if { $prop(type_of_run) == "UHF" } {
	    $can create rectangle \
		    [expr $X3 - $off] $Y3 [expr $X4 - $off] $Y4 \
		    -fill #ffffff \
		    -outline #ffffff \
		    -tags "graf"	
	    $can create rectangle \
		    [expr $X3 - $off] $Y3 [expr $X4 - $off] $Y4 \
		    -fill    $grafsize($gID,bar${i}_fill,2) \
		    -outline $grafsize($gID,bar${i}_outline,2) \
		    -stipple $grafsize($gID,bar${i}_stipple,2) \
		    -width   $grafsize($gID,bar${i}_width,2) \
		    -stipple $grafsize($gID,bar${i}_stipple,2) \
		    -tags    "graf bar b2,$i"
	}
    }

    GrapherDummyRect $gID $can
    GrapherArbitraryText $gID $can
    XYAxis $gID $can
    GraphTickLines $gID $can
}

proc GrapherArbitraryText {gID can} {
    global grafdata grafsize

    #
    # display arbitrary texts
    #
    for {set i 1} {$i <= $grafdata($gID,N_text)} {incr i} {
	# display text
	if [info exists grafdata($gID,text${i}_X)] {	    
	    if ![info exists grafdata($gID,text${i}_font)] {
		set grafdata($gID,text${i}_font) $grafdata($gID,deffont)
	    }	    
	    xcDebug "font:: $grafdata($gID,text${i}_font)"
	    $can create text \
		    [expr [XnMRel2Graph $gID $grafdata($gID,text${i}_X)] + \
		    $grafsize($gID,margin_Y1)] \
		    [expr [YnMRel2Graph $gID $grafdata($gID,text${i}_Y)] + \
		    $grafsize($gID,margin_X2)] \
		    -anchor nw \
		    -text $grafdata($gID,text$i) \
		    -tags "graf text t$i" \
		    -font $grafdata($gID,text${i}_font)
	    # display bbox ???
	    if { [info exists grafdata($gID,text${i}_bbox) ] } {
		if $grafdata($gID,text${i}_bbox) {
		    set bbox [$can bbox t$i]
		    eval {$can create rectangle} $bbox \
			    {-outline #000 \
			    -tags "graf text_bbox t_bb$i"}
		} else {
		    $can delete t_bb$i
		}
	    }
	}
    }
}


proc GrapherDummyRect {gID can} {
    global graf grafsize

    ######################################################################
    # this will be used for dummy_rectangles -> look below !!!
    set cw  $grafsize($gID,canW)
    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	set ch $grafsize($gID,Yscroll)
    } else {
	set ch $grafsize($gID,canH)
    }

    # left rectangle
    $can create rectangle 0 0 $graf($gID,X1) $ch \
	    -fill      "#ffffff" \
	    -outline   "#ffffff" \
	    -tags      "graf dummy_rectangles"
    # top rectangle
    $can create rectangle 0 0 $cw $graf($gID,Y1) \
	    -fill      "#ffffff" \
	    -outline   "#ffffff" \
	    -tags      "graf dummy_rectangles"
    # bottom rectangle
    $can create rectangle 0 $graf($gID,Y2) $cw $ch \
	    -fill      "#ffffff" \
	    -outline   "#ffffff" \
	    -tags      "graf dummy_rectangles"
    # right rectangle
    $can create rectangle $graf($gID,X2) 0 $cw $ch \
	    -fill      "#ffffff" \
	    -outline   "#ffffff" \
	    -tags      "graf dummy_rectangles"

}
		

proc GraphTickLines {gID can} {
    global grafdata grafsize graf graf_setgrid graf_setaxe system

    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	# display Xline
	if ![info exists grafdata($gID,XYline_textX,$g)] {
	    set grafdata($gID,XYline_textX,$g)   $grafdata(XYline_textX)
	    set grafdata($gID,XYline_textY,$g)   $grafdata(XYline_textY)
	}
	    
	set xt [expr $graf($gID,Xstart,$g) + \
		[XnMRel2Graph $gID [expr $graf($gID,RorigX,$g) + \
		(1.0 - $grafdata($gID,XYline_textX,$g)) * \
		$graf($gID,RsizeX,$g)]]]
	
	set yt [expr $graf($gID,Ystart,$g) - \
		[YnMRel2Graph $gID [expr $graf($gID,RorigY,$g) + \
		(1.0 - $grafdata($gID,XYline_textY,$g)) * \
		$graf($gID,RsizeX,$g)]]]

	if ![info exists grafdata($gID,legendfont)] {
	    set grafdata($gID,legendfont) $grafdata($gID,deffont)
	}
	
	for {set i 1} {$i <= $grafdata($gID,N_Xline,$g)} {incr i} { 
	    set Xline [expr $graf($gID,Xstart,$g) + \
		    ($grafdata($gID,Xline$i,$g) - \
		    $grafdata($gID,Xmin,$g)) * $graf($gID,Xratio,$g)]
	    puts stdout "Xline: $Xline"

	    if { $Xline > $graf($gID,Xstart,$g) && \
		    $Xline < $graf($gID,Xend,$g) } {
		$can create line \
			$Xline $graf($gID,Ystart,$g) \
			$Xline $graf($gID,Yend,$g) \
			-width   $grafsize($gID,Xline${i}_width,$g) \
			-fill    $grafsize($gID,Xline${i}_fill,$g) \
			-stipple $grafsize($gID,Xline${i}_stipple,$g) \
			-tags    "graf Xline Xl$g,$i" 
	    
		# is there any label for Xline$i
		if [info exists grafdata($gID,Xline${i}_text,$g)] {
		    $can create line \
			    [expr $xt - 20] $yt $xt $yt \
			    -width   $grafsize($gID,Xline${i}_width,$g) \
			    -fill    $grafsize($gID,Xline${i}_fill,$g) \
			    -stipple $grafsize($gID,Xline${i}_textstipple,$g) \
			    -tags    "graf legend l$g Xline_textline Xl_tl$g,$i"
		    $can create text \
			    [expr $xt - 25] $yt \
			    -text $grafdata($gID,Xline${i}_text,$g) \
			    -anchor e \
			    -fill $grafsize($gID,Xline${i}_fill,$g) \
			    -tags "graf legend l$g Xline_text Xl_t$g,$i" \
			    -font $grafdata($gID,legendfont)
		    
		    # t.k - this is strange
		    set yt [expr $yt + 20]; #MAKE THIS BETTER IN THE FUTURE
		}
	    }
	}

	# display Yline
	for {set i 1} {$i <= $grafdata($gID,N_Yline,$g)} {incr i} { 
	    set Yline [expr $graf($gID,Ystart,$g) - \
		    ($grafdata($gID,Yline$i,$g) - \
		    $grafdata($gID,Ymin,$g)) * $graf($gID,Yratio,$g)]
	    if { $Yline < $graf($gID,Ystart,$g) && \
		    $Yline > $graf($gID,Yend,$g) } {
		$can create line \
			$graf($gID,Xstart,$g) $Yline \
			$graf($gID,Xend,$g) $Yline \
			-width   $grafsize($gID,Yline${i}_width,$g) \
			-fill    $grafsize($gID,Yline${i}_fill,$g) \
			-stipple $grafsize($gID,Yline${i}_stipple,$g) \
			-tags    "graf Yline Yl$g,$i"
		# is there any label for Yline$i
		if [info exists grafdata($gID,Yline${i}_text,$g)] {
		    $can create line \
			    [expr $xt - 20] $yt $xt $yt \
			    -width   $grafsize($gID,Yline${i}_width,$g) \
			    -fill    $grafsize($gID,Yline${i}_fill,$g) \
			    -stipple $grafsize($gID,Yline${i}_stipple,$g) \
			    -tags    "graf legend l$g Yline_textline Yl_tl$g,$i"
		    $can create text [expr $xt - 25] $yt \
			    -text $grafdata($gID,Yline${i}_text,$g) \
			    -anchor e \
			    -fill $grafsize($gID,Yline${i}_fill,$g) \
			    -tags "graf legend l$g Yline_text Yl_t$g,$i" \
			    -font $grafdata($gID,legendfont)
		    set yt [expr $yt + 20]; #MAKE THIS BETTER IN THE FUTURE
		}
	    }
	}
	
	# display X ticks
	if ![info exists grafdata($gID,Xtick_draw,$g)] {
	    set grafdata($gID,Xtick_draw,$g)      $grafdata(Xtick_draw)
	}
	if ![info exists graf_setaxe($gID,mirrorXaxe,$g)] {
	    set graf_setaxe($gID,mirrorXaxe,$g)   $graf_setaxe(mirrorXaxe)
	}
	if ![info exists grafdata($gID,Xtick_text_draw,$g)] {
	    set grafdata($gID,Xtick_text_draw,$g) $grafdata(Xtick_text_draw)
	}
	if ![info exists grafdata($gID,Xtickfont)] {
	    set grafdata($gID,Xtickfont) $grafdata($gID,deffont)
	}

	for {set i 1} {$i <= $grafdata($gID,N_MXtick,$g)} {incr i} {
	    set Xth  [XGraphValue $gID $grafdata($gID,Xtick$i,$g) $g]
	    set Yth  [expr $graf($gID,Ystart,$g) + $grafsize(tt_offset)]
	    set Yth1 [expr $graf($gID,Ystart,$g) - $grafsize(Mtick_size)]
	    set Yth2 $graf($gID,Ystart,$g) 
	    # mirror axis:
	    set mYth  [expr $graf($gID,Yend,$g) - $grafsize(tt_offset)]
	    set mYth1 [expr $graf($gID,Yend,$g) + $grafsize(Mtick_size)]
	    set mYth2 $graf($gID,Yend,$g) 
	    
	    if { $Xth >= [XGraphValue $gID $grafdata($gID,Xmin,$g) $g] && \
		    $Xth <= [XGraphValue $gID $grafdata($gID,Xmax,$g) $g] } {
		if $grafdata($gID,Xtick_draw,$g) {
		    $can create line \
			    $Xth $Yth1 $Xth $Yth2 \
			    -width 2 \
			    -fill "#000" \
			    -tags "graf XMtick XMt$g"
		    # display X mirror tics ???
		    if $graf_setaxe($gID,mirrorXaxe,$g) {
			$can create line \
				$Xth $mYth1 $Xth $mYth2 \
				-width 2 \
				-fill "#000" \
				-tags "graf XMtick mXaxe mXMt$g mXa$g"
			    }
		}
		if $grafdata($gID,Xtick_text_draw,$g) {
		    $can create text \
			    $Xth $Yth \
			    -text $grafdata($gID,Xtick${i}_text,$g) \
			    -anchor n \
			    -tags "graf tickXValue tXV$g" \
			    -font $grafdata($gID,Xtickfont)		    
		}
		#
		# display X major grid if graf_setgrid is set
		#
		if ![info exists graf_setgrid($gID,Xmajor,$g)] {
		    set graf_setgrid($gID,Xmajor,$g) $graf_setgrid(Xmajor) 
		    set graf_setgrid($gID,Xminor,$g) $graf_setgrid(Xminor) 
		    set graf_setgrid($gID,Ymajor,$g) $graf_setgrid(Ymajor) 
		    set graf_setgrid($gID,Yminor,$g) $graf_setgrid(Yminor)
		    set graf_setgrid($gID,Xmajor_exist,$g) \
			    $graf_setgrid(Xmajor_exist) 
		    set graf_setgrid($gID,Ymajor_exist,$g) \
			    $graf_setgrid(Ymajor_exist) 
		    set graf_setgrid($gID,Xminor_exist,$g) \
			    $graf_setgrid(Xminor_exist) 
		    set graf_setgrid($gID,Yminor_exist,$g) \
			    $graf_setgrid(Yminor_exist) 
		    set graf_setaxe($gID,mirrorXaxe,$g) \
			    $graf_setaxe(mirrorXaxe)
		    set graf_setaxe($gID,mirrorYaxe,$g) \
			    $graf_setaxe(mirrorYaxe)
		}
		if $graf_setgrid($gID,Xmajor,$g) {
		    set graf_setgrid($gID,Xmajor_exist,$g) 1
		    $can create line \
			    $Xth $graf($gID,Yend,$g) \
			    $Xth $graf($gID,Ystart,$g) \
			    -width 0.5 \
			    -fill "#999" \
			    -tags "graf XMgrid XMg$g" \
			    -stipple @$system(BMPDIR)/dot1x3V.bmp
		    # gridlines must be the lowest ones
		    $can lower XMgrid graf
		} elseif { $graf_setgrid($gID,Xmajor,$g) == 0 && \
			$graf_setgrid($gID,Xmajor_exist,$g) == 1 } {
		    # delete gridlines
		    $can delete XMg$g
		    set graf_setgrid($gID,Xmajor_exist,$g) 0
		}
	    }
	    
	    if { $i < $grafdata($gID,N_MXtick,$g) } {
		set ii [expr $i + 1]
		set dx [expr ($grafdata($gID,Xtick$ii,$g) - \
			$grafdata($gID,Xtick$i,$g)) \
			* $graf($gID,Xratio,$g) / \
			($grafdata($gID,N_mXtick,$g) + 1)]
		for {set j 1} {$j <= $grafdata($gID,N_mXtick,$g)} {incr j} {
		    set xth [expr $Xth + $j * $dx]
		    set yth1 [expr $graf($gID,Ystart,$g) - \
			    $grafsize(mtick_size)]
		    set myth1 [expr $graf($gID,Yend,$g) + \
			    $grafsize(mtick_size)]
		    if { $xth >= \
			    [XGraphValue $gID $grafdata($gID,Xmin,$g) $g] \
			    && $xth <= \
			    [XGraphValue $gID $grafdata($gID,Xmax,$g) $g] } {
			if $grafdata($gID,Xtick_draw,$g) {
			    $can create line \
				    $xth $yth1 $xth $Yth2 \
				    -width 1 \
				    -fill "#000" \
				    -tags "graf Xmtick Xmt$g"
			    # display X mirror tics ???
			    if $graf_setaxe($gID,mirrorXaxe,$g) {
				$can create line \
					$xth $myth1 $xth $mYth2 \
					-width 1 \
					-fill "#000" \
					-tags "graf Xmtick mXaxe mXmt$g mXa$g"
			    }
			}
			#
			# display X minor grid if graf_setgrid is set
			#
			if $graf_setgrid($gID,Xminor,$g) {
			    set graf_setgrid($gID,Xminor_exist,$g) 1
			    $can create line \
				    $xth $graf($gID,Yend,$g) \
				    $xth $graf($gID,Ystart,$g) \
				    -width 0.5 \
				    -fill "#999" \
				    -tags "graf Xmgrid Xmg$g" \
				    -stipple @$system(BMPDIR)/dot1x3V.bmp
			    # gridlines must be the lowest ones
			    $can lower Xmg$g graf
			} elseif { $graf_setgrid($gID,Xminor,$g) == 0 && \
				$graf_setgrid($gID,Xminor_exist,$g) == 1 } {
			    # delete gridlines
			    $can delete Xmg$g
			    set graf_setgrid($gID,Xminor_exist,$g) 0
			}
		    }
		}
	    }
	}
	# display Y ticks
	if ![info exists grafdata($gID,Ytick_draw,$g)] {
	    set grafdata($gID,Ytick_draw,$g)      $grafdata(Ytick_draw)
	}
	if ![info exists graf_setaxe($gID,mirrorYaxe,$g)] {
	    set graf_setaxe($gID,mirrorYaxe,$g)   $graf_setaxe(mirrorYaxe)
	}
	if ![info exists grafdata($gID,Ytick_text_draw,$g)] {
	    set grafdata($gID,Ytick_text_draw,$g) $grafdata(Ytick_text_draw)
	}
	if ![info exists grafdata($gID,Ytickfont)] {
	    set grafdata($gID,Ytickfont) $grafdata($gID,deffont)
	}

	for {set i 1} {$i <= $grafdata($gID,N_MYtick,$g)} {incr i} {
	    set Xth  [expr $graf($gID,Xstart,$g) - $grafsize(tt_offset)]
	    set Yth  [YGraphValue $gID $grafdata($gID,Ytick$i,$g) $g] 
	    set Xth1 $graf($gID,Xstart,$g)
	    set Xth2 [expr $graf($gID,Xstart,$g) + $grafsize(Mtick_size)]
	    # mirror Axis
	    set mXth  [expr $graf($gID,Xend,$g) + $grafsize(tt_offset)]
	    set mXth1 $graf($gID,Xend,$g)
	    set mXth2 [expr $graf($gID,Xend,$g) - $grafsize(Mtick_size)]
	
	    #
	    # WARNING: point (0,0) -> upper left corner of canvas
	    #             >>> take care of that <<<
	    #
	    if { $Yth <= [YGraphValue $gID $grafdata($gID,Ymin,$g) $g] && \
		    $Yth >= [YGraphValue $gID $grafdata($gID,Ymax,$g) $g] } {
		if $grafdata($gID,Ytick_draw,$g) {
		    $can create line \
			    $Xth1 $Yth $Xth2 $Yth \
			    -width 2 \
			    -fill "#000" \
			    -tags "graf YMtick YMt$g"	    
		    # display mirror tics ???
		    if $graf_setaxe($gID,mirrorYaxe,$g) {
			$can create line \
				$mXth1 $Yth $mXth2 $Yth \
				-width 2 \
				-fill "#000" \
				-tags "graf YMtick mYaxe mYMt$g mXa$g"
		    }
		}
		if $grafdata($gID,Ytick_text_draw,$g) {
		    $can create text $Xth $Yth \
			    -text $grafdata($gID,Ytick${i}_text,$g) \
			    -anchor e \
			    -tags "graf tickYValue tYV$g" \
			    -font $grafdata($gID,Ytickfont)
		    

		}
		#
		# display Y major grid if graf_setgrid is set
		#
		if $graf_setgrid($gID,Ymajor,$g) {
		    set graf_setgrid($gID,Ymajor_exist,$g) 1
		    $can create line \
			    $graf($gID,Xstart,$g) $Yth \
			    $graf($gID,Xend,$g) $Yth \
			    -width 0.5 \
			    -fill "#999" \
			    -tags "graf YMgrid YMg$g" \
			    -stipple @$system(BMPDIR)/dot1x3H.bmp
		    # gridlines must be the lowest ones
		    $can lower YMgrid graf
		} elseif { $graf_setgrid($gID,Ymajor,$g) == 0 && \
			$graf_setgrid($gID,Ymajor_exist,$g) == 1 } {
		    # delete gridlines
		    $can delete YMg$g
		    set graf_setgrid($gID,Ymajor_exist,$g) 0
		}		
	    }
	    if { $i < $grafdata($gID,N_MYtick,$g) } {
		set ii [expr $i + 1]
		set dy [expr ($grafdata($gID,Ytick$ii,$g) - \
			$grafdata($gID,Ytick$i,$g)) \
			* $graf($gID,Yratio,$g) / \
			($grafdata($gID,N_mYtick,$g) + 1)]	    
		for {set j 1} {$j <= $grafdata($gID,N_mYtick,$g)} {incr j} {
		    set yth [expr $Yth - $j * $dy]
		    set xth2 [expr $graf($gID,Xstart,$g) + \
			    $grafsize(mtick_size)]
		    # mirror tics
		    set mxth2 [expr $graf($gID,Xend,$g) - \
			    $grafsize(mtick_size)]
		    if { $yth <= [YGraphValue \
			    $gID $grafdata($gID,Ymin,$g) $g] \
			    && \
			    $yth >= [YGraphValue \
			    $gID $grafdata($gID,Ymax,$g) $g] } {
			if $grafdata($gID,Ytick_draw,$g) {
			$can create line \
				$Xth1 $yth $xth2 $yth \
				-width 1 \
				-fill "#000" \
				-tags "graf Ymtick Ymt$g"
			    if $graf_setaxe($gID,mirrorYaxe,$g) {
				$can create line \
					$mXth1 $yth $mxth2 $yth \
					-width 1 \
					-fill "#000" \
					-tags "graf Ymtick mYaxe mYmt$g mYa$g"
			    }
			}
			puts stdout "mticks: $Xth1 $yth $xth2 $yth"
			#
			# display Y minor grid if graf_setgrid is set
			#
			if $graf_setgrid($gID,Yminor,$g) {
			    set graf_setgrid($gID,Yminor_exist,$g) 1
			    $can create line \
				    $graf($gID,Xstart,$g) $yth \
				    $graf($gID,Xend,$g) $yth \
				    -width 0.5 \
				    -fill "#999" \
				    -tags "graf Ymgrid Ymg$g" \
				    -stipple @$system(BMPDIR)/dot1x3H.bmp
			    # gridlines must be the lowest ones
			    $can lower Ymgrid graf
			} elseif { $graf_setgrid($gID,Yminor,$g) == 0 && \
				$graf_setgrid($gID,Yminor_exist,$g) == 1 } {
			    # delete gridlines
			    $can delete Ymg$g
			    set graf_setgrid($gID,Yminor_exist,$g) 0
			}
		    }
		}
	    }
	}
    }
}


proc XGraphValue {gID value g} {
    global grafdata graf
    return [expr $graf($gID,Xstart,$g) + \
		($value - $grafdata($gID,Xmin,$g)) * $graf($gID,Xratio,$g)]
}


proc YGraphValue {gID value g} {
    global grafdata graf
    return [expr $graf($gID,Ystart,$g) - \
		($value - $grafdata($gID,Ymin,$g)) * $graf($gID,Yratio,$g)]
}


proc XDataValue {gID value g} {
    global grafdata graf
    return [expr $grafdata($gID,Xmin,$g) + \
	    double($value - $graf($gID,Xstart,$g)) / $graf($gID,Xratio,$g)]
}


proc YDataValue {gID value g} {
    global grafdata graf
    return [expr $grafdata($gID,Ymin,$g) - \
	    double($value - $graf($gID,Ystart,$g)) / $graf($gID,Yratio,$g)]
}


proc XPix2Coor {gID value g} {
    global graf

    return [expr double($value) / $graf($gID,Xratio,$g)]
}


proc YPix2Coor {gID value g} {
    global graf

    return [expr -1.0 * double($value) / $graf($gID,Yratio,$g)]
}


proc XRel2Graph {gID value} {
    global grafsize
    
    return [expr $value * double($grafsize($gID,canW))]
}
    

proc YRel2Graph {gID value} {
    global grafsize

    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	return [expr $value * double($grafsize($gID,Yscroll))]
    } else {
	return [expr $value * double($grafsize($gID,canH))]
    }
}


proc XGraph2Rel {gID value} {
    global grafsize
 
    return [expr double($value) / double($grafsize($gID,canW))]
}
    

proc YGraph2Rel {gID value} {
    global grafsize

    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	return [expr double($value) / double($grafsize($gID,Yscroll))]
    } else {
	return [expr double($value) / double($grafsize($gID,canH))]
    }
}


proc XnMRel2Graph {gID value} {
    global grafsize

    return [expr $value * double($grafsize($gID,canW) - \
	    $grafsize($gID,margin_Y1) - $grafsize($gID,margin_Y2))]
}


proc YnMRel2Graph {gID value} {
    global grafsize

    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	return [expr $value * double($grafsize($gID,Yscroll) - \
		$grafsize($gID,margin_X1) - $grafsize($gID,margin_X2))]
    } else {
	return [expr $value * double($grafsize($gID,canH) - \
		$grafsize($gID,margin_X1) - $grafsize($gID,margin_X2))]
    }
}

proc XGraph2nMRel {gID value} {
    global grafsize

    return [expr double($value) / double($grafsize($gID,canW) - \
	    $grafsize($gID,margin_Y1) - $grafsize($gID,margin_Y2))]
}


proc YGraph2nMRel {gID value} {
    global grafsize

    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	return [expr double($value) / double($grafsize($gID,Yscroll) - \
		$grafsize($gID,margin_X1) - $grafsize($gID,margin_X2))]
    } else {
	return [expr double($value) / double($grafsize($gID,canH) - \
		$grafsize($gID,margin_X1) - $grafsize($gID,margin_X2))]
    }
}


proc XYAxis {gID can} {
    global grafdata grafsize graf

    set cw  $grafsize($gID,canW)
    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	set ch $grafsize($gID,Yscroll)
    } else {
	set ch $grafsize($gID,canH)
    }

    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	# Y1-axe
	$can create line $graf($gID,Xstart,$g) $graf($gID,Yend,$g) \
		$graf($gID,Xstart,$g) $graf($gID,Ystart,$g) \
		-width $grafsize(X1axe_width) \
		-tags "graf xyaxis Y1axe Y1a$g"
	# X1-axe
	$can create line $graf($gID,Xstart,$g) $graf($gID,Ystart,$g) \
		$graf($gID,Xend,$g) $graf($gID,Ystart,$g) \
		-width $grafsize(X1axe_width) \
		-tags "graf xyaxis X1axe X1a$g"
	# X2-axe
	$can create line $graf($gID,Xend,$g) $graf($gID,Ystart,$g) \
		$graf($gID,Xend,$g) $graf($gID,Yend,$g) \
		-width $grafsize(X2axe_width) \
		-tags "graf xyaxis X2axe X2a$g"
	# Y2-axe
	$can create line $graf($gID,Xend,$g) $graf($gID,Yend,$g) \
		$graf($gID,Xstart,$g) $graf($gID,Yend,$g) \
		-width $grafsize(X2axe_width) \
		-tags "graf xyaxis Y2axe Y2a$g"
	
	#
	# Y-title
	#
	if ![info exists grafdata($gID,Ytitlefont)] {
	    set grafdata($gID,Ytitlefont) $grafdata($gID,deffont)
	}
	if [info exists grafdata($gID,Y_title,$g)] {
	    if ![info exists grafdata($gID,Y_title_X,$g)] {
		if ![info exists grafsize($gID,Yaxe_title_offset)] {
		    set grafdata($gID,Y_title_X,$g) \
			    -$grafsize(axe_title_offset)
		} else {
		    set grafdata($gID,Y_title_X,$g) \
			    -$grafsize($gID,Yaxe_title_offset)
		}
		set grafdata($gID,Y_title_Y,$g) 0.5
	    }
	    $can create text \
		    [expr $graf($gID,Xstart,$g) + \
		    $grafdata($gID,Y_title_X,$g)] \
		    [expr $grafsize($gID,margin_X2) + \
		    [YnMRel2Graph $gID [expr $graf($gID,RorigY,$g) + \
		    $graf($gID,RsizeY,$g) * $grafdata($gID,Y_title_Y,$g)]]] \
		    -anchor e \
		    -text $grafdata($gID,Y_title,$g) \
		    -tags "graf xyaxis Ytitle Yti$g" \
		    -font $grafdata($gID,Ytitlefont)
	}
	#
	# X-title
	#
	if ![info exists grafdata($gID,Xtitlefont)] {
	    set grafdata($gID,Xtitlefont) $grafdata($gID,deffont)
	}
	if [info exists grafdata($gID,X_title,$g)] {
	    if ![info exists grafdata($gID,X_title_X,$g)] {
		set grafdata($gID,X_title_X,$g) 0.5		
		if ![info exists grafsize($gID,Xaxe_title_offset)] {
		    set grafdata($gID,X_title_Y,$g) \
			    $grafsize(axe_title_offset)
		} else {
		    set grafdata($gID,X_title_Y,$g) \
			    $grafsize($gID,Xaxe_title_offset)
		}
	    }
	    $can create text \
		    [expr $grafsize($gID,margin_Y1) + \
		    [XnMRel2Graph $gID [expr $graf($gID,RorigX,$g) + \
		    $graf($gID,RsizeX,$g) * $grafdata($gID,X_title_X,$g)]]] \
		    [expr $graf($gID,Ystart,$g) + \
		    $grafdata($gID,X_title_Y,$g)] \
		    -anchor n \
		    -text $grafdata($gID,X_title,$g) \
		    -tags "graf xyaxis Xtitle Xti$g" \
		    -font $grafdata($gID,Xtitlefont)
	}
    }
}


proc GetGraphWorld {gID can} {
    global grafdata grafsize graf
    # variable "type" is obsolete now; its is here just for compatibility

    #########################################################################
    # axis square
    set cw $grafsize($gID,canW)
    if { $grafsize($gID,canH) < $grafsize($gID,Yscroll) } {
	set ch $grafsize($gID,Yscroll)
    } else {
	set ch $grafsize($gID,canH)
    }

    if ![info exists grafsize($gID,margin_X1)] {
	set grafsize($gID,margin_X1) $grafsize(margin_X1)
	set grafsize($gID,margin_X2) $grafsize(margin_X2)
	set grafsize($gID,margin_Y1) $grafsize(margin_Y1)
	set grafsize($gID,margin_Y2) $grafsize(margin_Y2)
    }

    set graf($gID,X1) $grafsize($gID,margin_Y1)
    set graf($gID,Y1) $grafsize($gID,margin_X2)
    set graf($gID,X2) [expr $cw - $grafsize($gID,margin_Y2)]
    set graf($gID,Y2) [expr $ch - $grafsize($gID,margin_X1)]
    set graf($gID,dX) [expr abs($graf($gID,X2) - $graf($gID,X1))]
    set graf($gID,dY) [expr abs($graf($gID,Y2) - $graf($gID,Y1))]

    if ![info exists grafdata($gID,N_graf)] {
	set grafdata($gID,N_graf) 1
    }
    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	# Xoffset is 0.1 --> 1.1 
	set graf($gID,Xrange,$g) \
		[expr $grafdata($gID,Xmax,$g) - $grafdata($gID,Xmin,$g)]
	set graf($gID,Yrange,$g) \
		[expr $grafdata($gID,Ymax,$g) - $grafdata($gID,Ymin,$g)]

	set graf($gID,Xratio,$g) \
		[expr $graf($gID,dX) \
		* $graf($gID,RsizeX,$g) / $graf($gID,Xrange,$g)]
	set graf($gID,Xstart,$g) \
		[expr $graf($gID,X1) + \
		[XnMRel2Graph $gID $graf($gID,RorigX,$g)]]

	set graf($gID,Yratio,$g) \
		[expr $graf($gID,dY) \
		* $graf($gID,RsizeY,$g) / $graf($gID,Yrange,$g)]
	set graf($gID,Ystart,$g) \
		[expr $graf($gID,Y2) - \
		[YnMRel2Graph $gID $graf($gID,RorigY,$g)]]
	
	set graf($gID,Xend,$g) \
		[expr $graf($gID,X1) + \
		[XnMRel2Graph $gID [expr $graf($gID,RorigX,$g) + \
		$graf($gID,RsizeX,$g)]]]
	set graf($gID,Yend,$g) \
		[expr $graf($gID,Y2) - \
		[YnMRel2Graph $gID [expr $graf($gID,RorigY,$g) + \
		$graf($gID,RsizeY,$g)]]]

    }
}



proc TickFormat {x {dx 1}} {
    # return formated tick's number
    if [IsEqual 1e-7 $x 0.0] {
	return [format "%6.3f" $x]
    } elseif { abs($x) < 0.01 || abs($x) >= 100 } {
	return [format "%6.3E" $x]
    } else {
	return [format "%6.3f" $x]
    }
}


proc GraphInit {} {
    global grafsize grafdata properties graf_setgrid \
	    graf_setaxe grafselection

    xcDebug "GraphInit"
    # 
    # some initializations for Grapher
    # 
    set grafsize(canW)              600
    set grafsize(canH)              450
    set grafsize(margin_X1)         60
    set grafsize(margin_X2)         60
    set grafsize(margin_Y1)         90
    set grafsize(margin_Y2)         60   
    set grafsize(X1axe_width)       2
    set grafsize(X2axe_width)       2
    set grafsize(axe_title_offset)  30
    set grafsize(curve_width)       0.5
    set grafsize(curve_color)       #000
    set grafsize(Yline_width)       1
    set grafsize(Xline_width)       1
    set grafsize(Mtick_size)        8
    set grafsize(mtick_size)        4
    set grafsize(tt_offset)         4

    set grafdata(N_text)            0
    set grafdata(XYline_textX)      0.05
    set grafdata(XYline_textY)      0.05
    set grafdata(Xtick_draw)        1
    set grafdata(Xtick_text_draw)   1
    set grafdata(Ytick_draw)        1
    set grafdata(Ytick_text_draw)   1

    #set properties(N,1) -1
    #set properties(NDM1,1) 1
    set properties(TICK1) GAMMA
    set properties(TICK2) M
    set properties(TICK3) K
    set properties(TICK4) GAMMA
    set properties(TICK5) K
    set properties(TICK6) L
    #
    # grid definition: Default: grid off
    #
    set graf_setgrid(Xmajor)       1
    set graf_setgrid(Ymajor)       1
    set graf_setgrid(Xminor)       1
    set graf_setgrid(Yminor)       1
    set graf_setgrid(Xmajor_exist) 1
    set graf_setgrid(Ymajor_exist) 1
    set graf_setgrid(Xminor_exist) 1
    set graf_setgrid(Yminor_exist) 1
    #
    # axis definition
    #
    set graf_setaxe(mirrorXaxe) 1
    set graf_setaxe(mirrorYaxe) 1
    #
    # selection
    #
    set grafselection(ButtonPressed) 0
    set grafselection(active)        0
    set grafselection(allow_move)    0
    set grafselection(it_is_curve)   0
    set grafselection(item_allowed_to_move) 0

    global prop
    if ![info exists prop(type_of_run)] {
	# due to stupid programming XCrySDen->CRYSTAL
	set prop(type_of_run) RHF
    }
}


proc GrapherSelect {gID can x y {graphcommand {}} {doubleclick {}} \
	{graphcommand {}}} {
    global grafdata grafselection graf

    set selection 0

    #
    # some cleaning up
    #
    if [info exists grafselection($gID,it_is_curve)] {
	if $grafselection($gID,it_is_curve) {
	    # make the previous selected curve back to normal
	    $can itemconfigure c$grafselection($gID,selcurve_num) \
		    -width $grafselection($gID,selcurve_width) \
		    -fill  $grafselection($gID,selcurve_fill) \
		    -stipple {}
	    set grafselection($gID,it_is_curve) 0
	}
    }
    $can dtag graf selected
    $can delete selection
    if [info exists grafselection($gID,item)] { 
	unset grafselection($gID,item) 
    }
    if [info exists grafselection($gID,id)] {
	unset grafselection($gID,id)
    }
    if [info exists grafselection($gID,X_textline_i)] {
	unset grafselection($gID,X_textline_i)
    }
    if [info exists grafselection($gID,Y_textline_i)] {
	unset grafselection($gID,Y_textline_i)
    }
    if [info exists grafselection($gID,tickXValue_g)] {
	unset grafselection($gID,tickXValue_g)
    }
    if [info exists grafselection($gID,tickYValue_g)] {
	unset grafselection($gID,tickYValue_g)
    }

    set grafselection($gID,item_allowed_to_move) 0
    set grafselection($gID,allow_move)           0
    set grafselection($gID,active)               0
    #
    # tags of Grapher items
    #

    # curve c$i
    # text  t$i
    # dummy_rectangles
    # Xline          Xl$i
    # Xline_textline Xl_tl$i
    # Xline_text     Xl_t$i
    # Yline          Yl$i
    # Yline_textline Yl_tl$i
    # Yline_text     Yl_t$i
    # xyaxis X1axe  X1$g
    # xyaxis X2axe  X2$g
    # xyaxis Y1axe  Y1$g
    # xyaxis Y2axe  Y2$g
    # xyaxis Xtitle Xti$g
    # xyaxis Ytitle Xti$g

    set id [$can find closest $x $y]
    xcDebug "GrapherSelect id:: $id"

    #
    # definitions for some items - what can be done with them
    #
    # bboxitems - bbox will be drawn around them
    # moveitems - items that are alowed to be moved
    #
    # if item is bboxitem --> we have selection
    set item_found 0
    set bboxitems {
	xyaxis text Xline_text Yline_text
	tickXValue tickYValue bar
    }
    foreach tagOrld $bboxitems {
	foreach item [$can find withtag $tagOrld] {
	    if { $item == $id } {
		$can addtag selected withtag $id
		for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g } {
		    # check if tickXValue or tickYValue
		    foreach itemid [$can find withtag tXV$g] {
			if { $itemid == $id } {
			    set grafselection($gID,tickXValue_g) $g
			    set bbox [$can bbox tXV$g]
			    $can addtag selected withtag tXV$g
			    break
			}
		    }
		    foreach itemid [$can find withtag tYV$g] {
			if { $itemid == $id } {
			    set grafselection($gID,tickYValue_g) $g
			    set bbox [$can bbox tYV$g]
			    $can addtag selected withtag tYV$g
			    break
			}
		    }		    
		    # check if item is [X|Y]line_text
		    for {set i 1} {$i <= $grafdata($gID,N_Xline,$g)} {incr i} {
			if { [$can find withtag Xl_t$g,$i] == $id } {
			    set grafselection($gID,X_textline_i) $i
			    set grafselection($gID,X_textline_g) $g
			    set bbox [$can bbox Xl_t$g,$i Xl_tl$g,$i]
			    $can addtag selected withtag l$g
			    break
			}
		    }
		    for {set i 1} {$i <= $grafdata($gID,N_Yline,$g)} {incr i} {
			if { [$can find withtag Yl_t$g,$i] == $id } {
			    set grafselection($gID,Y_textline_i) $i
			    set grafselection($gID,Y_textline_g) $g
			    set bbox [$can bbox Yl_t$g,$i Yl_tl$g,$i]
			    $can addtag selected withtag l$g
			    break
			}
		    }
		}
		# check if it is bar
		# bar can be very small and its selection will not be visible
		# - thatwhy enlarge the bounding box
		if { $grafdata($gID,type) == "BAR" } {
		    for {set i 1} {$i <= $grafdata($gID,N_point,1)} {incr i} {
			for {set j 1} {$j <= $grafdata($gID,barn,1)} {incr j} {
			    if { [$can find withtag b$j,$i] == $id } {
				if $grafdata($gID,barshadow,1) {
				    set bbox [$can bbox b$j,$i sb$j,$i]
				    $can addtag selected withtag sb$j,$i
				}
				set bbox [$can bbox b$j,$i]
				set x1 [expr [lindex $bbox 0] - 2]
				set y1 [expr [lindex $bbox 1] - 5]
				set x2 [expr [lindex $bbox 2] + 2]
				set y2 [expr [lindex $bbox 3] + 5]
				set bbox [list $x1 $y1 $x2 $y2]
				$can addtag selected withtag b$j,$i
				break
			    }
			}
		    }
		}
		# check if item is an arbitrary text
		for {set i 1} {$i <= $grafdata($gID,N_text)} {incr i} {
		    if { [$can find withtag t$i] == $id } {
			$can addtag selected withtag t_bb$i
		    }
		}
		
		if ![info exists bbox] { set bbox [$can bbox $id] }
		set x1 [lindex $bbox 0]
		set y1 [lindex $bbox 1]
		set x2 [lindex $bbox 2]
		set y2 [lindex $bbox 3]
		$can create line \
			$x1 $y1 $x1 $y2    $x1 $y2 $x2 $y2 \
			$x2 $y2 $x2 $y1    $x2 $y1 $x1 $y1 \
			-fill "#f00" \
			-width 2 \
			-stipple gray25 \
			-tags "graf selection"
		set selection 1
		set grafselection(active) 1
		if { $doubleclick != {} } {
		    BlinkingRectangle $can 3 $x1 $y1 $x2 $y2
		}
		set item_found 1
		break
	    }
	}
	if $item_found { break }
    }
    
    set move 0
    set moveitems {Xtitle Ytitle text Xline_text Yline_text}
    foreach tagOrld $moveitems {
	foreach item [$can find withtag $tagOrld] {
	    if { $item == $id } {
		set move 1
	    }
	}
    }
    if $move {
	set grafselection($gID,allow_move) 1
	set grafselection($gID,item_allowed_to_move) 1
    } 
    
    #
    # goes trough items and locate the item's ID
    # 

    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	#
	# maybe curve-tag
	for {set i 1} {$i <= $grafdata($gID,N_segment,$g)} {incr i} {
	    # one curve might be composed from several lines
	    foreach tmp_id [$can find withtag "c$g,$i"] {
		#xcDebug "GrapherSelect tmp_id:: $tmp_id"
		if { $id == $tmp_id } {
		    set grafselection($gID,id)             $id
		    set grafselection($gID,item)           c$g,$i
		    set grafselection($gID,it_is_curve)    1
		    set grafselection($gID,selcurve_num)   $g,$i
		    set grafselection($gID,selcurve_fill)  \
			    [$can itemcget $tmp_id -fill]
		    set grafselection($gID,selcurve_width) \
			    [$can itemcget $tmp_id -width]
		    $can itemconfigure c$g,$i -width 3 -fill "#f00"
		    $can raise c$g,$i curve
		    $can addtag selected withtag c$g,$i
		    if { $doubleclick != {} } {
			foreach aaa {1 2 3} {
			    $can itemconfigure c$g,$i -width 0
			    xcPause 0.1
			    $can itemconfigure c$g,$i -width 3
			    xcPause 0.1
			}
			SetCurveOptions $gID $can $i $g
		    }
		}
	    }
	}
	#
	# maybe Xtitle
	if { [$can find withtag Xti$g] == $id } {
	    if { $doubleclick != {} } {
		xcDebug "Text_fontname:: $grafdata($gID,Xtitlefont)"
		xcSetTextAtrib "Set X Title #$g" .graph$gID "Set X Title #$g" \
			grafdata($gID,X_title,$g) $can \
			grafdata($gID,Xtitlefont) Xti$g \
			grafdata($gID,Xtitle_fontsize) \
			grafdata($gID,Xtitle_fontfamily) \
			grafdata($gID,Xtitle_fontweight) \
			grafdata($gID,Xtitle_fontslant) \
			grafdata($gID,Xtitle_fontunderline) \
			grafdata($gID,Xtitle_fontoverstrike) \
			grafdata($gID,Xtitle_bbox)
		# update font
		set grafdata($gID,Xtitlefont) \
			[ModifyFont $grafdata($gID,Xtitlefont) $can \
			-size       $grafdata($gID,Xtitle_fontsize) \
			-family     $grafdata($gID,Xtitle_fontfamily) \
			-weight     $grafdata($gID,Xtitle_fontweight) \
			-slant      $grafdata($gID,Xtitle_fontslant) \
			-underline  $grafdata($gID,Xtitle_fontunderline) \
			-overstrike $grafdata($gID,Xtitle_fontoverstrike)]
	    }
	    set grafselection($gID,id)   $id
	    set grafselection($gID,item) Xti$g
	}
	#
	# maybe Ytitle
	if { [$can find withtag Yti$g] == $id } {	
	    if { $doubleclick != {} } {
		xcDebug "Text_fontname:: $grafdata($gID,Ytitlefont)"
		xcSetTextAtrib "Set Y Title #$g" .graph$gID "Set Y Title #$g" \
			grafdata($gID,Y_title,$g) $can \
			grafdata($gID,Ytitlefont) Yti$g \
			grafdata($gID,Ytitle_fontsize) \
			grafdata($gID,Ytitle_fontfamily) \
			grafdata($gID,Ytitle_fontweight) \
			grafdata($gID,Ytitle_fontslant) \
			grafdata($gID,Ytitle_fontunderline) \
			grafdata($gID,Ytitle_fontoverstrike) \
			grafdata($gID,Ytitle_bbox)
		# update font
		set grafdata($gID,Ytitlefont) \
			[ModifyFont $grafdata($gID,Ytitlefont) $can \
			-size       $grafdata($gID,Ytitle_fontsize) \
			-family     $grafdata($gID,Ytitle_fontfamily) \
			-weight     $grafdata($gID,Ytitle_fontweight) \
			-slant      $grafdata($gID,Ytitle_fontslant) \
			-underline  $grafdata($gID,Ytitle_fontunderline) \
			-overstrike $grafdata($gID,Ytitle_fontoverstrike)]
	    }
	    set grafselection($gID,id)   $id
	    set grafselection($gID,item) Yti$g	    
	}
	#
	# maybe Xaxe
	if { [$can find withtag X1a$g] == $id || \
		[$can find withtag X2a$g] == $id } {
	    if { $doubleclick != {} } {
		set graf($gID,active_graf) $g
		GrapherActivate $gID $can
		ConfigGrapher $gID $graphcommand axe
	    }
	    set grafselection($gID,id) $id
	    if { [$can find withtag X1a$g] == $id } {
		set grafselection($gID,item) X1a$g
	    } else {
		set grafselection($gID,item) X2$g
	    }
	}
	#
	# maybe Yaxe
	if { [$can find withtag Y1a$g] == $id || \
		[$can find withtag Y2a$g] == $id } {
	    if { $doubleclick != {} } {
		set graf($gID,active_graf) $g
		GrapherActivate $gID $can
		ConfigGrapher $gID $graphcommand axe
	    }
	    set grafselection($gID,id) $id
	    if { [$can find withtag Y1a$g] == $id } {
		set grafselection($gID,item) Y1a$g
	    } else {
		set grafselection($gID,item) Y2a$g
	    }
	}
    }
    #
    # maybe bar
    if { $grafdata($gID,type) == "BAR" } {
	for {set i 1} {$i <= $grafdata($gID,N_point,1)} {incr i} {
	    for {set j 1} {$j <= $grafdata($gID,barn,1)} {incr j} {
		if { [$can find withtag b$j,$i] == $id } {
		    set grafselection($gID,id)   $id
		    set grafselection($gID,item) b$j,$i
		    if { $doubleclick != {} } {
			SetBarAtrib $gID $j $i
		    }
		}
	    }
	}
    }
    #
    # maybe Xline_text
    if [info exists grafselection($gID,X_textline_i)] {
	set i $grafselection($gID,X_textline_i)
	set g $grafselection($gID,X_textline_g)
	set grafselection($gID,id)   $id
	set grafselection($gID,item) Xl_t$g,$i
	if { $doubleclick != {} } {
	    xcSetTextAtrib "Set Legend #$i,$g" \
		    .graph$gID "Set Legend #$i,$g" \
		    grafdata($gID,Xline${i}_text,$g) $can \
		    grafdata($gID,legendfont) Xl_t$g,$i \
		    grafdata($gID,legend_fontsize) \
		    grafdata($gID,legend_fontfamily) \
		    grafdata($gID,legend_fontweight) \
		    grafdata($gID,legend_fontslant) \
		    grafdata($gID,legend_fontunderline) \
		    grafdata($gID,legend_fontoverstrike) \
		    grafdata($gID,legend_bbox)
	    # update font
	    set grafdata($gID,legendfont) \
		    [ModifyFont $grafdata($gID,legendfont) $can \
		    -size       $grafdata($gID,legend_fontsize) \
		    -family     $grafdata($gID,legend_fontfamily) \
		    -weight     $grafdata($gID,legend_fontweight) \
		    -slant      $grafdata($gID,legend_fontslant) \
		    -underline  $grafdata($gID,legend_fontunderline) \
		    -overstrike $grafdata($gID,legend_fontoverstrike)]
	    xcDebug "Text_fontname:: $grafdata($gID,legendfont)"
	}
    }
    #
    # maybe Yline_text
    if [info exists grafselection($gID,Y_textline_i)] {
	set i $grafselection($gID,Y_textline_i)
	set g $grafselection($gID,Y_textline_g)
	set grafselection($gID,id)   $id
	set grafselection($gID,item) Yl_t$g,$i
	if { $doubleclick != {} } {
	    xcSetTextAtrib "Set Legend #$i,$g" \
		    .graph$gID "Set Legend #$i,$g" \
		    grafdata($gID,Yline${i}_text,$g) $can \
		    grafdata($gID,legendfont) Yl_t$g,$i \
		    grafdata($gID,legend_fontsize) \
		    grafdata($gID,legend_fontfamily) \
		    grafdata($gID,legend_fontweight) \
		    grafdata($gID,legend_fontslant) \
		    grafdata($gID,legend_fontunderline) \
		    grafdata($gID,legend_fontoverstrike) \
		    grafdata($gID,legend_bbox)
	    # update font
	    set grafdata($gID,legendfont) \
		    [ModifyFont $grafdata($gID,legendfont) $can \
		    -size       $grafdata($gID,legend_fontsize) \
		    -family     $grafdata($gID,legend_fontfamily) \
		    -weight     $grafdata($gID,legend_fontweight) \
		    -slant      $grafdata($gID,legend_fontslant) \
		    -underline  $grafdata($gID,legend_fontunderline) \
		    -overstrike $grafdata($gID,legend_fontoverstrike)]
	    xcDebug "Text_fontname:: $grafdata($gID,legendfont)"
	}
    }
    #
    # maybe tickXValue
    if [info exists grafselection($gID,tickXValue_g)] {
	set g $grafselection($gID,tickXValue_g)
	set grafselection($gID,id)   $id
	set grafselection($gID,item) tXV$g
	if { $doubleclick != {} } {
	    xcSetTextAtrib "Set Font for X-Tics #$g" \
		    .graph$gID "Set Font for X-Tics #$g" \
		    {} $can \
		    grafdata($gID,Xtickfont) tXV$g \
		    grafdata($gID,Xtick_fontsize) \
		    grafdata($gID,Xtick_fontfamily) \
		    grafdata($gID,Xtick_fontweight) \
		    grafdata($gID,Xtick_fontslant) \
		    grafdata($gID,Xtick_fontunderline) \
		    grafdata($gID,Xtick_fontoverstrike) \
		    grafdata($gID,Xtick_bbox)
	    # update font
	    set grafdata($gID,Xtickfont) \
		    [ModifyFont $grafdata($gID,Xtickfont) $can \
		    -size       $grafdata($gID,Xtick_fontsize) \
		    -family     $grafdata($gID,Xtick_fontfamily) \
		    -weight     $grafdata($gID,Xtick_fontweight) \
		    -slant      $grafdata($gID,Xtick_fontslant) \
		    -underline  $grafdata($gID,Xtick_fontunderline) \
		    -overstrike $grafdata($gID,Xtick_fontoverstrike)]
	    xcDebug "Text_fontname:: $grafdata($gID,Xtickfont)"
	}	    
    }
    #
    # maybe tickYValue
    if [info exists grafselection($gID,tickYValue_g)] {
	set g $grafselection($gID,tickYValue_g)
	set grafselection($gID,id)   $id
	set grafselection($gID,item) tYV$g
	if { $doubleclick != {} } {
	    xcSetTextAtrib "Set Font for Y-Tics #$g" \
		    .graph$gID "Set Font for Y-Tics #$g" \
		    {} $can \
		    grafdata($gID,Ytickfont) tYV$g \
		    grafdata($gID,Ytick_fontsize) \
		    grafdata($gID,Ytick_fontfamily) \
		    grafdata($gID,Ytick_fontweight) \
		    grafdata($gID,Ytick_fontslant) \
		    grafdata($gID,Ytick_fontunderline) \
		    grafdata($gID,Ytick_fontoverstrike) \
		    grafdata($gID,Ytick_bbox)
	    # update font
	    set grafdata($gID,Ytickfont) \
		    [ModifyFont $grafdata($gID,Ytickfont) $can \
		    -size       $grafdata($gID,Ytick_fontsize) \
		    -family     $grafdata($gID,Ytick_fontfamily) \
		    -weight     $grafdata($gID,Ytick_fontweight) \
		    -slant      $grafdata($gID,Ytick_fontslant) \
		    -underline  $grafdata($gID,Ytick_fontunderline) \
		    -overstrike $grafdata($gID,Ytick_fontoverstrike)]
	    xcDebug "Text_fontname:: $grafdata($gID,Ytickfont)"
	}	    
    }
    #
    # maybe arbitrary text
    for {set i 1} {$i <= $grafdata($gID,N_text)} {incr i} {
	if { [$can find withtag t$i] == $id } {
	    set grafselection($gID,id)   $id
	    set grafselection($gID,item) t$i
	    if { $doubleclick != {} } {
		xcSetTextAtrib "Set Text #$i" .graph$gID "Set Text #$i" \
			grafdata($gID,text$i) $can \
			grafdata($gID,text${i}_font) t$i \
			grafdata($gID,text${i}_fontsize) \
			grafdata($gID,text${i}_fontfamily) \
			grafdata($gID,text${i}_fontweight) \
			grafdata($gID,text${i}_fontslant) \
			grafdata($gID,text${i}_fontunderline) \
			grafdata($gID,text${i}_fontoverstrike) \
			grafdata($gID,text${i}_bbox)
		# update font
		set grafdata($gID,text${i}_font) \
			[ModifyFont $grafdata($gID,text${i}_font) $can \
			-size       $grafdata($gID,text${i}_fontsize) \
			-family     $grafdata($gID,text${i}_fontfamily) \
			-weight     $grafdata($gID,text${i}_fontweight) \
			-slant      $grafdata($gID,text${i}_fontslant) \
			-underline  $grafdata($gID,text${i}_fontunderline) \
			-overstrike $grafdata($gID,text${i}_fontoverstrike)]
		xcDebug "Text_fontname:: $grafdata($gID,text${i}_font)"
	    }
	}
    }

    #
    # this is used with conjuction of GrapherSelectMotion
    #
    set grafselection($gID,ButtonPressed) 1
    set grafselection($gID,lastX) $x
    set grafselection($gID,lastY) $y

    #
    #
    # if it was double-click do the following
    if { $doubleclick != {} } {
	set grafselection($gID,ButtonPressed) 0
	set grafselection($gID,allow_move)    0
	# now update the Grapher
	UpdateGrapher $gID $can $graphcommand
	$can raise position dummy_rectangles
    }
}


proc GrapherSelectMotion {gID can x y} {
    global grafselection grafselection grafdata graf
    
    if ![info exists grafselection($gID,allow_move)] {
	set grafselection($gID,allow_move)    0
	set grafselection($gID,ButtonPressed) 0
    }

    xcDebug "GrapherSelectMotion:: $grafselection($gID,allow_move)"
    #
    # in the upper left corner of Grapher the current position of 
    # mouse pointer will be displayed
    #
    if { ![info exists grafselection($gID,exist_position_text)] } {

	set grafselection($gID,position_text) \
		[format "(%f, %f)" \
		[XDataValue $gID $x $graf($gID,active_graf)] \
		[YDataValue $gID $y $graf($gID,active_graf)]]
	set grafselection($gID,label_window) [label $can.l -relief sunken \
		-textvariable grafselection($gID,position_text) \
		-font $grafdata($gID,deffont) \
		-width 20 \
		-bd 2]
	$can create window 10 10 \
		-anchor nw \
		-window $grafselection($gID,label_window) \
		-tags "graf position"
	set grafselection($gID,exist_position_text) 1
    } elseif { [$can find withtag position] == {} } {
	#
	# the "position" item was deleted; recreate
	#
	$can create window 10 10 \
		-anchor nw \
		-window $grafselection($gID,label_window) \
		-tags "graf position"
    } else {
	set grafselection($gID,position_text) \
		[format "(%f, %f)" \
		[XDataValue $gID $x $graf($gID,active_graf)] \
		[YDataValue $gID $y $graf($gID,active_graf)]]	
    }
    
    if { $grafselection($gID,ButtonPressed) == 1 && \
	    [info exists grafselection($gID,id)] && \
	    $grafselection($gID,allow_move) } {
	if ![info exist grafselection($gID,lastX)] {
	    set grafselection($gID,lastX) $x
	    set grafselection($gID,lastY) $y
	    return
	}
	
	set dx [expr $x - $grafselection($gID,lastX)]
	set dy [expr $y - $grafselection($gID,lastY)]
	xcDebug "(dx,dy):: $dx , $dy"
	# move selected item and its bounding box
	# maybe there are more th
	$can move selected $dx $dy
	$can move selection $dx $dy
	
	set grafselection($gID,lastX) $x
	set grafselection($gID,lastY) $y

	# for some items the current coordinates must be updated
	GrapherMoveCoor $gID $dx $dy
    }
}


proc GrapherMoveCoor {gID dx dy} {
    global grafdata grafselection

    #
    # is it X/Ytitle item
    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
	if { $grafselection($gID,item) == "Xti$g" } {	    
	    set grafdata($gID,X_title_X,$g) \
		    [expr $grafdata($gID,X_title_X,$g) + \
		    [XGraph2nMRel $gID $dx]]
	    set grafdata($gID,X_title_Y,$g) \
		    [expr $grafdata($gID,X_title_Y,$g) + $dy]
	    break
	}
	if { $grafselection($gID,item) == "Yti$g" } {
	    set grafdata($gID,Y_title_X,$g) \
		    [expr $grafdata($gID,Y_title_X,$g) + $dx]
	    set grafdata($gID,Y_title_Y,$g) \
		    [expr $grafdata($gID,Y_title_Y,$g) + \
		    [YGraph2nMRel $gID $dy]]
	    break
	}
    }
    #
    # check if item is [X|Y]line_text
    set xyline 0
    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g } { 
	for {set i 1} {$i <= $grafdata($gID,N_Xline,$g)} {incr i} {
	    if { $grafselection($gID,item) == "Xl_t$g,$i" } {
		set xyline 1
		break
	    }
	}
	for {set i 1} {$i <= $grafdata($gID,N_Yline,$g)} {incr i} {
	    if { $grafselection($gID,item) == "Yl_t$g,$i" } {
		set xyline 1
		break
	    }
	}
	if $xyline {
	    set grafdata($gID,XYline_textX,$g) \
		    [expr $grafdata($gID,XYline_textX,$g) - \
		    [XGraph2nMRel $gID $dx]]
	    set grafdata($gID,XYline_textY,$g) \
		    [expr $grafdata($gID,XYline_textY,$g) + \
		    [YGraph2nMRel $gID $dy]]
	    break
	}	
    }
    #
    # check for "text" items
    if [string match t\[0-9\]* $grafselection($gID,item)] {
	for {set i 1} {$i <= $grafdata($gID,N_text)} {incr i} {
	    if { $grafselection($gID,item) == "t$i" } {
		if [info exists grafdata($gID,text${i}_X)] {
		    xcDebug "--text recognized--"
		    set grafdata($gID,text${i}_X) \
			    [expr $grafdata($gID,text${i}_X) + \
			    [XGraph2nMRel $gID $dx]]
		    set grafdata($gID,text${i}_Y) \
			    [expr $grafdata($gID,text${i}_Y) + \
			    [YGraph2nMRel $gID $dy]]
		}
	    }   
	}
    }
}


proc GrapherSelectRelease gID {
    global grafselection
    
    set grafselection($gID,ButtonPressed) 0
    set grafselection($gID,allow_move)    0
}


proc GrapherKeyBindings {gID can what} {
   global grafselection grafdata

    xcDebug "GrapherKeyBindings"

    if ![info exists grafselection($gID,id)] { return }
    set dx 0
    set dy 0
    set id [$can find withtag selected]

    xcDebug "GrapherKeyBindings $can $what; $grafselection($gID,item)"

    switch -exact -- $what {
	left   {set dx -1}
	right  {set dx  1}
	up     {set dy -1}
	down   {set dy  1}
	delete {
	    # only arbitrary texts can be deleted so far
	    for {set i 1} {$i <= $grafdata($gID,N_text)} {incr i} {
		if { [$can find withtag t$i] == $grafselection($gID,id) } {
		    $can delete selected selection t_bb$i
		    unset grafdata($gID,text$i)
		    unset grafdata($gID,text${i}_X)
		    unset grafdata($gID,text${i}_Y)
		}
	    }
	}
	default { return }
    }   
    
    if $grafselection($gID,item_allowed_to_move) {
	$can move selected $dx $dy
	$can move selection $dx $dy
	GrapherMoveCoor $gID $dx $dy
    }
}


proc BlinkingRectangle {can ntime x1 y1 x2 y2} {

    xcDebug "BlinkingRectangle:: $can $ntime $x1 $y1 $x2 $y2"
    for {set i 0} {$i < $ntime} {incr i} {
	$can create rectangle $x1 $y1 $x2 $y2 \
		-fill "#f00" \
		-tags "blinking_rectangle" \
		-stipple gray50
	update
	xcPause 0.1
	$can delete blinking_rectangle
	update
	xcPause 0.1
    }
}
    
    
#lappend auto_path "/home/tone/src/xcrysden0.0/src"
#lappend auto_path [pwd]
#
#set system(TOPDIR) [pwd]
#set system(SCRDIR) /tmp
#set system(PWD)    [pwd]
#set system(PID)    [pid]
#set system(FORDIR) $system(TOPDIR)/F
#set system(TCLDIR) $system(TOPDIR)
#set system(BINDIR) $system(TOPDIR)
#set system(BMPDIR) $system(TOPDIR)/bitmap
#GraphInit
#
#button .b -text TONE
#entry  .e -text TONE
#set xcFonts(normal)       [lindex [.b configure -font] end]
#set xcFonts(normal_entry) [lindex [.e configure -font] end]
#set xcFonts(small)        [ModifyFontSize .b 10 \
#	    {-family helvetica -slant r -weight bold}]
#set xcFonts(small_entry)  [ModifyFontSize .e 10 \
#	    {-family helvetica -slant r -weight normal}]
#destroy .b .e
#
#trace variable grafdata(X_title_X) w PrintTrace
#trace variable grafdata(X_title_Y) w PrintTrace
#trace variable grafdata(Y_title_X) w PrintTrace
#trace variable grafdata(Y_title_Y) w PrintTrace
#
#proc PrintTrace {array elem op} {
#    set var "${array}(${elem})"
#    upvar $var varn
#
#    xcDebug "PrintTrace:: $var == $varn"
#}
#
#set prop(type_of_run) RHF
#set prop(N,1) 1
#set prop(NDM,1) 1
#set prop(N,2) 1
#set prop(NDM,2) 2
#set prop(N,3) 1
#set prop(NDM,3) 3
#
#set prop(dir)  /home/tone/pt/bulk/scanbulk/band_doss
#
##set prop(file) doss_spd.f
##DOSSGraph 3
##Grapher XYGraph   
#
#set prop(file) bwid.out
#BWIDGraph $prop(dir)/$prop(file)
#Grapher  BARGraph
#
##set prop(file) LDAVWN_pband111.f
##BANDGraph 3
##Grapher XYGraph
#
