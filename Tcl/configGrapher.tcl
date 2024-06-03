#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/configGrapher.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc ConfigGrapher {gID graphcommand what {args {}}} {
    global grafdata graf grafsize graf_setticX graf_setticY \
	    configGrapher graf_setgrid graf_setaxe graf_setdos \
	    graf_setmar system font fillEntries


    set canvas $graf($gID,can)
    set g      $graf($gID,active_graf)

    if [info exists configGrapher($gID,done)] {
	unset configGrapher($gID,done)
    }
    Prepare_and_Map_TickVar $gID
    
    switch -exact -- $what {
	ranges {set title Ranges}
	tics   {set title Tics}
	grid   {set title Grid}
	axe    {set title Axis}
	margin {set title Margin}
	dos    {set title "DOS Plot"}
	init_text {
	    set n [expr $grafdata($gID,N_text) + 1]
	    set title "Set Text #$n"
	    OneEntryToplevel [WidgetName] $title $title "${title}:" \
		    30 grafdata($gID,text$n) text 100 100
	    if { $grafdata($gID,text$n) != {} } {
		incr grafdata($gID,N_text)
		set grafdata($gID,text${n}_X) 0.5
		set grafdata($gID,text${n}_Y) 0.5
	    }
	    # now update the Grapher
	    UpdateGrapher $gID $canvas $graphcommand	    
	    return
	}
    }
    
    set t [xcToplevel [WidgetName] "Set $title for Graph #$g" \
	    "Set $title" . 100 100 1]
    # f1-frame should be flat so far (we may partition the f1-frame)
    set f1  [frame $t.f1]
    set f2  [frame $t.f2 -relief raised -bd 2]	    
    set f1l [frame $f1.l -relief raised -bd 2]
    set f1r [frame $f1.r -relief raised -bd 2]
    set ok  [DefaultButton $f2.ok -text OK \
	    -command [list ConfigGrapherOK $gID $what]]
    set can [button $f2.can -text Cancel -command [list CancelProc $t]]

    set configGrapher($gID,apply_to_all) 0
    if { $what != "margin" && $what != "dos" && $grafdata($gID,N_graf) > 1 } {
	set ckb [checkbutton $f1.cb \
		-text "Apply to all Graphs" \
		-relief raised -bd 2 \
		-variable configGrapher($gID,apply_to_all)]
    }
    switch -exact -- $what {
	ranges {
	    # set the X & Y ranges	    
	    if { $grafdata($gID,type) == "BAR" } {
		set el1 [FillEntries $f1l "{First Bar:}" \
			graf_setticX($gID,$g,min) \
			8 10 top left]
		set el2 [FillEntries $f1l "{Last Bar:}" \
			graf_setticX($gID,$g,max) \
			8 10 top left]
		set type [list intrange $grafdata($gID,firstbar,1) \
			    $grafdata($gID,lastbar,1)]	    
	    } else {
		set el1 [FillEntries $f1l "{X min:}" \
			graf_setticX($gID,$g,min) \
			5 10 top left]
		set el2 [FillEntries $f1l "{X max:}" \
			graf_setticX($gID,$g,max) \
			5 10 top left]
		set type real
	    }
	    set er1 [FillEntries $f1r "{Y min:}" \
		    graf_setticY($gID,$g,min) \
		    5 10 top left]
	    set er2 [FillEntries $f1r "{Y max:}" \
		    graf_setticY($gID,$g,max) \
		    5 10 top left]
	    set configGrapher($gID,check_num_list) \
		[list \
		     [list \
			  [list graf_setticX($gID,$g,min) $type] \
			  [list graf_setticX($gID,$g,max) $type] \
			  [list graf_setticY($gID,$g,min) real] \
			  [list graf_setticY($gID,$g,max) real] \
			 ] \
		     [list $el1 $el2 $er1 $er2] ]
	}
	tics {
	    # X tics
	    set f1lf [frame $f1l.f -relief flat -class StressText]
	    set xl [label $f1lf.l -text "X Tics" \
		    -relief groove \
		    -bd 2]
	    pack $f1lf -side top -fill both
	    pack $xl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set configGrapher($gID,check_num_listX) \
		    [SpecifyTick $gID $g $f1l graf_setticX X]
	    # Y tics
	    set f1rf [frame $f1r.f -relief flat -class StressText]
	    set yl [label $f1rf.l -text "Y Tics" \
		    -relief groove \
		    -bd 2]
	    pack $f1rf -side top -fill both
	    pack $yl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set configGrapher($gID,check_num_listY) \
		    [SpecifyTick $gID $g $f1r graf_setticY Y]
	}
	grid {
	    # X grid
	    set f1lf [frame $f1l.f -relief flat -class StressText]
	    set xl [label $f1lf.l -text "X Grid" \
		    -relief groove \
		    -bd 2]
	    pack $f1lf -side top -fill both
	    pack $xl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set XM [checkbutton $f1l.b1 -text "X major grid" \
		    -variable graf_setgrid($gID,Xmajor,$g)]
	    set Xm [checkbutton $f1l.b2 -text "X minor grid" \
		    -variable graf_setgrid($gID,Xminor,$g)]
	    pack $XM $Xm -side top -padx 5 -pady 5 -fill x
	    if { [WhatIsTickMode $gID $g X] == "position-label" } {
		set graf_setgrid($gID,Xminor,$g) 0
		$Xm configure -state disabled
	    }
	    # Y grid
	    set f1rf [frame $f1r.f -relief flat -class StressText]
	    set yl [label $f1rf.l -text "Y Grid" \
		    -relief groove \
		    -bd 2]
	    pack $f1rf -side top -fill both
	    pack $yl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set YM [checkbutton $f1r.b1 -text "Y major grid" \
		    -variable graf_setgrid($gID,Ymajor,$g)]
	    set Ym [checkbutton $f1r.b2 -text "Y minor grid" \
		    -variable graf_setgrid($gID,Yminor,$g)]
	    pack $YM $Ym -side top -padx 5 -pady 5 -fill x
	    if { [WhatIsTickMode $gID $g Y] == "position-label" } {
		set graf_setgrid($gID,Yminor,$g) 0
		$Ym configure -state disabled
	    }
	}
	axe {
	    # X axe
	    set f1lf [frame $f1l.f -relief flat -class StressText]
	    set xl [label $f1lf.l -text "X Axe" \
		    -relief groove \
		    -bd 2]
	    pack $f1lf -side top -fill both
	    pack $xl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set mX [checkbutton $f1l.b1 -text "mirror X axe" \
		    -variable graf_setaxe($gID,mirrorXaxe,$g)]
	    # Y axe
	    set f1rf [frame $f1r.f -relief flat -class StressText]
	    set yl [label $f1rf.l -text "Y Axe" \
		    -relief groove \
		    -bd 2]
	    pack $f1rf -side top -fill both
	    pack $yl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set mY [checkbutton $f1r.b1 -text "mirror Y axe" \
		    -variable graf_setaxe($gID,mirrorYaxe,$g)]
	    pack $mX $mY -side top -padx 5 -pady 5 -fill x
	}
	margin {
	    # X margins
	    set f1lf [frame $f1l.f -relief flat -class StressText]
	    set xl [label $f1lf.l -text "Set X Margins" \
		    -relief groove \
		    -bd 2]
	    pack $f1lf -side top -fill both
	    pack $xl -side top -expand 1 \
		-ipadx 10 -ipady 2 \
		-padx 10 -pady 5
	    set e1 [FillEntries $f1l {
		"Bottom X Margin:" 
		"Top X Margin:" } \
			[list graf_setmar($gID,margin_X1) \
			     graf_setmar($gID,margin_X2)] 17 10]
	    set entryList $fillEntries

	    # Y margin
	    set f1rf [frame $f1r.f -relief flat -class StressText]
	    set xl [label $f1rf.l -text "Set Y Margins" \
		    -relief groove \
		    -bd 2]
	    pack $f1rf -side top -fill both
	    pack $xl -side top -expand 1 \
		    -ipadx 10 -ipady 2 \
		    -padx 10 -pady 5
	    set e2 [FillEntries $f1r {
		"Left Y Margin:"
		"Right Y Margin:" } \
			[list graf_setmar($gID,margin_Y1) \
			     graf_setmar($gID,margin_Y2)] 17 10]
	    set entryList [concat $entryList $fillEntries]

	    #set e1t [string trimright $e1 entry1]
	    #set e2t [string trimright $e2 entry1]
	    set configGrapher($gID,check_num_list) \
		[list \
		     [list \
			  [list graf_setmar($gID,margin_X1) real] \
			  [list graf_setmar($gID,margin_X2) real] \
			  [list graf_setmar($gID,margin_Y1) real] \
			  [list graf_setmar($gID,margin_Y2) real]] \
		     $entryList]
	}
	dos {
	    #set f1lf [frame $f1l.f -relief flat -class StressText]
	    #set xl [label $f1lf.l -text "Magnify DOS Plots" \
	    #	     -relief groove \
	    #	     -bd 2]
	    #pack $f1lf -side top -fill both
	    #pack $xl -side top -expand 1 \
	    #	     -ipadx 10 -ipady 2 \
	    #	     -padx 10 -pady 5
	    pack $f1 $f2 -side top -expand 1 -fill both
	    pack $f1l -side top -expand 1 -fill both -ipadx 5 -ipady 5

	    set configGrapher($gID,check_num_list) \
		    [ScrollEntries $f1l $grafdata($gID,N_graf) \
		    {Magnify DOS Plot} \
		    {{Magnifying Factor:}} $gID,magnify posreal \
		    5 graf_setdos 0 4]
	}
	    
    }

    if { $what != "dos" } {
	pack $f1 $f2 -side top -expand 1 -fill both
	if { $what != "margin" && $grafdata($gID,N_graf) > 1 } {
	    pack $ckb -side top -expand 1 -fill both
	}
	pack $f1l $f1r -side left -expand 1 -fill both -ipadx 5 -ipady 5
    }
    if { $what == "text" } { destroy $f1r }
    pack $ok $can -side left -expand 1 -padx 10 -pady 5
    
    
    ######################################
    #                                    #
    # OK button was pressed successfully #
    #                                    #
    ######################################
    tkwait variable configGrapher($gID,done)
    xcDebug "ConfigGrapher:: DefaultButton was pressed"
    destroy $t

    #
    # taking into account apply_to_all flag 
    #
    if $configGrapher($gID,apply_to_all) {		
	for {set i 1} {$i <= $grafdata($gID,N_graf)} {incr i} {
	    append glist "$i "
	}
    } else {
	set glist $g
    }
    
    switch -exact -- $what {
	ranges {
	    foreach gg $glist {
		if { $grafdata($gID,type) == "BAR" } {
		    set b1 $graf_setticX($gID,$g,min)
		    set b2 $graf_setticX($gID,$g,max)
		    if { $b1 > $b2 } {
			set b  $b1
			set b2 $b1
			set b1 $b
		    }
		    set dx [expr 1.0 / \
			    double($grafdata($gID,N_point,1))]
		    set grafdata($gID,Xmin,$gg) \
			    [expr $dx * (double($b1) - 1.0)]
		    set grafdata($gID,Xmax,$gg) [expr $dx * double($b2)]
		} else {
		    set grafdata($gID,Xmin,$gg) $graf_setticX($gID,$g,min)
		    set grafdata($gID,Xmax,$gg) $graf_setticX($gID,$g,max)
		}
		set grafdata($gID,Ymin,$gg) $graf_setticY($gID,$g,min)
		set grafdata($gID,Ymax,$gg) $graf_setticY($gID,$g,max)
	    }
	}
	tics {
	    foreach gg $glist {
		set grafdata($gID,Xtick_draw,$gg)      $graf_setticX($gID,$g,tick_draw)     
		set grafdata($gID,Ytick_draw,$gg)      $graf_setticY($gID,$g,tick_draw)     
		set grafdata($gID,Xtick_text_draw,$gg) $graf_setticX($gID,$g,tick_text_draw)
		set grafdata($gID,Ytick_text_draw,$gg) $graf_setticY($gID,$g,tick_text_draw)
	    }
	    #
	    # X tics
	    #
	    if { $graf_setticX($gID,$g,mode) == "start-increment-end" } {
		foreach gg $glist {
		    set wdone 1
		    set i 0
		    set grafdata($gID,N_MXtick,$gg) 0
		    while $wdone {
			set x [expr $graf_setticX($gID,$g,start) + \
				double($i) * \
				$graf_setticX($gID,$g,increment)]
			if [IsLEQ 0.0001 $x $graf_setticX($gID,$g,end)] {
			    incr i
			    incr grafdata($gID,N_MXtick,$gg)
			    set grafdata($gID,Xtick$i,$gg) $x
			    set grafdata($gID,Xtick${i}_text,$gg) \
				    [TickFormat $x \
				    $graf_setticX($gID,$g,increment)]
			} else {
			    set wdone 0
			    break
			}
		    }
		    set grafdata($gID,N_mXtick,$gg) \
			    $graf_setticX($gID,$g,N_mtick)
		}
	    } else { 
		# if label-position are manually specified, 
		# the number of minor tics is 0
		foreach gg $glist {
		    set grafdata($gID,N_mXtick,$gg) 0
		    set grafdata($gID,N_MXtick,$gg) \
			    $graf_setticX($gID,$g,N_Mtick)
		    for {set i 1} {$i <= $graf_setticX($gID,$g,N_Mtick)} \
			    {incr i} {
			set grafdata($gID,Xtick$i,$gg) \
				$graf_setticX($gID,$g,tick,$i)     
			set grafdata($gID,Xtick${i}_text,$gg) \
				$graf_setticX($gID,$g,tick_text,$i)
		    }
		}
	    }
	    #
	    # Y tics
	    #
	    if { $graf_setticY($gID,$g,mode) == "start-increment-end" } {
		foreach gg $glist {
		    set wdone 1
		    set i 0
		    set grafdata($gID,N_MYtick,$gg) 0
		    while $wdone {
			set y [expr $graf_setticY($gID,$g,start) + \
				double($i) * \
				$graf_setticY($gID,$g,increment)]
			if [IsLEQ 0.0001 $y $graf_setticY($gID,$g,end)] {
			    incr i
			    incr grafdata($gID,N_MYtick,$gg)
			    set grafdata($gID,Ytick$i,$gg) $y
			    set grafdata($gID,Ytick${i}_text,$gg) \
				    [TickFormat $y \
				    $graf_setticY($gID,$g,increment)]
			} else {
			    set wdone 0
			    break
			}
		    }
		    set grafdata($gID,N_mYtick,$gg) \
			    $graf_setticY($gID,$g,N_mtick)
		}
	    } else {
		# if label-position are manually specified, 
		# the number of minor tics is 0
		foreach gg $glist {
		    set grafdata($gID,N_mYtick,$gg) 0
		    set grafdata($gID,N_MYtick,$gg) \
			    $graf_setticY($gID,$g,N_Mtick)
		    for {set i 1} {$i <= $graf_setticY($gID,$g,N_Mtick)} \
			    {incr i} {
			set grafdata($gID,Ytick$i,$gg) \
				$graf_setticY($gID,$g,tick,$i)     
			set grafdata($gID,Ytick${i}_text,$gg) \
				$graf_setticY($gID,$g,tick_text,$i)
		    }
		}
	    }
	}
	grid {
	    if $configGrapher($gID,apply_to_all) {
		foreach gg $glist {
		    set graf_setgrid($gID,Xmajor,$gg) \
			    $graf_setgrid($gID,Xmajor,$g)
		    set graf_setgrid($gID,Xminor,$gg) \
			    $graf_setgrid($gID,Xminor,$g)
		    set graf_setgrid($gID,Ymajor,$gg) \
			    $graf_setgrid($gID,Ymajor,$g)
		    set graf_setgrid($gID,Yminor,$gg) \
			    $graf_setgrid($gID,Yminor,$g)
		}
	    }
	}
	margin {
	    set grafsize($gID,margin_X1) $graf_setmar($gID,margin_X1)
	    set grafsize($gID,margin_X2) $graf_setmar($gID,margin_X2)
	    set grafsize($gID,margin_Y1) $graf_setmar($gID,margin_Y1)
	    set grafsize($gID,margin_Y2) $graf_setmar($gID,margin_Y2)
	}
	dos {
	    for {set g 1} {$g <= $grafdata($gID,N_graf)} {incr g} {
		for {set i 1} {$i <= $grafdata($gID,N_segment,$g)} {incr i} {
		    xcDebug "Magnify:: $graf_setdos($gID,magnify,$g)"
		    for {set j 1} {$j <= $grafdata($gID,N_point,$g)} {incr j} {
			set grafdata($gID,$j,$i,$g) \
				[expr $grafdata($gID,$j,$i,$g) * \
				$graf_setdos($gID,magnify,$g) / \
				$grafdata($gID,magnify$g)]
		    }
		}
		set grafdata($gID,magnify$g) $graf_setdos($gID,magnify,$g)
	    }	    	    
	}
	axe {
	    foreach gg $glist {
		set graf_setaxe($gID,mirrorXaxe,$gg) \
			$graf_setaxe($gID,mirrorXaxe,$g)
		set graf_setaxe($gID,mirrorYaxe,$gg) \
			$graf_setaxe($gID,mirrorYaxe,$g)
	    }
	}
    }

    # now update the Grapher
    UpdateGrapher $gID $canvas $graphcommand
}


proc ConfigGrapherOK {gID what} {
    global configGrapher grafdata graf_setticX graf_setticY graf

    set g   $graf($gID,active_graf)    

    if [string match *${what}* {grid axe text}] {
	set configGrapher($gID,done) 1
	return
    }

    switch -exact -- $what {
	tics {
	    # we have configGrapher(check_num_listX) & 
	    # configGrapher(check_num_listY)
	    # the format of configGrapher(check_num_listX) is the following:
	    #
	    #      configGrapher(check_num_listX) -> {varlist2 foclist2}
	    #
	    #      varlist2 -> {varlist_a varlist_b}
	    #      foclist2 -> {foclist_a foclist_b}
	    set Xvarlist2 [lindex $configGrapher($gID,check_num_listX) 0]
	    set Xfoclist2 [lindex $configGrapher($gID,check_num_listX) 1]
	    set Yvarlist2 [lindex $configGrapher($gID,check_num_listY) 0]
	    set Yfoclist2 [lindex $configGrapher($gID,check_num_listY) 1]
	    # X list
	    if { $graf_setticX($gID,$g,mode) == "start-increment-end" } {
		set indX 0
	    } else {
		set indX 1
	    }
	    set varlist [lindex $Xvarlist2 $indX]
	    set foclist [lindex $Xfoclist2 $indX]
	    if ![check_var $varlist $foclist] {
		return
	    }
	    # Y list
	    if { $graf_setticY($gID,$g,mode) == "start-increment-end" } {
		set indY 0
	    } else {
		set indY 1
	    }
	    set varlist [lindex $Yvarlist2 $indY]
	    set foclist [lindex $Yfoclist2 $indY]
	    if ![check_var $varlist $foclist] {
		return
	    }	 
	}
	default {
	    # the format of configGrapher(check_num_list) is the following:
	    #
	    #      configGrapher(check_num_list) -> {varlist foclist}
	    if ![eval check_var $configGrapher($gID,check_num_list)] {
		return
	    }
	}
    }

    set configGrapher($gID,done) 1
}


proc Prepare_and_Map_TickVar gID {
    global grafdata graf_setticX graf_setticY graf_setmar grafsize \
	    graf_setdos graf

    # WARNING::
    # because of some limitation the syntax of graf_setXXX variable will be:
    # graf_setXXX($gID,$g,*) !!!!!!!!!!!!!
    #
    # with some exceptions: grafdata($gID,magnify$g)

    xcDebug "Prepare_and_Map_TickVar #1"
    
    set g $graf($gID,active_graf)

    if { $grafdata($gID,type) == "BAR" } {
	set graf_setticX($gID,$g,min) \
		[expr round(double($grafdata($gID,N_point,1)) * \
		$grafdata($gID,Xmin,$g)) + 1]
	set graf_setticX($gID,$g,max) \
		[expr round(double($grafdata($gID,N_point,1)) * \
		$grafdata($gID,Xmax,$g))]
    } else {
	set graf_setticX($gID,$g,min)   $grafdata($gID,Xmin,$g)
	set graf_setticX($gID,$g,max)   $grafdata($gID,Xmax,$g)
    }
    set graf_setticY($gID,$g,min)   $grafdata($gID,Ymin,$g)
    set graf_setticY($gID,$g,max)   $grafdata($gID,Ymax,$g)

    set graf_setticX($gID,$g,N_Mtick) $grafdata($gID,N_MXtick,$g)
    set graf_setticX($gID,$g,N_mtick) $grafdata($gID,N_mXtick,$g)	
    set graf_setticY($gID,$g,N_Mtick) $grafdata($gID,N_MYtick,$g)	
    set graf_setticY($gID,$g,N_mtick) $grafdata($gID,N_mYtick,$g)

    set graf_setticX($gID,$g,tick_draw)      $grafdata($gID,Xtick_draw,$g)
    set graf_setticY($gID,$g,tick_draw)      $grafdata($gID,Ytick_draw,$g)
    set graf_setticX($gID,$g,tick_text_draw) $grafdata($gID,Xtick_text_draw,$g)
    set graf_setticY($gID,$g,tick_text_draw) $grafdata($gID,Ytick_text_draw,$g)
    
    set graf_setmar($gID,margin_X1) $grafsize($gID,margin_X1)
    set graf_setmar($gID,margin_X2) $grafsize($gID,margin_X2)
    set graf_setmar($gID,margin_Y1) $grafsize($gID,margin_Y1)
    set graf_setmar($gID,margin_Y2) $grafsize($gID,margin_Y2)

    for {set i 1} {$i <= $grafdata($gID,N_graf)} {incr i} {
	if ![info exists grafdata($gID,magnify$i)] {
	    set grafdata($gID,magnify$i) 1
	}
	set graf_setdos($gID,magnify,$i) $grafdata($gID,magnify$i)
    }

    WhatIsTickMode $gID $g X
    WhatIsTickMode $gID $g Y	    
    xcDebug "Prepare_and_Map_TickVar"

}


proc SpecifyTick {gID g w array who} {
    # w ....... name of parent
    # array ... $array == global arrayname
    # who ..... X or Y
    global graf_setticX graf_setticY
    if { $array == "graf_setticX" } {
	set mode                graf_setticX($gID,$g,mode)
	set start               graf_setticX($gID,$g,start)
	set increment           graf_setticX($gID,$g,increment)
	set end                 graf_setticX($gID,$g,end)
	set N_Mtick             graf_setticX($gID,$g,N_Mtick)
	set N_mtick             graf_setticX($gID,$g,N_mtick)
	set tick_draw           graf_setticX($gID,$g,tick_draw)
	set tick_text_draw      graf_setticX($gID,$g,tick_text_draw)
    } elseif { $array == "graf_setticY"} {
	set mode                graf_setticY($gID,$g,mode)
	set start               graf_setticY($gID,$g,start)
	set increment           graf_setticY($gID,$g,increment)
	set end                 graf_setticY($gID,$g,end)
	set N_Mtick             graf_setticY($gID,$g,N_Mtick)
	set N_mtick             graf_setticY($gID,$g,N_mtick)
	set tick_draw           graf_setticY($gID,$g,tick_draw)
	set tick_text_draw      graf_setticY($gID,$g,tick_text_draw)
    }

    xcDebug "SpecifyTick:: $array"

    set f1lu [frame $w.u -relief groove -bd 2]
    set r1 [radiobutton $f1lu.r1 \
	    -text "\"start-increment-end\"" \
	    -variable $mode \
	    -value start-increment-end \
	    -anchor w]
    pack $r1 -side top
    set elu1 [FillEntries $f1lu "{Start:}" \
	    $start \
	    17 10 top left]
    set elu2 [FillEntries $f1lu "{Increment:}" \
	    $increment \
	    17 10 top left]
    set elu3 [FillEntries $f1lu "{End:}" \
	    $end \
	    17 10 top left]
    set elu4 [FillEntries $f1lu "{No. of Minor Tics:}" \
	    $N_mtick \
	    17 10 top left]

    set f1lb [frame $w.b -relief groove -bd 2]
    set r2 [radiobutton $f1lb.r1 \
	    -text "\"position & label\"" \
	    -variable $mode \
	    -value position-label \
	    -anchor w]
    pack $r2 -side top
    set elb1 [FillEntries $f1lb "{Number of Tics:}" \
	    $N_Mtick \
	    15 10 top left]
    set blb1 [button $f1lb.b -text "Specify Tics & Labels" \
	    -command [list SpecifyTickLabels $gID $g $who]]

    pack $f1lu $f1lb -side top -fill both -expand 1 \
	    -padx 5 -pady 5 -ipadx 5 -ipady 5
    pack $blb1 -side top -expand 1 -padx 5 -pady 5
    
    set f1lbb [frame $w.bb -relief groove -bd 2]
    set cb1 [checkbutton $f1lbb.cb1 \
	    -text "display tics" \
	    -variable $tick_draw \
	    -onvalue 1 \
	    -offvalue 0 \
	    -anchor w]
    set cb2 [checkbutton $f1lbb.cb2 \
	    -text "display tics labels" \
	    -variable $tick_text_draw \
	    -onvalue 1 \
	    -offvalue 0 \
	    -anchor w]
    
    pack $f1lbb -side top -fill both -expand 1 \
	    -padx 5 -pady 5 -ipadx 5 -ipady 5
    pack $cb1 $cb2 -side top -fill both -expand 1

    #
    # make diable/enable entries control
    #
    $r1 configure -command \
	    [list SpecifyTicsState [list $elu1 $elu2 $elu3 $elu4] \
	    $elb1 $blb1 0]
    $r2 configure -command \
	    [list SpecifyTicsState $elb1 [list $elu1 $elu2 $elu3 $elu4] \
	    $blb1 1]

    xcDebug [list SpecifyTicsState [list $elu1 $elu2 $elu3 $elu4] $elb1 $blb1]
    upvar $mode m
    if { $m == "position-label" } {
	SpecifyTicsState $elb1 [list $elu1 $elu2 $elu3 $elu4] $blb1 1
    } else {
	SpecifyTicsState [list $elu1 $elu2 $elu3 $elu4] $elb1 $blb1 0
    }

    #
    # prepare for checking the variables
    #
    set varlist [list \
	    [list\
	    [list $start     real] \
	    [list $increment real] \
	    [list $end       real] \
	    [list $N_mtick   int]] \
	    [list \
	    [list $N_Mtick   int]]]

    set foclist [list \
	    [list $elu1 $elu2 $elu3 $elu4] \
	    [list $elb1]]

    return [list $varlist $foclist]
}


proc SpecifyTicsState {active disabled but butstate} {

    if $butstate  {
	$but configure -state normal
    } else {
	$but configure -state disabled
    }
    set disabled_color [lindex \
	    [GetWidgetConfig button -disabledforeground] end]
    set enabled_color  [lindex \
	    [GetWidgetConfig button -foreground] end]
    
    foreach wid $active {
	xcDebug "SpecifyTicsState_active:: $wid"
	# enable entry widget
	$wid configure -relief sunken -state normal
	# "enable label"
	set lab [string trimright $wid entry1]lab1
	$lab configure -foreground $enabled_color
    }

    foreach wid $disabled {
	xcDebug "SpecifyTicsState_passive:: $wid"
	# disable entry widget
	$wid configure -relief flat -state disabled
	# "disenable label"
	set lab [string trimright $wid entry1]lab1
	$lab configure -foreground $disabled_color
    }
}


proc SpecifyTickLabels {gID g who} {
    global graf_settic graf_setticX graf_setticY grafdata

    xcDebug "SpecifyTickLabels"
    set t [xcToplevel [WidgetName] "${who}-Set Tics & Labels" \
	    "${who}-Set Tics & Labels" . 100 100 1]
    set ft [frame $t.ft -relief raised -bd 2]
    set fb [frame $t.fb -relief raised -bd 2]
    pack $ft $fb -side top -expand 1 -fill both 

    if { $who == "X" } {
	if { $graf_setticX($gID,$g,N_Mtick) == 0 } { 
	    destroy $t
	    return 
	}
	set graf_settic($gID,check_num_list) \
		[ScrollEntries $ft $graf_setticX($gID,$g,N_Mtick) \
		"X-Tic & Label n.:" "${who}-Position: Label:" \
		[list $gID,$g,tick $gID,$g,tick_text] {real text} \
		10 graf_setticX 0 4]
    } elseif { $who == "Y" } {
	if { $graf_setticY($gID,$g,N_Mtick) == 0 } { 
	    destroy $t
	    return 
	}
	set graf_settic($gID,check_num_list) \
		[ScrollEntries $ft $graf_setticY($gID,$g,N_Mtick) \
		"Y-Tic & Label n.:" "${who}-Position: Label:" \
		[list $gID,$g,tick $gID,$g,tick_text] {real text} \
		10 graf_setticY 0 4]
    }

    set ok  [DefaultButton $fb.ok -text OK \
	    -command [list SpecifyTickLabelsOK $gID $t]]
    set can [button $fb.can -text Cancel -command [list CancelProc $t]]
    pack $ok $can -side left -expand 1 -padx 10 -pady 5
}


proc SpecifyTickLabelsOK {gID t} {
    global graf_settic
    
    if ![eval check_var $graf_settic($gID,check_num_list)] {
	return
    }
    
    destroy $t
}


proc WhatIsTickMode {gID g axis} {
    global grafdata graf_setticX graf_setticY
    
    if { $axis == "X" } {
	if [info exists graf_setticX($gID,$g,mode)] {
	    return $graf_setticX($gID,$g,mode)
	} else {
	    # now find out if tics are specified in start-increment-end mode
	    
	    #
	    # first X tics
	    #
	    set start_increment_end 1    
	    # maybe grafdata(N_mXtics) == 0, but still the tics are 
	    # specified in start-increment-end mode
	    if { $grafdata($gID,N_MXtick,$g) > 1 } {
		set dXold [expr \
			$grafdata($gID,Xtick2,$g) - $grafdata($gID,Xtick1,$g)]
		for {set i 3} {$i <= $grafdata($gID,N_MXtick,$g)} {incr i} {
		    set ii [expr $i - 1]
		    set dXnew \
			    [expr $grafdata($gID,Xtick$i,$g) - \
			    $grafdata($gID,Xtick$ii,$g)]
		    if ![IsEqual 0.0001 $dXnew $dXold] {
			set start_increment_end 0
			break
		    }
		    set dXold $dXnew
		}
		# 
		# make one more check
		#
		for {set i 1} {$i <= $grafdata($gID,N_MXtick,$g)} {incr i} {
		    # tics value and text must be the same
		    # for position-increment-end mode
		    if ![catch {expr abs($grafdata($gID,Xtick${i}_text,$g))}] {
			# $grafdata($gID,Xtick${i}_text,$g) is a number
			if ![IsEqual 0.001 $grafdata($gID,Xtick$i,$g) \
				$grafdata($gID,Xtick${i}_text,$g)] {
			    set start_increment_end 0
			    break
			}
		    }
		}
		
		
		set graf_setticX($gID,$g,mode) "position-label"
		if $start_increment_end {
		    set graf_setticX($gID,$g,mode)  "start-increment-end"
		    set graf_setticX($gID,$g,start) $grafdata($gID,Xtick1,$g)
		    set graf_setticX($gID,$g,increment) \
			    [expr $grafdata($gID,Xtick2,$g) - \
			    $grafdata($gID,Xtick1,$g)]
		    set last $grafdata($gID,N_MXtick,$g)
		    set graf_setticX($gID,$g,end) $grafdata($gID,Xtick$last,$g)
		}
		for {set i 1} {$i <= $grafdata($gID,N_MXtick,$g)} {incr i} {
		    set graf_setticX($gID,$g,tick,$i) \
			    $grafdata($gID,Xtick$i,$g)
		    set graf_setticX($gID,$g,tick_text,$i) \
			    $grafdata($gID,Xtick${i}_text,$g)
		}
	    } else {
		set graf_setticX($gID,$g,mode)  "position-label"
		if [info exists grafdata($gID,Xtick1,$g)] {
		    set graf_setticX($gID,$g,tick,1) $grafdata($gID,Xtick1,$g)
		    set graf_setticX($gID,$g,tick_text,1) \
			    $grafdata($gID,Xtick1_text,$g)
		}
	    }
	}
	return $graf_setticX($gID,$g,mode)
    }

    if { $axis == "Y" } {
	if [info exists graf_setticY($gID,$g,mode)] {
	    return $graf_setticY($gID,$g,mode)
	} else {
	    # now find out if tics are specified in start-increment-end mode
	    
	    #
	    # now, Y tics
	    #
	    set start_increment_end 1    
	    # maybe grafdata(N_mYtics) == 0, but still the tics are 
	    # specified in start-increment-end mode
	    if { $grafdata($gID,N_MYtick,$g) > 1 } {
		set dYold [expr \
			$grafdata($gID,Ytick2,$g) - $grafdata($gID,Ytick1,$g)]
		for {set i 3} {$i <= $grafdata($gID,N_MYtick,$g)} {incr i} {
		    set ii [expr $i - 1]
		    set dYnew [expr $grafdata($gID,Ytick$i,$g) - \
			    $grafdata($gID,Ytick$ii,$g)]
		    if ![IsEqual 0.0001 $dYnew $dYold] {
			set start_increment_end 0
			break
		    }
		    set dYold $dYnew
		}
		# 
		# make one more check
		#
		for {set i 1} {$i <= $grafdata($gID,N_MYtick,$g)} {incr i} {
		    # tics value and text must be the same
		    # for position-increment-end mode
		    if ![catch {expr abs($grafdata($gID,Ytick${i}_text,$g))}] {
			# $grafdata($gID,Ytick${i}_text,$g) is a number
			if ![IsEqual 0.001 $grafdata($gID,Ytick$i,$g) \
				$grafdata($gID,Ytick${i}_text,$g)] {
			    set start_increment_end 0
			    break
			}
		    }
		}
		
		set graf_setticY($gID,$g,mode) "position-label"
		if $start_increment_end {
		    set graf_setticY($gID,$g,mode)  "start-increment-end"
		    set graf_setticY($gID,$g,start) $grafdata($gID,Ytick1,$g)
		    set graf_setticY($gID,$g,increment) \
			    [expr $grafdata($gID,Ytick2,$g) - \
			    $grafdata($gID,Ytick1,$g)]
		    set last $grafdata($gID,N_MYtick,$g)
		    set graf_setticY($gID,$g,end) $grafdata($gID,Ytick$last,$g)
		}
		for {set i 1} {$i <= $grafdata($gID,N_MYtick,$g)} {incr i} {
		    set graf_setticY($gID,$g,tick,$i) \
			    $grafdata($gID,Ytick$i,$g)
		    set graf_setticY($gID,$g,tick_text,$i) \
			    $grafdata($gID,Ytick${i}_text,$g)
		}
	    } else {
		set graf_setticY($gID,$g,mode)  "position-label"
		if [info exists grafdata($gID,Ytick1)] {
		    set graf_setticY($gID,$g,tick,1) $grafdata($gID,Ytick1,$g)
		    set graf_setticY($gID,$g,tick_text,1) \
			    $grafdata($gID,Ytick1_text,$g)
		}
	    }
	}
	return $graf_setticY($gID,$g,mode)
    }
}


proc SetCurveOptions {gID can i g} {
    global grafsize done grafselection graf
    # i ... curve number

    if ![info exists grafsize($gID,curve_width_$i,$g)] {
	set grafsize($gID,curve_width_$i,$g) $grafsize(curve_width)
    }
    if ![info exists grafsize($gID,curve_color_$i,$g)] {
	set grafsize($gID,curve_color_$i,$g) $grafsize(curve_color)
    }

    set tcan [winfo toplevel $can]
    set t [xcToplevel [WidgetName] \
	    "Settings of Curve #$i" "Settings of Curve #$i" $tcan 0 0 1]
    set f1 [frame $t.up -relief raised -bd 2]
    set f2 [frame $t.dn -relief raised -bd 2]
    FillEntries $f1 \
	    [list "Width of Curve #$i:" "Color of Curve #$i in #RGB:"] \
	    [list grafsize($gID,curve_width_$i,$g) \
	    grafsize($gID,curve_color_$i,$g)] \
	    24 10 top left
    set ok [button $f2.ok -text OK \
	    -command [list SetCurveOptionsOK $gID $g $i]]
    bind $ok <Return> [list SetCurveOptionsOK $gID $g $i]
    bind $f1.f1.2.entry2 <Return> [list SetCurveOptionsOK $gID $g $i]

    pack $f1 $f2 -side top -fill both -ipady 5
    pack $ok -expand 1 -padx 5

    tkwait variable done
    destroy $t
    $can itemconfigure c$i \
	    -width $grafsize($gID,curve_width_$i,$g) \
	    -fill  $grafsize($gID,curve_color_$i,$g)

    # update grafselection(curve_*)
    set grafselection($gID,selcurve_fill)  $grafsize($gID,curve_color_$i,$g)
    set grafselection($gID,selcurve_width) $grafsize($gID,curve_width_$i,$g)
}


proc SetCurveOptionsOK {gID g i} {
    global grafsize done

    set err 1
    if ![number grafsize($gID,curve_width_$i,$g) posreal] {
	set err 0
    }

    set lenok 0
    set length [string length $grafsize($gID,curve_color_$i,$g)]
    foreach len {4 7 10 13} {
	if { $len == $length} {
	    set lenok 1
	}
    }
    set cok 1
    for {set ii 1} {$ii < $length} {incr ii} {
	set c [string index $grafsize($gID,curve_color_$i,$g) $ii]
	if ![string match \[0-9a-f\] $c] {
	    set cok 0
	}
    }
    if { $lenok == 0 || $cok == 0 } {
	tk_dialog [WidgetName] "ERROR" "ERROR: You have badly specified the color of curve $i. Should be specified as #RGB" error 0 OK
	set err 0
    }
    
    if $err {
	set done 1
    }
}


proc xcSetTextAtrib { title master labeltext entryvar canvas \
	 font_name font_item font_size font_family \
	 font_weight font_slant font_underline font_overstrike bboxvar} {
    global system stA

    # can be more xcSetTextAtrib widows open at ones
    if ![info exists stA(id)] {
	 set stA(id) 1
    } else {
	 incr stA(id)
    }

    set id $stA(id)

    set tplw [xcToplevel [WidgetName] $title {Set Text} $master 100 100]

    set f1   [frame $tplw.f1 -relief raised -bd 2]
    set f2   [frame $tplw.f2 -relief raised -bd 2]
    pack $f1 $f2 -side top -expand 1 -fill both

    set f1lf [frame $f1.f1 -relief flat -class StressText]
    set xl [label $f1lf.l -text $labeltext -relief groove -bd 2]
    pack $f1lf -side top -fill both
    pack $xl -side top -expand 1 -ipadx 10 -ipady 2 -padx 10 -pady 5
    
    set f1l1 [frame $f1.1]
    set f1l2 [frame $f1.2]
    set f1l3 [frame $f1.3]
    pack $f1l1 $f1l2 $f1l3 -side top -fill both
    
    #
    # query font
    #
    set stA($id,fontname) [$canvas itemcget $font_item -font]
    xcDebug "stA($id,fontname):: $stA($id,fontname) $font_item"
    set stA($id,fontsize) [GetFontAtribute $stA($id,fontname) $canvas -size]


    #
    # frame #1
    #
    if { $entryvar != {} } {
	set e [Entries $f1l1 Text: $entryvar 30]
	focus $e
    }

    #
    # frame #2
    #
    xcMenuEntry $f1l2 {Font Size:} 10 \
	     stA($id,fontsize) \
	     {6 8 10 12 14 16 18 20 24} \
	     -entrystate disabled \
	     -labelwidth 11 \
	     -labelanchor w
    set ck [checkbutton $f1l2.b1 \
	     -text "Bounding Box" \
	     -variable $bboxvar]
	     
    pack $ck -side top -padx 5 -pady 5 -fill x
    
    #
    # frame #3
    #
    set stA($id,fontfamily) [GetFontAtribute $stA($id,fontname) $canvas -family]

    xcMenuEntry $f1l3 {Font Family:} 10 \
	     stA($id,fontfamily) \
	     {times helvetica courier} \
	     -entrystate disabled \
	     -labelwidth 11 \
	     -labelanchor w
    
    if { [GetFontAtribute $stA($id,fontname) $canvas -weight] == "normal" } {
	 set stA($id,fontbold) 0
    } else {
	 set stA($id,fontbold) 1
    }
    if { [GetFontAtribute $stA($id,fontname) $canvas -slant] == "roman" } {
	 set stA($id,fontitalic) 0
    } else {
	 set stA($id,fontitalic) 1
    }	    
    set stA($id,fontunderline)  \
	     [GetFontAtribute $stA($id,fontname) $canvas -underline]  
    set stA($id,fontoverstrike) \
	     [GetFontAtribute $stA($id,fontname) $canvas -overstrike]

    xcCheckButtonRow $f1l3 4 [list \
	     @$system(BMPDIR)/bold.xbm \
	     @$system(BMPDIR)/italic.xbm \
	     @$system(BMPDIR)/underline.xbm \
	     @$system(BMPDIR)/overstrike.xbm] \
	     [list \
	     stA($id,fontbold) \
	     stA($id,fontitalic) \
	     stA($id,fontunderline) \
	     stA($id,fontoverstrike)] \
	     [list \
	     xcCheckButtonDummy xcCheckButtonDummy \
	     xcCheckButtonDummy xcCheckButtonDummy]

    set ok  [DefaultButton $f2.ok -text OK \
	     -command [list xcSetTextAtribBut $id ok]]
    set can [button $f2.can -text Cancel \
	     -command [list xcSetTextAtribBut $id can $tplw]]
    pack $ok $can -side left -expand 1 -padx 10 -pady 5

    tkwait variable stA($id,done)
    upvar #0 \
	     $font_name       fontname      \
	     $font_size       fontsize      \
	     $font_family     fontfamily    \
	     $font_weight     fontweight    \
	     $font_slant      fontslant     \
	     $font_underline  fontunderline \
	     $font_overstrike fontoverstrike

    if $stA($id,fontbold) {
	set weight bold
    } else {
	set weight normal
    }
    if $stA($id,fontitalic) {
	set slant italic
    } else {
	set slant roman
    }

    set fontname       $stA($id,fontname)
    set fontsize       $stA($id,fontsize)     
    set fontfamily     $stA($id,fontfamily)    
    set fontweight     $weight              
    set fontslant      $slant     
    set fontunderline  $stA($id,fontunderline) 
    set fontoverstrike $stA($id,fontoverstrike)
    #
    # now unset the stA($id,*) array elements
    #
    foreach elem [array names stA $id,*] {
	unset stA($elem)
    }
    destroy $tplw
}


proc xcSetTextAtribBut {id what {top {}}} {
    global stA

    if { $what == "can" } {
	foreach elem [array names stA $id,*] {	    
	    unset stA($elem)
	}
	destroy $top
    }
    
    if { $what == "ok" } {
	set stA($id,done) 1
    }
}


proc SetBarAtrib {gID j i} {
    global grafdata grafsize bar_atrib_color mody_col
    # j ... alpha/beta (for open-shell)
    # i ... which bar

    set t [xcToplevel [WidgetName] "Set Bar #$i Attributes" "Set Bar" \
	    .graph$gID 100 100 1]
    catch { grab $t }

    set bar_atrib_color($gID,apply_to_all) 0
    set ckb [checkbutton $t.cb \
	    -text "Apply to all bars" \
	    -relief raised -bd 2 \
	    -variable bar_atrib_color($gID,apply_to_all)]


    set f1 [frame $t.f1 -relief flat -class StressText]
    set xl [label $f1.l -text "Set Bar Attributes" \
	    -relief groove \
	    -bd 2]
	
    pack $ckb  -side top -fill both -expand 1
    pack $f1 -side top -fill both -expand 1
    pack $xl -side top -expand 1 -padx 10 -pady 5 -ipadx 10 -ipady 2

    set f2 [frame $t.f2]
    set f3 [frame $t.f3]
    set f4 [frame $t.f4]
    set f5 [frame $t.f5]
    set f6 [frame $t.f6]
    pack $f2 $f3 $f4 $f5 $f6 -side top -fill both -expand 1

    # FRAME #2
    set command [list \
	    [list {Fill Color:   } [list SetBarAtribColor $gID $j $i FILL]] \
	    [list {Outline Color:} [list SetBarAtribColor $gID $j $i OUTLINE]] \
	    [list {Shadow Color: } [list SetBarAtribColor $gID $j $i SHADOW]]]
    MultiWidget $f2 -b_height 2 -testbutton 1 \
	    -create_tplw 0 \
	    -command $command
    
    # FRAME #3
    xcMenuEntry $f3 {Outline Width:} 11 \
	    grafsize($gID,bar${i}_width,1) \
	    {0 1 2 3 4 5 6 7 8 9} \
	    -entrystate disabled \
	    -labelwidth 13 \
	    -labelanchor w

    if { $grafsize($gID,bar${i}_stipple,$j) == {} } {
	set grafsize($gID,bar${i}_stipple,$j) none
    }

    # FRAME #4
    xcMenuEntry $f4 {Stipple Pattern:} 11 \
	    grafsize($gID,bar${i}_stipple,$j) \
	    { none gray75 gray50 gray25 gray12 } \
	    -entrystate disabled \
	    -labelwidth 13 \
	    -labelanchor w

    # FRAME #5    
    set bar_atrib_color($gID,xoffset) \
	    [expr round((1.0 - $grafdata($gID,Xoffset,1)) * 100)]
    xcMenuEntry $f5 {Bar Width (%):} 11 \
	    bar_atrib_color($gID,xoffset) \
	    { 30 40 50 60 65 70 75 80 85 90 95 100 } \
	    -entrystate disabled \
	    -labelwidth 13 \
	    -labelanchor w
    set ck [checkbutton $f5.cb1 \
	    -text "Display Shadow" \
	    -variable grafdata($gID,barshadow,1)]

    pack $ck -side top


    # FRAME #6
    set ok  [DefaultButton $f6.ok -text OK -done_var bar_atrib_color(done)]
    pack $ok -side top -expand 1 -padx 10 -pady 15

    tkwait variable bar_atrib_color(done)
    set ilist $i
    if $bar_atrib_color($gID,apply_to_all) {
	for {set n 1} {$n <= $grafdata($gID,N_point,1)} {incr n} {
	    append ilist "$n "
	}
    }
    foreach ii $ilist {
	if [info exists bar_atrib_color($gID,fill_cID)] {
	    set cID $bar_atrib_color($gID,fill_cID)
	    set grafsize($gID,bar${ii}_fill,$j) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	}

	if [info exists bar_atrib_color($gID,outline_cID)] {
	    set cID $bar_atrib_color($gID,outline_cID)
	    set grafsize($gID,bar${ii}_outline,$j) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	}
	
	if [info exists bar_atrib_color($gID,shadow_cID)] {
	    set cID $bar_atrib_color($gID,shadow_cID)
	    set grafsize($gID,bar${ii}_shadow,$j) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	}
	
	set grafsize($gID,bar${ii}_width,$j) $grafsize($gID,bar${i}_width,$j)
	if { $grafsize($gID,bar${ii}_width,$j) == 0 } {
	    set grafsize($gID,bar${ii}_outline,$j) "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
	}
	set grafsize($gID,bar${ii}_stipple,$j) \
		$grafsize($gID,bar${i}_stipple,$j)
    }

    if [info exists bar_atrib_color($gID,fill_cID)] {
	unset bar_atrib_color($gID,fill_cID)
    }
    if [info exists bar_atrib_color($gID,outline_cID)] {
	unset bar_atrib_color($gID,outline_cID)
    }
    if [info exists bar_atrib_color($gID,shadow_cID)] {
	unset bar_atrib_color($gID,shadow_cID)
    }

    set grafdata($gID,Xoffset,1) \
	    [expr -1.0 * double($bar_atrib_color($gID,xoffset)) / 100.0 + 1.0]

    destroy $t
}


proc SetBarAtribColor {gID j i what f {test {}}} {
    global grafsize multi_widget_list bar_atrib_color
    
    switch -exact -- $what {
	FILL {
	    if ![info exists bar_atrib_color($gID,fill_cID)] {
		set bar_atrib_color($gID,fill_cID) [xcModifyColorID]
	    }
	    set w [xcModifyColor $f "Fill Color:" \
		    $grafsize($gID,bar${i}_fill,$j) \
		    groove left left 100 100 70 5 15 \
		    $bar_atrib_color($gID,fill_cID)]
	}
	OUTLINE {
	    if ![info exists bar_atrib_color($gID,outline_cID)] {
		set bar_atrib_color($gID,outline_cID) [xcModifyColorID]
	    }
	    set w [xcModifyColor $f "Outline Color:" \
		    $grafsize($gID,bar${i}_outline,$j) \
		    groove left left 100 100 70 5 15\
		    $bar_atrib_color($gID,outline_cID)]
	}
	SHADOW {
	    if ![info exists bar_atrib_color($gID,shadow_cID)] {
		set bar_atrib_color($gID,shadow_cID) [xcModifyColorID]
	    }
	    set w [xcModifyColor $f "Shadow Color:" \
		    $grafsize($gID,bar${i}_shadow,$j) \
		    groove left left 100 100 70 5 15 \
		    $bar_atrib_color($gID,shadow_cID)]
	}
    }
    set multi_widget_list(post) $w
    if { $test == "test" } {
	tkwait visibility $w
	set wid [winfo width  $w]
	set hig [winfo height $w]
	set hlt [$w cget -highlightthickness]
	$f config -width  [expr $wid + 2 * $hlt + 6]
	$f config -height [expr $hig + 2 * $hlt + 6]
	pack propagate $f false
    }
}

