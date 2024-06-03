#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/multiWid.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc MultiWidget {t args} {
    global multi_widget_list

    set i  0
    set ic 0
    set ib 0
    set create_tplw 0
    foreach option $args {
	incr i
        # odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
            switch -- $tag {
		"-create_tplw" {set create_tplw $option}
		"-tplw_args"   {set tplw_args $option}
		"-testbutton"  {set testbutton [expr $option - 1]}
                "-width"       {set width $option}
                "-height"      {set height $option}
		"-b_height"    {set b_height $option}
		"-b_width"     {set b_width  $option}
		"-command" {
		    foreach item $option {
			set b_text($ic) [lindex $item 0]
			set b_com($ic)  [lindex $item 1]
			incr ic
		    }
		}
		"-bottom_button" {
		    foreach item $option {
			set botb_text($ib) [lindex $item 0]
			set botb_com($ib)  [lindex $item 1]
			incr ib
		    }
		}
		default { 
		    tk_dialog .mb_error Error \
			    "ERROR: Bad FillEntries configure option $tag" \
			    error 0 OK 
		    return 0
		}
		
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error \
		"ERROR: You called MultiWidget with an odd number of args !" \
		error 0 OK
	return 0    
    }
    
    ########################################
    # is it needed to create toplevel ?????
    if $create_tplw {
	set t [eval $tplw_args]
    }
    ########################################

    set multi_widget_list(post) {}
    set mas [frame $t.mas]
    set top [frame $mas.top -highlightthickness 0]
    set bot [frame $mas.bot -relief raised -bd 3 \
	    -highlightthickness 0 ]
    if [info exists width]  { $bot configure -width  $width }
    if [info exists height] { $bot configure -height $height}
   
    set multi_widget_list(buttons) {}
    for {set i 0} {$i < $ic} {incr i} {	
	set b($i) [button $top.b$i -text $b_text($i) \
		-relief raised -bd 1 \
		-highlightthickness 0 \
		-command [list PostMultiWidget $bot $top.b$i $b_com($i)]]
	lappend multi_widget_list(buttons) $b($i)
	if [info exists b_height] { $b($i) config -height $b_height }
	if [info exists b_width]  { $b($i) config -width  $b_width }	
	pack $b($i) -side left -expand 1 -fill both	
    }
    
    pack $mas -padx 10 -pady 10 -fill both
    pack $top $bot -side top -fill both -padx 0 -pady 0 \
	    -ipadx 0 -ipady 0 -expand 1

    # BOTTOM-BUTTONS
    # get toplevel window
    set tplw [winfo toplevel $t]
    for {set i 0} {$i < $ib} {incr i} {
	button $mas.b$i \
		-relief raised -bd 3 \
		-highlightthickness 0 \
		-text $botb_text($i) \
		-command [list eval $botb_com($i) $tplw]
	
	pack $mas.b$i -side left -expand 1 -fill both -ipady 2
    }

    # test geometry of largest window
    PostMultiWidget $bot $top.b$testbutton $b_com($testbutton) test
    # default Posted widget is first one
    PostMultiWidget $bot $top.b0 $b_com(0) 
    puts stdout $mas
}

proc PostMultiWidget {f b com_list {test {}}} {
    global multi_widget_list

    foreach but $multi_widget_list(buttons) {
	$but configure -bd 1
    }
    $b configure -bd 3

    xcDebug "PostMultiWidget:: $multi_widget_list(post)"
    foreach wid $multi_widget_list(post) {
	if [winfo exists $wid] { destroy $wid }
    }
    set multi_widget_list(post) {}

    if { $test == {} } {
	eval $com_list $f
    } else {
	eval $com_list $f $test
    }
}

