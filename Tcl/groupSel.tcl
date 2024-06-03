#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/groupSel.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc group_sel { w type groups } {
    global groupsel n_groupsel done but_press

    # this Procedure "select" and return a GROUP (point, line, plane, 
    #                                             space group)
    #
    # w  - name of toplevel widget
    # ^
    # type  - name for type of group (point, line, plane, space )
    #
    # groups  - list of groups
    # ^^^^^^
    # groupsel  - selected group
    # ^^^^^^^^

    # groupsel is global; proc group_sel wait unit $but_press is TRUE 
    # and then exit; gropusel is returned
    if ![info exists groupsel ] {set groupsel {}}
    set but_press 0
    xcToplevel $w "Select a $type group" "Selection" . 120 70 1
    
    # make a top_frame & bottom_frame
    frame $w.top -relief raised -bd 1
    pack $w.top -side top -fill both 
    frame $w.bot -relief raised -bd 1
    pack $w.bot -side bottom -fill both 
    
    # name of widget goes in the top frame
    label $w.top.lab -text "Select a $type group" 
    pack $w.top.lab -side top -expand 1 -fill both -padx 5m -pady 5m
    
    # in the top of BOTTOM frame will be an ENTRY for selected group;
    # bellow will be "selection listbox"
    label $w.top.label -text "Group:" -padx 0

    
    #------------------------------------------------------
    # take care of this ENTRY and groupsel textvariable
    #------------------------------------------------------
    entry $w.top.entry -relief sunken -textvariable groupsel
    focus $w.top.entry

    # make LEFT & RIGHT frame in bottom frame;
    # in LEFT goes LISTBOX;
    # in RIGHT go OK & CANCEL buttons
    frame $w.bot.left
    frame $w.bot.right
    frame $w.bot.right.ok
    set ok [button $w.bot.right.ok.ok -text OK -command [list group_sel_ok $groups]]
    set can [button $w.bot.right.can -text Cancel \
	    -command [list group_sel_cancel $w]]
   
    # now we'll create LISTBOX; it's ScrolledListbox2 from the book of
    # Brent B. Welch
    ScrolledListbox2 $w.bot.left.lb -width 20 -height 20 -setgrid true

    # now we'll pack what we create so far
    pack $w.top.label -side left -pady 10 -padx 5
    pack $w.top.entry -side left -fill x -expand true -pady 10 -padx 5
    pack $w.bot.left -side left -fill both
    pack $w.bot.right -side left -fill both    
    pack $w.bot.right.ok $can -side top -padx 10 -pady 5
    pack $ok -padx 4 -pady 4
    
    # now we'll create a bindings for ENTRY
    bind $w.top.entry <Return> "set done 1" 
    # now we'll create BINDINGs for this LISTBOX
    bind $w.bot.left.lb.list <ButtonPress-1> {GroupSelectStart %W %y}
    bind $w.bot.left.lb.list <ButtonRelease-1> \
	    [list GroupSelectEnd %W %y]

    # now we'll insert the GROUPS in the listbox;
    # GROUPS variable holds all groups 
    eval {$w.bot.left.lb.list insert 0} $groups 

 
    # wait unit $done is specefied
    tkwait variable but_press
    destroy $w
    return $groupsel
}
#-----------------------------
# END OF GROUP_SEL PROC
#-----------------------------


#------------------------------------
# here are auxiliary procs
#------------------------------------


proc GroupSelectStart { w y } {
    $w select anchor [$w nearest $y]
}

proc GroupSelectEnd { w y } {
    global groupsel

    # $w is the name of listbox widget who contains groups
    # $y is the vertical position of "selection"
    
    $w select anchor [$w nearest $y]
    # nline is index of selected line
    set nline [$w curselection]
    set groupsel [$w get $nline]
    # now we must get rid of seq_num: (ex. 20:   GROUP)
    #                        ^^^^^^^       ^^^
    # i=position/index of the first caracter og group which is 4 more than
    # index of :
    # j=lenght_of$groupsel
    set i [expr [string first : $groupsel] + 4]
    set j [string length $groupsel]
    # assing "pure" group to groupsel 
    set groupsel [string range $groupsel $i $j]
    }


proc group_sel_cancel { w } {
    global groupsel but_press
    set groupsel {}
    set but_press 1
}

proc group_sel_ok { groups } {
    global groupsel n_groupsel done but_press err
    # $groups is a list of groups

    # err is used by contrl_var variable
    set err 0

    # translate $groupsel to UPPER CASE
    set groupsel [string toupper $groupsel]

    # check if specified group is OK
    # $groups contanins "inpurities" (ex.{  5:   P 1 1 2        }), 
    # but $groupsel could also contains white_space inpurity

    set ok 0
    set n 1
    # pure the groupsel
    # NOTE: this REGEXP is extremly complex, but it looks that it works
    regexp {(([A-Z0-9] )|[A-Z0-9\/\-])+[A-Z0-9]} $groupsel groupsel

    foreach word $groups {
	# pure the $word
	set last [ string length $word ]
	set word [ string range $word 5 $last ]
	regexp {(([A-Z0-9] )|[A-Z0-9\/\-])+[A-Z0-9]} $word word	
	if { $groupsel == $word } {
	    set ok 1
	    set n_groupsel $n
	}
	incr n
    }
    
    # if "ok" is still is still 0, it maybe not a standard group; 
    # but that is just for crystals
    if { $ok == 0 } { 
	set ok [IsNotStandardGroup $groupsel] 
	set n_groupsel 999
    }

    if { $ok == 0} {
	# specefied group is WRONG
	dialog .group_sel_warning Warning \
		"Group \"$groupsel\" is false. You probably mistype the group.\
		Try again or select one from the list" warning 0 OK
	set err 1
    } else {
	# group is OK
	# we have selected; so $done will be set to 1
	set done 1	
	#this is used for GROUP_SEL PROC
	set but_press 1
    }
}



#############################################################################
# this proc is used when non-standard group is used that is not on 
# $space_group list; we check if it is valid and which parameters are needed
# for it!!!
proc IsNotStandardGroup {group} {
    global class inp system
    
    set par(1) 11.00111
    set par(2) 22.00222
    set par(3) 33.00333
    set par(4) 44.00444
    set par(5) 55.00555
    set par(6) 66.00666

    set cla    {}
    set cl(1)  A
    set cl(2)  B
    set cl(3)  C
    set cl(4)  ALFA
    set cl(5)  BETA
    set cl(6)  GAMMA
 
    set    input "EOF"
    append input "XCrySDen 1.0\n"
    append input "CRYSTAL\n"
    append input "1 0 0\n"
    append input "$group\n"
    append input "$par(1) $par(2) $par(3) $par(4) $par(5) $par(6)\n"
    append input "1\n"
    append input "1 0.0 0.0 0.0\n"
    append input "STOP\n"
    append input "XCrySDen 1.0\n"
    append input "EOF"
    set nlat 0
    set num  ""

    cd $system(SCRDIR)
    set output [RunAndGetC95Output $system(c95_integrals) {} $input]
    xcDebug "output: \n$output"
	
    # check if error occured

    set is_error 0
    if { [string match *ERROR* $output] } { 	
	set is_error 1
	if { ($system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14" ) && [string match "*STOP KEYWORD - EXECUTION STOPS*" $output] } {
	    # OK, we don't have error
	    set is_error 0
	}
    }
    
    if { $is_error } {
	return 0
    } else {
	set output [split $output \n]
	foreach line $output {
	    if [string match "*CRYSTAL FAMILY*" $line] {
		set inp(CRY_FAM)    [lrange $line 3 end]
	    }
	    if [string match "*CRYSTAL CLASS*" $line] {
		set inp(CRY_CLASS)  [lrange $line 6 end]
	    }
	    # LATTICE PARAMETERS
	    if { $nlat == 2 } {
		for {set i 1} {$i <= 6} {incr i} {
		    set ii [expr $i - 1]
		    set re($i) [lindex $line $ii]
		    for {set j 1} {$j <= 6} {incr j} {
			if { $re($i) == $par($j) } {
			    append cla "$cl($i) "			
			    # delete element that was found from "par" array
			    set par($j) ""
			}
		    }
		}
		incr nlat
	    }
	    if { $nlat == 1 } {
		incr nlat
	    }
	    if { [string match "*ICE PAR*" $line] && $nlat == 0 } {
		incr nlat
	    }
	}

	if { ! [info exists cla] } {
	    tk_messageBox -message "Ups! This is a bug in the program. Please
report it to: tone.kokalj@ijs.si with the detailed explanation what yoy were doing and if possible also Email the CRYSTAL input file you have been working on when the error occured !\n\nThe application will exit now." -type ok -icon info
	    exit_pr immediately
	}

	foreach item $cla {
	    if { $item == "A" || $item == "B" || $item == "C" } {
		append class1 "$item "
	    } else {
		append class2 "$item "
	    }
	}
	# get rid of last space-character in class1 & class2 
	set class1 [string trimright $class1 " "]
	if [info exists class2] { set class2 [string trimright $class2 " "] }
	set class ""
	lappend class $class1
	if [info exists class2] { lappend class $class2 }

	#set class [list $class]
	puts stdout "GROUP: $group CLASS: $class"
	return 1
    }
}


