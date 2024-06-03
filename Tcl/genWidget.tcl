#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/genWidget.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# in this file are GENERAL WIDGET procedures like DIALOGS, ... !!!

proc RadioButtons { parent varname side args } {
    # side -- value for -side option
    set f [frame $parent.choices -relief groove -borderwidth 2]
    set b 0
    foreach item $args {
	radiobutton $f.$b -variable $varname \
		-text $item -value $item -anchor sw 
	pack $f.$b -side $side -fill both -padx 10 -pady 5
	incr b
    }
    pack $f -side top -ipadx 3 -ipady 3 -padx 5 -pady 5
}


proc RadioBut { parent labeltext varname lside rside ring \
	expand args } {
    # labeltext -- text to display in label
    # varname   -- name of variable
    # lside     -- value for -side option for label
    # rside     -- value for -side option for radiobutton
    # ring      -- wheather $f's relief is groove or not !!!
    # expand    -- value for -expand option

    set wlist {}

    if { $ring == 1 } {
	set f [frame $parent.f -relief groove -borderwidth 2]
    } else {
	set f [frame $parent.f  -borderwidth 0]
    }
    set f1 [frame $f.lbl -bd 0]
    set f2 [frame $f.f2 -bd 0]
    set lbl [label $f1.lbl -text $labeltext -anchor center]
    pack $f1 $f2 $lbl -side $lside -padx 0 -padx 0 \
	    -ipadx 0 -ipady 0 -expand $expand -fill both

    set wlist $lbl
    set b 0
    foreach item $args {
	radiobutton $f2.$b -variable $varname \
		-text $item -value $item -anchor sw 
	pack $f2.$b -side $rside -fill both -padx 0 -pady 0 \
		-ipadx 0 -ipady 0 -expand $expand
	lappend wlist $f2.$b
	incr b
    }
    pack $f -side $lside -ipadx 0 -ipady 0 -padx 2 -pady 2 \
	    -fill x -expand $expand
    return $wlist
}


proc RadioButCmd { parent labeltext varname cmd lside rside ring \
	expand padx args } {
    global radio_but_cmd_frame
    # labeltext -- text to display in label
    # varname   -- name of variable
    # cmd       -- command to execute--> a value of $item is passed to cmd
    #                                               ^^^^^
    # lside     -- value for -side option for label
    # rside     -- value for -side option for radiobutton
    # ring      -- wheather $f's relief is groove or not !!!
    # padx      -- padx value for groove frame
    # expand    -- value for -expand option

    set wlist {}

    if { $ring == 1 } {
	set f [frame [WidgetName $parent] -relief groove -borderwidth 2]
    } else {
	set f [frame [WidgetName $parent]  -borderwidth 0]
    }
    set radio_but_cmd_frame $f

    set f1 [frame $f.lbl -bd 0]
    set f2 [frame $f.f2 -bd 0]
    set lbl [label $f1.lbl -text $labeltext -anchor center]
    pack $f1 $f2 $lbl $lbl -side $lside -padx 0 -pady 0 \
	    -ipadx 0 -ipady 0 -expand $expand

    set wlist $lbl

    set b 0
    foreach item $args {
	radiobutton $f2.$b -variable $varname \
		-text $item -value $item \
		-anchor center \
		-command [list $cmd $item]  
	pack $f2.$b -side $rside -fill both -padx 0 -pady 0 \
		-ipadx 0 -ipady 0 -expand $expand
	lappend wlist $f2.$b	
	incr b
    }
    pack $f -side $lside -ipadx 0 -ipady 0 -padx $padx -pady 2 \
	    -fill x -expand $expand 

    return $wlist
}


proc RadioButVarCmd { parent labeltext varname cmd lside rside ring \
	expand args } {
    # labeltext -- text to display in label
    # varname   -- name of variable
    # cmd       -- command to execute--> a value of $varname is passed to cmd
    #                                               ^^^^^^^^
    # lside     -- value for -side option for label
    # rside     -- value for -side option for radiobutton
    # ring      -- wheather $f's relief is groove or not !!!
    # expand    -- value for -expand option

    if { $ring == 1 } {
	set f [frame $parent.f -relief groove -borderwidth 2]
    } else {
	set f [frame $parent.f  -borderwidth 0]
    }
    set f1 [frame $f.lbl -bd 0]
    set f2 [frame $f.f2 -bd 0]
    set lbl [label $f1.lbl -text $labeltext -anchor center]
    pack $f1 $f2 $lbl $lbl -side $lside -padx 0 -padx 0 \
	    -ipadx 0 -ipady 0 
    set b 0
    foreach item $args {
	radiobutton $f2.$b -variable $varname \
		-text $item -value $item \
		-command [list $cmd $varname] -anchor sw 
	pack $f2.$b -side $rside -fill both -padx 0 -pady 0 \
		-ipadx 0 -ipady 0
	incr b
    }
    pack $f -side $lside -ipadx 0 -ipady 0 -padx 2 -pady 2 \
	    -fill x -expand $expand
}


proc CheckButVarCmd { parent varname cmd side ring expand args } {
    # varname   -- name of variable
    # cmd       -- command to execute--> a value of $varname is passed to cmd
    #                                               ^^^^^^^^
    # side      -- value for -side option for radiobutton
    # ring      -- wheather $f's relief is groove or not !!!
    # expand    -- value for -expand option

    if { $ring == 1 } {
	set f [frame $parent.f -relief groove -borderwidth 2]
    } else {
	set f [frame $parent.f  -borderwidth 0]
    }
    set b 0
    foreach item $args {
	checkbutton $f.$b -variable $varname \
		-onvalue On -offvalue Off -text $item \
		-command [list $cmd $varname] -anchor sw 
	pack $f.$b -side $side -fill both -padx 0 -pady 0 \
		-ipadx 0 -ipady 0
	incr b
    }
    pack $f -side $side -ipadx 0 -ipady 0 -padx 2 -pady 2 \
	    -fill x -expand $expand
    return $f
}


proc CheckButtons { parent args } {
    set f [frame $parent.booleans -borderwidth 5]
    set b 0
    foreach item $args {
	checkbutton $f.$b -text $item -variable $item
	pack $f.$b -side left
	incr b
    }
    pack $f -side top
}


proc CheckVarButtons { parent 
		       labellist 
		       varlist 
		       side 
		       {onvalue 1} 
		       {offvalue 0}
		   } {
    set f [frame $parent.booleans -borderwidth 5]
    set b 0
    foreach label $labellist var $varlist {
	frame $f.$b
	checkbutton $f.$b.cb -text $label -variable $var \
	    -onvalue $onvalue -offvalue $offvalue
	pack $f.$b -side $side -fill x -expand 1
	pack $f.$b.cb -side left
	incr b
    }
    pack $f -side top
}


proc Entries { w lablist entrylist width {expand {1}} {side {left}} \
	{args {}}} {
    # w - widget
    # lablist  - list of labels
    # entrylist - list of entries variables
    # width - width of entry
    # expand - parameter for -expand option

    frame $w.frame 
    pack $w.frame -expand $expand
    
    set n 1
    foreach ebl $entrylist {
	set m [expr $n - 1] 
	label $w.frame.lab$n -text [lindex  $lablist $m] 
	entry $w.frame.entry$n -relief sunken -width $width \
		-textvariable $ebl
	eval {pack $w.frame.lab$n $w.frame.entry$n -side $side \
		-padx 5 -pady 5 -anchor w} $args
	incr n 
    }
    return $w.frame.entry1
}


proc OneEntries { w lablist entrylist labelwidth ewidth {pady 5} {args {}}} {
    # w - widget
    # lablist  - list of labels
    # entrylist - list of entries variables
    # labelwidth - width of label
    # ewidth - width of entry
    
    set n 1
    foreach ebl $entrylist {
	frame $w.frame$n 
	if { $args != {} } {
	    eval {pack $w.frame$n} $args
	} else {
	    pack $w.frame$n
	}
	set m [expr $n - 1] 
	label $w.frame${n}.lab$n -text [lindex  $lablist $m] \
		-width $labelwidth \
		-anchor w
	entry $w.frame${n}.entry$n -relief sunken \
		-textvariable $ebl \
		-width $ewidth

	lappend foclist $w.frame${n}.entry$n

	pack $w.frame${n}.lab$n -side left -padx 5 -pady $pady
	pack $w.frame${n}.entry$n -side right -padx 5 -pady $pady \
		-fill x -expand 1
	incr n 
    }

    return $foclist
}


proc dialog {win title text bitmap default args} {
    global button
    
    #win - name of top level window
    #title - title of toplevel window
    #text/bitmap - text/bitmap to be dispayed in the Dialog
    #default - index of default button; -1 if none
    #args - strings to be displayed in  the buttons

    # from where do we came here, that we'll be able to set 
    # the grab back where it was before

    set oldgrab [grab current]
    # 1. create TOP_LEVEL & divdide into TOP & BOTTOM
    
    toplevel $win -class Dialog
    wm title $win $title
    wm iconname $win Dialog
    frame $win.top -relief raised -bd 1
    pack $win.top -side top -fill both
    frame $win.bot -relief raised -bd 1
    pack $win.bot -side bottom -fill both

    xcPlace . $win 200 200
    # 2. fill TOP with bitmap & message    
    message $win.top.msg -text $text -aspect 500
    set font [ModifyFont fixed $win.top.msg \
	    -family helvetica -weight bold -size 12]
    $win.top.msg config -font $font
    pack $win.top.msg -side right -expand 1 -fill both \
	    -padx 5m -pady 5m
    if {$bitmap != "" } {
	label $win.top.bitmap -bitmap $bitmap
	pack $win.top.bitmap -side left -padx 5m -pady 5m
    }

    # 3. create a row of buttons at the BOTTOM

    set i 0
    foreach but $args {	
	if {$i == $default} {
	    frame $win.bot.default -relief sunken -bd 2
	    button $win.bot.default.button$i -text $but \
		    -command "set button $i"
	    pack $win.bot.default -side left -expand 1 \
		    -padx 5m -pady 2m
	    pack $win.bot.default.button$i -side left \
		    -padx 2m -pady 2m -ipadx 2m -ipady 1m
	    focus $win.bot.default.button$i
	} else {
	    frame $win.bot.rest -bd 10
	    button $win.bot.rest.button$i -text $but \
		    -command "set button $i"
            pack $win.bot.rest -side left -expand 1 \
		    -padx 5m -pady 2m
	    pack $win.bot.rest.button$i -side left -expand 1 \
		    -padx 2m -pady 2m -ipadx 2m -ipady 1m
	}
	incr i
    }

    # 4. Set up binding for <Return>

    if {$default > 0} {
	bind $win.bot.default.button$default <Return> \
		"$win.bot.default.button$default flash; \
		set button $default"
    }


    # set a grab
    
    tkwait visibility $win
    catch { grab $win }
    
    # 5. Wait for the user to respond, then release the grab
    # and return the index of the selected button.

    tkwait variable button
    destroy $win
    catch { grab release $win }

    # set grab to "old one"
    if { $oldgrab != {} } {
	catch { grab $oldgrab }
    }
    return $button
}


proc xcToplevel {w title iconname {master {.}} {x {0}} {y {0}} {transient 1}} {
    # w............name of toplevel
    # title........title of toplevel
    # iconname
    # master.......name of widow that will be used to place toplevel
    # x,y..........where to place toplevel
    
    #if { [winfo exist $w] } {
    #	xcDebug -stderr "\n\n\n DEBUG> toplevel \"$w\" already exist!!!!\n\
    #		ERROR: please report to autor: Tone.Kokalj@ijs.si\n"
    #	return 
    #} 

    
    if { [winfo exist $w] } {
	# toplevel already exists; return from the calling proc
	return -code return
    } 

    toplevel $w
    if { $master != "" } {
	xcPlace $master $w $x $y
	raise $w
    }
    wm title $w $title
    wm iconname $w $iconname
    
    if { $transient } { 
	wm transient $w [winfo toplevel [winfo parent $w]] 
    }
    return $w
}


# make text widget with xscrollball & yscrollbar and insert text
proc DispText {f text w h {update 0}} {
    # f...      window (YET TO BE CREATED)
    # text...   text to be displayed
    # w...      width of text widget
    # h...      height of text widget
    # update    if $f elready exists && update=1 -> just update the text
    #
    # PROC RETURNS name of text widget or 0 if it fails!!!!!!!

    # frame $f may already exists
    if { ![winfo exists $f] } {
	xcDebug -debug "#1"
	frame $f
	pack $f -side top -expand true -fill both
	set fb [frame $f.bottom	]
	set font [SetFont text -family courier -size 14]
	set t [text $f.t -setgrid true -wrap none -width $w -height $h \
		-font $font \
		-yscrollcommand "$f.sy set" -xscrollcommand "$fb.sx set"]
	puts stderr "TEXT-WIDGET: $t"

	scrollbar $f.sy -orient vert -command "$f.t yview"
	scrollbar $fb.sx -orient hori -command "$f.t xview"
	xcDebug -debug "#2"
	#set tplw .[lindex [split $f .] 1]; # whatfore is that used ????
	# Create padding based on the scrollbar width and border
	set pad [expr [$f.sy cget -width] + 2 * \
		([$f.sy cget -bd] + \
		 [$f.sy cget -highlightthickness])]
	frame $fb.pad -width $pad -height $pad
	xcDebug -debug "#3"
	pack $fb -side bottom -fill x
	pack $fb.pad -side right
	pack $fb.sx -side bottom -fill x
	pack $f.sy -side right -fill y
	pack $f.t -side left -fill both -expand true
	xcDebug -debug "#4"
	$f.t insert end $text
	$f.t config -state disabled
	return $f.t
    } elseif $update {
	# just update text
	set dis 0
	if { [$f.t cget -state] == "disabled" } {
	    set dis 1
	    $f.t configure -state normal
	}
	$f.t delete 1.0 end
	$f.t insert 1.0 $text

	if $dis {
	    $f.t configure -state disabled
	}
	return $f.t
    }
    return 0
}


proc OneEntryToplevel {w title iconname text width varname vartype x y} {
    global done oneentry
    upvar $varname var 
    
    if ![info exist var] { set var {} }
    set oneentry $var
    update
    xcDebug "OneEntryToplevel:: $oneentry"
    xcToplevel $w $title $iconname . $x $y
    set f1 [frame $w.f1 -relief raised -bd 2]
    set f2 [frame $w.f2 -relief raised -bd 2]
    set l1 [label $f1.l1 -text $text]
    set e1 [entry $f1.e1 -relief sunken -width $width -textvariable oneentry]
    focus $e1
    set varlist [list "oneentry $vartype"]
    set foclist $e1
    puts stdout "varlist:: $varlist"
    puts stdout "foclist:: $foclist"
    set b1 [button $f2.ok -text "OK" \
	    -command [list OneEntryOK $varlist $foclist]]
    set b2 [button $f2.can -text "Cancel" \
	    -command [list CancelProc $w]]
    pack $f1 $f2 -side top -fill both -padx 0 -pady 0
    pack $l1 $e1 -side top -expand 1 -padx 10 -pady 5
    pack $b1 $b2 -side left -expand 1 -padx 5 -pady 5

    bind $e1 <Return> [list OneEntryOK $varlist $foclist]
    bind $b1 <Return> [list OneEntryOK $varlist $foclist]
    tkwait visibility $w
    # check if there is some window grabed
    set oldgrab [grab current]	
    catch { grab $w }
    
    tkwait variable done
    destroy $w
    if { $oldgrab != ""} {
	catch { grab $oldgrab }
    }
    set var $oneentry
    xcDebug "OneEntryToplevel:: $oneentry"
    return $varname
}


proc OneEntryOK {varlist foclist} {
    global err done

    check_var $varlist $foclist
    if $err {return}
    set done 1
}
     
    
# proc makes Scrolled Entries on a Canvas
proc ScrollEntries { parent nn label labellist arraylist arraytypelist \
	width globvar buttonlist cheight } {
    global varlist foclist

    puts stdout "GLOBVAR NAME:: $globvar"

    # nn ....        number of Entries
    # label ........ top label
    # labellist .... list of labels
    # arraylist .... list of array elements
    #             expamle:  set arraylist "LB, NA,"
    #             name of variables is completed as: 
    #                                            $globvar(${varitem},$i)
    # arraytypelist . type of variable in array
    # width ....      width of entries
    # globvar     name of global variable
    # buttonlist .... 0 -> button do not exists
    #                 1 -> "<text1> <command1> <args1>" -> 1 button exist
    #                 2 -> "<text1> <command1> <args1>" "<text2> <command2> <args2>" -> 2 buttons exists
    #                 n -> "list #1" "list #2" ... "list #n"
    # cheight ....    height of canvas (Entries are units of height)


    # frame where canvas&scrollbar will be!!
    set ft [frame $parent.ft -relief sunken -bd 2]
    pack $ft -side top -expand true -fill y
	
    set c [canvas $ft.canv -yscrollcommand [list $ft.yscroll set]]
    set scb [scrollbar $ft.yscroll -orient vertical -command [list $c yview]]
    pack $scb -side right -fill y
    pack $c -side left -fill both -expand true
	
    # create FRAME to hold every LABEL&ENTRY
    set f [frame $c.f -bd 0]
    $c create window 0 0 -anchor nw -window $f -tags frame
    set varlist ""
    set foclist ""
    for {set i 1} {$i <= $nn} {incr i 1} {	    
	frame $f.fr$i -relief groove -bd 2 
	pack $f.fr$i -padx 5 -pady 5 -expand 1
	label $f.fr${i}.label$i -text "$label $i" 
	pack $f.fr${i}.label$i -anchor w -padx 7 -pady 7
	frame $f.fr${i}.frm$i 
	pack $f.fr${i}.frm$i -side top -anchor center 
	# coplite the varlist
	set tmplist ""
	set n 0
	foreach item $arraylist {	    
	    set var ${globvar}(${item},${i})
	    append tmplist " $var "
	    puts stdout "TMPLIST:: $var"
	    # make a varlist for PROC CHECK_VAR
	    lappend varlist "$var [lindex $arraytypelist $n]"
	    incr n
	}
	Entries $f.fr${i}.frm$i $labellist $tmplist $width
	set nb [lindex $buttonlist 0]
	for {set j 1} {$j <= $nb} {incr j} {
	    set com [lindex $buttonlist $j]
	    puts stdout "COM::: [list $com $i]"
	    set b [button $f.fr${i}.frm$i.b$j -text [lindex $com 0] \
		    -command [concat [lrange $com 1 end] $i]]
	    pack $b -side right -before $f.fr${i}.frm$i.frame -padx 10 -pady 5
	}
	# make a foclist for PROC CHECK_VAR
	lappend foclist \
		$f.fr${i}.frm$i.frame.entry1 $f.fr${i}.frm$i.frame.entry2
    }
	
    puts stdout "FOCLIST: $foclist\n\n"
    puts stdout "VARLIST: $varlist"
    set child [lindex [pack slaves $f] 0]
    
    # set the focus to first entry that upper FOR-LOOP create
    focus $f.fr1.frm1.frame.entry1
    
    tkwait visibility $child
    set width [winfo width $f]
    set height [winfo height $f]
    if { $nn < $cheight } {
    	 $c config -width $width -height $height 
    } else {
    	 $c config -width $width -height [expr $height / $nn * $cheight] \
    		 -scrollregion "0 0 $width $height"
    }

    return [list $varlist $foclist]
}


#-----------------------------
# SCROLLEDLISTBOX2
#-----------------------------
proc ScrolledListbox2 { parent args } {
	frame $parent
        pack $parent -side left -fill both -expand true -padx 5 -pady 10
	# Create listbox attached to scrollbars, pass thru $args
	eval {listbox $parent.list \
		-yscrollcommand [list $parent.sy set] \
		-xscrollcommand [list $parent.sx set]} $args
	scrollbar $parent.sy -orient vertical \
		-command [list $parent.list yview]
	# Create extra frame to hold pad and horizontal scrollbar
	frame $parent.bottom
	scrollbar $parent.sx -orient horizontal \
		-command [list $parent.list xview]
	# Create padding based on the scrollbar width and border
	set pad [expr [$parent.sy cget -width] + 2* \
		([$parent.sy cget -bd] + \
		 [$parent.sy cget -highlightthickness])]
	frame $parent.pad -width $pad -height $pad
	# Arrange everything in the parent frame
	pack $parent.bottom -side bottom -fill x
	pack $parent.pad -in $parent.bottom -side right
	pack $parent.sx -in $parent.bottom -side bottom -fill x
	pack $parent.sy -side right -fill y
	pack $parent.list -side left -fill both -expand true
	return $parent.list
}
#------------------------------------
# END OF SCROLLEDLISTBOX2
#------------------------------------



proc xcMenuEntry {parent l_text e_width e_var m_list {args {}}} {    
    global system
    # parent  ... parent widget
    # l_text  ... text for label
    # e_width ... width of entry
    # e_var   ... entry's textvariable
    # m_list  ... list of menu's entries
    # args    ... additional argumets that must be processed;
    #             -entryXXXXX means XXXXX atribute for entry

    set l [label $parent.l -text $l_text]
    set e [entry $parent.e -textvariable $e_var -width $e_width]
    set mb [menubutton $parent.mb \
	    -bitmap "@$system(BMPDIR)/xcMenuEntry_down.xbm" \
	    -menu $parent.mb.menu -relief raised]

    set menu [menu $mb.menu -tearoff 0 -relief raised]
    foreach word $m_list {
	$menu add command -label $word -command [list set $e_var "$word"]
    }
    
    # take care of options in $args
    if {$args == {}} { return 1 }
    set i 0    
    foreach option $args {
	incr i
	if { $i%2 } {
	    set tag $option
	} else {
	    switch -- $tag {
		"-labelrelief" {$l configure -relief $option}
		"-labelwidth"  {$l configure -width  $option}
		"-labelanchor" {$l configure -anchor $option}
		"-labelfont"   {$l configure -font   $option}
		"-entryrelief" {$e configure -relief $option}
		"-entryfont"   {$e configure -font   $option}
		"-entrystate"  {$e configure -state  $option}
		"-menubuttonrelief" {$mb configure -relief $option}
		"-menurelief"  {$menu configure -relief $option}
		default { tk_dialog .mb_error Error \
			"ERROR: Bad xcMenuEntry configure option $tag" \
			error 0 OK }
	    }
	}
    }

    if { $i%2 } {
	tk_dialog .mb_error1 Error "ERROR: You called xcMenuEntry with an odd number of args !" \
		error 0 OK
	return 0
    }

    pack $l -side left -padx 5 -pady 5
    pack $e -side left -fill x -pady 5
    pack $mb -side left -ipadx 2 -ipady 2 -pady 5 -padx 5

    return $e
}

	    
proc FillEntries { w lablist entrylist l_width e_width \
	{f_side top} {side left} {args {}}} {
    global xcFonts fillEntries
    # w         - parent widget
    # lablist   - list of labels
    # entrylist - list of entries variables
    # l_width   - width of label
    # e_width   - width of entry
    # f_side    - how to pack frame that holds frame & entry
    # side      - how to pack label & entry
    # args      - configuring options 
    set e_rel sunken
    set e_sta normal
    set e_bg  [GetWidgetConfig entry -background]
    set e_font $xcFonts(normal_entry)
    set l_font $xcFonts(normal)
    set i 0
    xcDebug "FillEntries Args:: $args"
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
	    xcDebug "FillEntries Options:: $tag $option"
            switch -- $tag {
                "-e_relief" {set e_rel  $option}
                "-e_state"  {set e_sta  $option}
		"-e_bg"     {set e_bg   $option}
		"-e_font"   {set e_font $option}
		"-l_font"   {set l_font $option}
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
		"ERROR: You called FillEntries with an odd number of args !" \
		error 0 OK
	return 0    
    }
    
    set i 1
    for {} {1} {incr i} {
	if ![winfo exists $w.f$i] {
	    set f $w.f$i
	    break
	}
    }

    frame $f     
    pack $f -expand 1 -fill both -side $f_side

    ##############################
    set n 1
    if { $e_width == {} } {
	set e_option [list -relief $e_rel \
		-state $e_sta \
		-bg $e_bg \
		-font $e_font]
    } else {
	set e_option [list -relief $e_rel \
		-width $e_width \
		-state $e_sta \
		-bg $e_bg \
		-font $e_font]
    }
    set fillEntries ""
    foreach ebl $entrylist {
	set m  [expr $n - 1]
	set fn [frame $f.$n]
	label $fn.lab$n -text [lindex  $lablist $m] \
		-width $l_width -anchor w -font $l_font
	lappend fillEntries $fn.entry$n
	eval {entry $fn.entry$n -textvariable $ebl} $e_option
	pack $fn -side $f_side -expand 1 -fill both -padx 5 -pady 2
	eval {pack $fn.lab$n -side $side} 
	eval {pack $fn.entry$n -side $side -fill x -expand 1} 
	incr n 
    }
    return $f.1.entry1
}



proc DisplayUpdateWidget {title text} {
    set t [xcToplevel [WidgetName] $title $title . 200 100 1]
    set m [message $t.m \
	    -text $text \
	    -aspect 500 \
	    -justify center\
	    -relief raised -bd 2 \
            -background "#f88" ]
    pack $m -expand 1 -ipadx 20 -ipady 20 -padx 0 -pady 0
    update
    #update idletask
    return $t
}


proc DefaultButton {name {args {}}} {
    
    set frame [frame $name -relief sunken -bd 2]
    
    xcDebug "DefaultButton Args:: $args"
    # args      - configuring options 
    set i 0
    set text ""
    set command ""
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
	    xcDebug "DefaultButton Options:: $tag $option"
            switch -- $tag {
                "-text"            {set text          $option}
                "-command"         {set command       $option}
		"-done_var"        {set done_var      $option}
		default { 
		    tk_dialog [WidgetName] Error \
			    "ERROR: Bad DefaultButton configure option $tag" \
			    error 0 OK 
		    return 0
		}
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error \
	"ERROR: You called DefaultButton with an odd number of args !" \
	error 0 OK
	return 0    
    }
    
    if { $command != "" } {
	button $name.b -text $text -command [list eval $command]
    } else {
	button $name.b -text $text -command [list set $done_var 1]
    }
    pack $name.b -side left -padx 1m -pady 1m
    focus $name.b

    return $frame
}


###############################################################################
# imitate pretty well the tk checkbutton, i.e takes the same option +
# option -image is possible, but does or have -indicatoron option
#
proc xcCheckButton {w args} {
    global checkButton
    
    
    #
    # get "-command" option out of $args
    #
    set id [button $w]

    set checkButton($id,pressed)     0
    set checkButton($id,is_variable) 0

    if {$args == {}} { 
	set args "-command xcCheckButtonDummy"
    }

    #
    # set default on/off value
    #
    set checkButton($id,offvalue) 0
    set checkButton($id,onvalue)  1

    set com 0
    set i   0
    set arg $args
    foreach option $arg {
	incr i
	if { $i%2 } {
	    set tag $option
	} else {
	    set j [lsearch $args $tag]		    
	    switch -- $tag {
		"-command" {
		    set com 1
		    set args [lreplace $args $j [expr $j + 1]]
		    set command $option		    
		}		
		"-offvalue"    { 
		    set args [lreplace $args $j [expr $j + 1]]
		    set checkButton($id,offvalue)    $option 
		}
		"-onvalue"     { 
		    set args [lreplace $args $j [expr $j + 1]]
		    set checkButton($id,onvalue)     $option 
		}
		"-selectcolor" { 
		    set args [lreplace $args $j [expr $j + 1]]
		    set checkButton($id,selectcolor) $option 
		}
		"-variable"    { 
		    set args [lreplace $args $j [expr $j + 1]]
		    set checkButton($id,is_variable) 1
		    set checkButton($id,variable)    $option
		    xcDebug "var:: $checkButton($id,variable)"

		}
	    }
	}
    }
    if !$com {
	set command xcCheckButtonDummy
    }

    puts stdout "args:: $args\n"
    flush stdout

    #
    # now configure the xcCheckButton
    #
    eval {$id configure} $args
    # this should be also tried with binding
    $id configure -command [concat xcCheckButtonCom $id $command]

    set checkButton($id,normalcolor) [$id cget -bg]
    if ![info exists checkButton($id,selectcolor)] {
	set checkButton($id,selectcolor) $checkButton($id,normalcolor)
    }

    # 
    # set the correct state
    #
    if $checkButton($id,is_variable) {
	upvar #0 $checkButton($id,variable) varn
	if { $varn == $checkButton($id,onvalue) } {
	    set checkButton($id,pressed) 1
	    $id configure -relief sunken \
		    -bg $checkButton($w,selectcolor)
	}
    }
	
    return $id
}


proc xcCheckButtonCom {w args} {
    global checkButton

    if $checkButton($w,is_variable) { 
	upvar #0 $checkButton($w,variable) varn
    }

    puts stdout "Com:: $w $varn $args"

    if !$checkButton($w,pressed) {
	# button was pressed
	set checkButton($w,pressed) 1
	if $checkButton($w,is_variable) {
	    set varn $checkButton($w,onvalue)
	}
	$w configure \
		-relief sunken \
		-bg $checkButton($w,selectcolor)
    } else {
	# button was releassed
	set checkButton($w,pressed) 0
	if $checkButton($w,is_variable) {
	    set varn $checkButton($w,offvalue)
	}
	$w configure \
		-relief raised \
		-bg $checkButton($w,normalcolor)
    }
    
    eval $args
}


proc xcCheckButtonDummy {} {
    return 0
}


proc xcCheckButtonRow {parent n bitmaplist varlist comlist \
	{fside left} {cbside left}} {
    #
    # n      ... number of xcCheckButtons
    # fside  ... side of frame to pack 
    # cbside ... side of checkbuttons to pack
    set f [frame $parent.f]
    pack $f -side $fside -expand 1

    for {set i 0} {$i < $n} {incr i } {
	set bmp [lindex $bitmaplist $i]
	set var [lindex $varlist    $i]
	set com [lindex $comlist    $i]
	puts stdout "$i: $bmp, $var, $com"
	set cb($i) [xcCheckButton $f.cb$i \
		-bitmap $bmp \
		-highlightthickness 0 \
		-selectcolor "#ffffff" \
		-command $com \
		-variable $var]	
	pack $cb($i) -side $cbside
    }
}


# ------------------------------------------------------------------------
# xcModifyColor and related routines
# ------------------------------------------------------------------------

proc xcModifyColorID {} {
    global mody_col
    
    if ![info exists mody_col(ID)] {
	set mody_col(ID) 1
    } else {
	incr mody_col(ID)
    }
    return $mody_col(ID)
}

proc xcModifyColorGetID {} {
    global mody_col

    if ![info exists mody_col(ID)] {
	return 0
    } else {
	return $mody_col(ID)
    }
}

proc xcModifyColor {parent labeltext init_color \
	frame_relief frame_side scale_side width height \
	scale_length scale_width slider_length {cID {}}} {
    global mody_col

    if { $cID == {} } {
	set cID [xcModifyColorID]
    }

    set f [frame [WidgetName $parent] -relief $frame_relief -bd 2]
    set l [label $f.l -text $labeltext -anchor w]
    pack $f -side $frame_side -padx 3 -pady 3 -ipadx 0 -ipady 0 \
	    -fill both -expand 1

    if ![info exists mody_col($cID,red)] {
	set color  [rgb_h2f $init_color]
	xcDebug "color:: $color"
	set mody_col($cID,red)   [lindex $color 0]
	set mody_col($cID,green) [lindex $color 1]
	set mody_col($cID,blue)  [lindex $color 2]
    }

    set fr  [frame $f.1 -relief sunken -bd 2]
    set mody_col($cID,col) [frame $fr.col  -bd 0 -width $width -height $height]
    _xcModifyColorSet $cID

    set f2 [frame $f.f2 -relief flat]
    scale $f2.red -from 0 -to 1 \
	    -length $scale_length \
	    -variable mody_col($cID,red) \
	    -orient horizontal -label "Red:" \
	    -digits 4 -resolution 0.001 -showvalue true \
	    -width $scale_width \
	    -sliderlength $slider_length \
	    -highlightthickness 0 \
	    -command [list _xcModifyColorSet $cID]
    scale $f2.green -from 0 -to 1 \
	    -length $scale_length \
	    -variable mody_col($cID,green) \
	    -orient horizontal -label "Green:" \
	    -digits 4 -resolution 0.001 -showvalue true \
	    -width $scale_width \
	    -sliderlength $slider_length \
	    -highlightthickness 0 \
	    -command [list _xcModifyColorSet $cID]
    scale $f2.blue -from 0 -to 1 \
	    -length $scale_length \
	    -variable mody_col($cID,blue) \
	    -orient horizontal -label "Blue:" \
	    -digits 4 -resolution 0.001 -showvalue true \
	    -width $scale_width \
	    -sliderlength $slider_length \
	    -highlightthickness 0 \
	    -command [list _xcModifyColorSet $cID]

    pack $l -side top -fill x -expand 1 -padx 10
    pack $fr $f2 -side $scale_side \
	    -fill both -expand 1 -padx 10 -pady 10 -ipadx 0 -ipady 0
    pack $mody_col($cID,col) -side top -fill both -expand 1 -padx 0 -pady 0 
    pack $f2.red $f2.green $f2.blue -side top -fill both -expand 1 \
	    -ipadx 0 -ipady 1  -pady 0

    return $f
}

proc _xcModifyColorSet {cID {dummy {}}} {
    global mody_col
    
    set mody_col($cID,hxred)   [d2h [expr round($mody_col($cID,red)   * 255)]]
    set mody_col($cID,hxgreen) [d2h [expr round($mody_col($cID,green) * 255)]]
    set mody_col($cID,hxblue)  [d2h [expr round($mody_col($cID,blue)  * 255)]]
    $mody_col($cID,col) configure \
	    -bg "#$mody_col($cID,hxred)$mody_col($cID,hxgreen)$mody_col($cID,hxblue)"
}

proc xcModifyColorSet {cID format type color} {
    global mody_col

    # NOTE: type is RGB or RGBA and is dummy for D and F (for
    #       compatibility with xcModifyColorGet)

    switch -glob -- $format {
	D - d* {
	    # D or decimal
	    set fc [rgb_d2f $color]
	    set mody_col($cID,red)   [lindex $fc 0]
	    set mody_col($cID,green) [lindex $fc 1]
	    set mody_col($cID,blue)  [lindex $fc 2]
	    _xcModifyColorSet $cID
	}
	F - f* {
	    # F or float
	    set mody_col($cID,red)   [lindex $color 0]
	    set mody_col($cID,green) [lindex $color 1]
	    set mody_col($cID,blue)  [lindex $color 2]
	    _xcModifyColorSet $cID
	}
	H - h* {
	    # H or hexadecimal
	    if { [string toupper $type] == "RGBA" } {
		set len   [expr {3 * ([string length $color] / 4)}]
		set color [string range $color $len]
	    }
	    set fc [rgb_h2f $color]
	    set mody_col($cID,red)   [lindex $fc 0]
	    set mody_col($cID,green) [lindex $fc 1]
	    set mody_col($cID,blue)  [lindex $fc 2]
	    _xcModifyColorSet $cID
	}
	deafult {
	    ErrorDialog "wrong format $format, must be one of D, F, or H"
	    return
	}
    }
}

proc xcModifyColorGet {cID format type} {
    global mody_col

    switch -glob -- $format {
	D - d* {
	    # D or decimal
	    set color [rgb_f2d [list $mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue)]]
	    if { [string toupper $type] == "RGBA" } {
		append color " 255"
	    }
	}
	F - f* {
	    # F or float
	    set color [list $mody_col($cID,red) $mody_col($cID,green) $mody_col($cID,blue)]
	    if { [string toupper $type] == "RGBA" } {
		append color " 1.0"
	    }
	}
	H - h* {
	    set color [list $mody_col($cID,hxred) $mody_col($cID,hxgreen) $mody_col($cID,hxblue)]
	    if { [string toupper $type] == "RGBA" } {
		append color "ff"
	    }
	}
	default {
	    ErrorDialog "wrong format $format, must be one of D, F, or H"
	    return ""
	}
    }
    return $color
}
# ------------------------------------------------------------------------
# END:: xcModifyColor
# ------------------------------------------------------------------------



# xcUpdate is toplevel window with Cancel, Update & Close button
proc xcUpdateWindow {{args {}}} {
    # options:
    #          -name
    #          -title
    #          -cancelcom
    #          -closecom
    #          -updatecom
    #          -frameside
    #          -buttonside
    #          -canceltext
    #          -closetext
    #          -updatetext

    # defaults
    set title      "Color Scheme"
    set updatecom  xcDummyProc
    set closecom   xcDummyProc
    set cancelcom  xcDummyProc
    set frameside  top
    set buttonside left
    set updatetext Update
    set closetext  Close
    set canceltext Cancel
    set name       [WidgetName]
    # parse args:
    set i 0    
    foreach option $args {
	incr i
	if { $i%2 } {
	    set tag $option
	} else {
	    switch -- $tag {
		"-name"       {set name $option}
		"-title"      {set title $option}
		"-updatecom"  {set updatecom $option}
		"-closecom"   {set closecom $option}
		"-cancelcom"  {set cancelcom $option}
		"-frameside"  {set frameside $option}
		"-buttonside" {set buttonside $option}
		"-updatetext" {set updatetext $option}
		"-closetext"  {set closetext $option}
		"-canceltext" {set canceltext $option}
		default { tk_dialog .mb_error Error \
			"ERROR: Bad xcUpdateWindow configure option $tag" \
			error 0 OK }
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error "ERROR: You called xcUpdateWindow with an odd number of args !" \
		error 0 OK
	return 0
    }

    set t [xcToplevel $name $title [lrange $title 0 2] . 0 0 1]

    set f1 [frame $t.f1 -class RaisedFrame]
    set f2 [frame $t.f2 -class RaisedFrame]
    pack $f1 $f2 -side $frameside -fill both -expand 1

    set can [button $f2.can -text $canceltext -command [list eval $cancelcom]]
    set upd [button $f2.upd -text $updatetext -command [list eval $updatecom]]
    set clo [button $f2.clo -text $closetext -command [list eval $closecom]]

    pack $can $upd $clo -side $buttonside \
	    -padx 10 -pady 10 -expand 1

    return $f1
}

proc xcMenuButton {w {args {}}} {
    # options:
    #   -labeltext
    #   -labelwidth
    #   -textvariable
    #   -menu {menutext1 menucom1 ...}

    # defaults:
    set labeltext    {}
    set labelwidth   {}
    set textv        xcMisc(dummy)
    set menu         {{} xcDummyProc}
    set side         left

    # parse args:
    set i   0    
    set wid 0
    foreach option $args {
	incr i
	if { $i%2 } {
	    set tag $option
	} else {
	    switch -- $tag {
		"-labeltext"    {set labeltext $option}
		"-labelwidth"   {set labelwidth $option}
		"-textvariable" {set textv $option}
		"-side"         {set side  $option}
		"-menu"         {
		    set nm 0
		    set wid 0
		    foreach {t c} $option {
			set wi [string length $t]
			if { $wid < $wi } { set wid $wi }
			set text($nm) $t
			set com($nm)  $c
			incr nm
		    }
		    incr wid 2
		}
		default { tk_dialog .mb_error Error \
			"ERROR: Bad xcMenuButton configure option $tag" \
			error 0 OK }
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error "ERROR: You called xcMenuButton with an odd number of args !" \
		error 0 OK
	return 0
    }

    set f [frame [WidgetName $w]]
    label $f.l -text $labeltext -relief flat -bd 0 -anchor w
    if { $labelwidth != {} } {
	$f.l config -width $labelwidth
    }
    upvar $textv value
    if { [info exists value] } {
	set len [string length $value]
	if { $len > $wid } {
	    set wid $len
	}
    }
    menubutton $f.mb \
	    -width $wid \
	    -textvariable $textv \
	    -menu $f.mb.menu \
	    -indicatoron 1 \
	    -relief raised

    set m [menu $f.mb.menu -tearoff 0]
    for {set i 0} {$i < $nm} {incr i} {
	$m add command -label $text($i) -command [list eval $com($i)]
    }
    pack $f.l $f.mb -side $side -padx 1 -fill x -anchor w
    return $f
}


#
# xcTextImageButton --
#      Create button with text+image but does not pack it
# 
# Arguments:
#      w      name of TextImageButton (it must be packed by the user)
#      text   text of textimagebutton
#      image  image of textimagebutton
#      side   how text and image is packed
#      args   arguments to button command
#
# Results:
#      Returns the name of the textimagebutton 
proc xcTextImageButton {w image side args} {
    button $w -highlightthickness 0
    $w config -state disabled
    set b [eval {button $w.b} $args {-bd 0 -highlightthickness 0}]
    set l [label $w.l -image $image -bg "#00f" \
	    -anchor c -bd 1 -highlightthickness 0]

    #foreach a [list $lt $li] b [list $li $lt] {
    #	bind $a <Enter>           +[list $b config -state active]
    #	bind $a <Leave>           +[list $b config -state normal]
    #}
    pack $b $l -side $side -fill both
    return $w
}


#
# special xcTextImageButton for "Hide"
proc xcHideButton {w image side args} {
    global xcFonts

    set font [SetFont button -size $xcFonts(small_size) -weight bold]
    eval {xcTextImageButton $w $image $side} $args {-bg "#00f" -fg "#fff" \
	    -activebackground "#88f" -activeforeground "#fff" -font $font}
}


#
# display content of a file in a separate toplevel window
# with scroll-text and Close widgets
#
proc xcDisplayFileText {file {title {Displayed Text}} \
	{w .} {x 0} {y 0} {transient 0}} {
    global system prop dispC95out unmapWin

    set f_content [ReadFile $file]

    return [xcDisplayVarText $f_content $title $w $x $y $transient]
}

proc xcDisplayVarText {varText {title {Displayed Text}} \
	{w .} {x 0} {y 0} {transient 0}} {

    set t  [xcToplevel [WidgetName] $title $title $w $x $y $transient]
    DispText $t.f1 $varText 80 20
    set f2 [frame $t.f2 -relief flat]
    pack $t.f1 -side top -expand 1 -fill both -padx 3 -pady 3
    pack $f2   -side top -fill x -padx 3 -pady 3


    button $f2.close -text "Close" -command [list destroy $t]
    pack $f2.close -side top -expand 1 -padx 3 -pady 3
    return $t
}



proc XCRYSDEN_Logo {file} {
    global xcMisc system

    #
    eval destroy [winfo children .]
    bind . <Destroy>   {}
    bind . <Configure> {}
    #/

    wm resizable . 0 0
    label .xcrysden_logo -image kpath -relief sunken -bd 2
    pack .xcrysden_logo -padx 2m -pady 2m -fill both -expand 1
    wm geometry  . +30+30
    wm deiconify .
    wm iconbitmap . @$system(BMPDIR)/xcrysden.xbm
    wm title . "*** XCrySDen *** "

    update
}
