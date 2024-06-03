#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/atomLabels.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


proc ModAtomLabels {} {
    global atomLabel

    #
    # make three-pages: 
    # page.1: for globalAtomLabel properties
    # page.2: for each-atom-label properties
    # page.3: advanced setting

    set t .atom_label
    if { [winfo exists $t] } {
	return
    }
    xcToplevel $t "Edit Atom labels & fonts" "Atom Labels"

    set nb [NoteBook $t.nb]
    pack $nb -expand 1 -fill both


    # --------------------------------------------------
    # page #.1: globalAtomLabel
    # --------------------------------------------------
    $nb insert 0 globalAtomLabel -text "Global Atom-Labels Font"
    set page [$nb getframe globalAtomLabel]

    ModAtomLabels:_font_and_color $t $page globalFont


    # --------------------------------------------------
    # page #.2: atomLabel
    # --------------------------------------------------
    $nb insert 1 atomLabel -text "Edit Custom Atom-Labels and Fonts"
    set page [$nb getframe atomLabel]

    set atom_frame [frame $page.atom  -relief groove -bd 2]
    pack $atom_frame -side top -expand 1 -fill both -padx 3 -pady 3 -ipadx 3 -ipady 3
    
    foreach i {1 2} text {	"Atom ID:" "ID's Atom Label:" } var {
	atomLabel(atomFont.id) atomLabel(atomFont.label)
    } ew { 5 25 } {
	set l [label $atom_frame.l$i -text $text]
	set e [entry $atom_frame.e$i -textvariable $var -width $ew]
	pack $l $e -side left -padx 1 -pady 3
	if { $i == 1 } {
	    bind $e <Enter> [list ModAtomLabels:atomFontUpdate $atomLabel(atomFont.id)]
	    set atomLabel(atomFont.entry) $e
	    set select [button $atom_frame.selectbutton \
			    -text "Select Atom" -command ModAtomLabels:selectAtom]
	    pack $select -side left -padx 0 -pady 5
	}
    }
    ModAtomLabels:_font_and_color $t $page atomFont


    # --------------------------------------------------
    # page #.3: advanced
    # --------------------------------------------------
    $nb insert 2 advanced -text "Advanced"
    set page [$nb getframe advanced]

    set check_frame    [frame $page.check  -relief groove -bd 2]
    set entry_frame(1) [frame $page.entry1 -relief groove -bd 2]
    set entry_frame(2) [frame $page.entry2 -relief groove -bd 2]
    set button_frame   [frame $page.button]
    pack $check_frame $entry_frame(1) $entry_frame(2) -side top \
	-fill x -padx 3 -pady 3 -ipady 2
    pack $button_frame -side bottom -fill x -padx 3 -pady 3 -ipady 10

    # check-buttons ...
    foreach i {cb1 cb2} text {
	"do not display all default atomic symbols"
	"do not display all custom atomic labels"
    } onvalue {0 0} offvalue {1 1} var {
	atomLabel(globalFont.do_display)
	atomLabel(atomFont.do_display)
    } command {
	{ModAtomLabels:advancedCheckButton default}
	{ModAtomLabels:advancedCheckButton custom}	
    } {
	set cb($i) [checkbutton $check_frame.$i -text $text \
			-onvalue $onvalue -offvalue $offvalue \
			-variable $var -command $command -anchor w]
	pack $cb($i) -side top -padx 3 -pady 3 -fill x
    }
    set b [button $check_frame.b1 -text "Clear all custom atomic labels" \
	       -command ModAtomLabels:advancedClearCustomLabels]
    pack $b -side left -padx 3 -pady 3

    # entries ..    
    foreach i {1 2} text {
	"do not display label for atoms (enter atom ID's):"
	"do not display labels for the following type of\natoms (enter atomic symbols):"
    } textvariable {
	atomLabel(atomIDs.do_not_display)
	atomLabel(atomTypes.do_not_display)
    } b1text {
	"Select Atom" "Select Atomic Type"
    } b2text {
	"De-select Atom" "De-select Atomic Type"
    } b1com {
	{ModAtomLabels:advancedSelectAtom select}
	{ModAtomLabels:advancedSelectType select}
    } b2com {
	{ModAtomLabels:advancedSelectAtom deselect}
	{ModAtomLabels:advancedSelectType deselect}
    } {
	set f1($i) [frame $entry_frame($i).1$i]
	set f2($i) [frame $entry_frame($i).2$i]
	set la($i) [label $f1($i).l$i -text $text -anchor w]
	set en($i) [entry $f1($i).e$i -textvariable $textvariable -width 40]	
	set b1($i) [button $f2($i).b1$i -text $b1text -width 14 -command $b1com]
	set b2($i) [button $f2($i).b2$i -text $b2text -width 21 -command $b2com]

	pack $f1($i) $f2($i) -side top -expand 1 -fill both -padx 2 -pady 2
	pack $la($i) -side left -padx 1
	pack $en($i) -side top  -padx 1 -fill x -expand 1
	pack $b1($i) $b2($i) -side left -padx 1
    }

    # buttons ...
    set cancel [button $button_frame.can -text "Cancel" \
		    -command [list CancelProc $t]]
    set close  [button $button_frame.close -text "Close" \
		    -command [list ModAtomLabels:advancedCloseUpdate $t close]]
    set update [button $button_frame.update -text "Update" \
		    -command [list ModAtomLabels:advancedCloseUpdate $t update]]
    pack $cancel $close $update -side left -expand 1 -padx 7 -pady 7 -ipadx 2 -ipady 2


    # --------------------------------------------------
    # raise the globalAtomLabel page
    # --------------------------------------------------
    $nb raise globalAtomLabel
    return $t
}

proc ModAtomLabels:_font_and_color {t page whichFont} {
    global atomLabel

    if { $whichFont != "atomFont" && $whichFont != "globalFont" } {
	ErrorDialog "wrong whichFont $whichFont for ModAtomLabels:_font_and_color"
	return
    }
    
    set font_frame   [frame $page.font  -relief groove -bd 2]
    set color_frame  [frame $page.color -relief groove -bd 2]
    set button_frame [frame $page.button]
    pack $font_frame $color_frame $button_frame -side top -expand 1 \
	-fill both -padx 3 -pady 3 -ipady 2

    # font ...
    
    FillEntries $font_frame {"Font:"} atomLabel($whichFont) 5 60
    set b [button $font_frame.browse -textvariable atomLabel(fontBrowser) \
	       -width 19 \
	       -command [list ModAtomLabels:fontBrowser $whichFont]]
    
    # TEMPORARY
    global tcl_platform
    if { $tcl_platform(platform) == "windows" } {
	# disable font-browser, as it does not work under windows
	$b config -state disabled
    }
    #/

    foreach i {s c} value {
	"Simple Font Browser" "Font Browser"
    } text {
	"use basic fonts only" "use all available fonts"
    } {
	set rb($i) [radiobutton $font_frame.rb$i \
			-text $text \
			-value $value \
			-variable atomLabel(fontBrowser)]
    }
    pack $b $rb(s) $rb(c) -side left -padx 5 -pady 5    
    $rb(c) configure -state disabled; # disable "use all available fonts" as it does not work Togl-font problem

    # font-colors ...

    #proc xcModifyColor {parent labeltext init_color \
    #			    frame_relief frame_side scale_side width height \
    #			    scale_length scale_width slider_length {cID {}}} 
    set f1 [frame $color_frame.1]
    set f2 [frame $color_frame.2]
    pack $f1 $f2 -side left -expand 1 -padx 5 -pady 5

    if { ! [info exists atomLabel($whichFont.brightID)] } {
	set atomLabel($whichFont.brightID) [xcModifyColorID]
	set atomLabel($whichFont.darkID)   [xcModifyColorID]    
    }
    set brightColor [rgb_f2h $atomLabel($whichFont.brightColor)]
    set darkColor   [rgb_f2h $atomLabel($whichFont.darkColor)]

    xcModifyColor $f1 "Set Bright Font Color:" \
	$brightColor groove left left 100 100 70 5 20 \
	$atomLabel($whichFont.brightID)
    
    xcModifyColor $f2 "Set Dark Font Color:" \
	$darkColor groove left left 100 100 70 5 20 \
	$atomLabel($whichFont.darkID)
    #if { $whichFont == "atomFont" } {
    #	if { [trace vinfo atomLabel(atomFont.brightColor)] != "" } {
    #	    trace variable atomLabel(atomFont.brightColor) w ModAtomLabels:trace
    #	}
    #	if { [trace vinfo atomLabel(atomFont.darkColor)] != "" } {
    #	    trace variable atomLabel(atomFont.darkColor) w ModAtomLabels:trace
    #	}
    #}

    # buttons ...
    set cancel [button $button_frame.can -text "Cancel" \
		    -command [list CancelProc $t]]
    set close  [button $button_frame.close -text "Close" \
		    -command [list ModAtomLabels:CloseUpdate $t $whichFont close]]
    set update [button $button_frame.update -text "Update" \
		    -command [list ModAtomLabels:CloseUpdate $t $whichFont update]]
    pack $cancel $close $update -side left -expand 1 -padx 7 -pady 7 -ipadx 2 -ipady 2
}


proc ModAtomLabels:selectAtom {} {
    global select atomLabel

    if { ! [info exists atomLabel(selecAtomWindow)] } {
	set atomLabel(selecAtomWindow) .atom_label_selectiton
    }
    if { [winfo exists $atomLabel(selecAtomWindow)] } {
	return
    }

    set select(done) 0
    PreSel $atomLabel(selecAtomWindow) .mesa "Select an Atom" \
	"For selecting an atom click on particular atom" AtomSel 1
    tkwait window $atomLabel(selecAtomWindow)
    if { !$select(done) } {
	return
    }

    set atomLabel(atomFont.id)    $select(Sqn1)
    set atomLabel(atomFont.label) $select(Sym1)
    ModAtomLabels:atomFontUpdate  $select(Sqn1)
}

proc ModAtomLabels:atomFontUpdate {id} {
    global atomLabel mody

    if { ! [string is integer $id] || $id == "" } {
	return
    }
    # query the label of atom $id
    
    set label [xc_getvalue $mody(GET_ATOMLABEL_LABEL) $id]
    if { $label != "" } {
	set atomLabel(atomFont.label) $label
    }
    
    # query the font-colors atom $id
    set atomLabel(atomFont.brightColor) [xc_getvalue $mody(GET_ATOMLABEL_BRIGHTCOLOR) $id]
    set atomLabel(atomFont.darkColor)   [xc_getvalue $mody(GET_ATOMLABEL_DARKCOLOR) $id]

    xcModifyColorSet $atomLabel(atomFont.brightID) float rgb $atomLabel(atomFont.brightColor)
    xcModifyColorSet $atomLabel(atomFont.darkID)   float rgb $atomLabel(atomFont.darkColor)
}


proc ModAtomLabels:fontBrowser {whichFont} {
    global atomLabel

    set t .atomlabelfontbrowser 
    if { [winfo exists $t] } {
	return
    }
    if { ! [info exists atomLabel(lastfont) ] } {
	#set atomLabel(lastfont) ""
        #set atomLabel(lastfont) "-adobe-helvetica-medium-r-normal--18-180-75-75-p-98-iso8859-1"
        set atomLabel(lastfont) [font create \
                                     -family Helvetica \
                                     -size   18 \
                                     -weight bold \
                                     -slant roman]         
    }
    set allFonts ""
    if { $atomLabel(fontBrowser) == "Font Browser" } {
	set allFonts 1
    }
    set font [fontToplevelWidget .atomlabelfontbrowser "Sample Text" $atomLabel(lastfont) $allFonts]
    #puts stderr "*** ModAtomLabels:fontBrowser :   font = $font"

    if { $font != "" } {
	set atomLabel(lastfont) $font
	set atomLabel($whichFont) [xcTkFontName2XLFD $font]
    }
}


proc ModAtomLabels:CloseUpdate {t type action} {
    global atomLabel mody_col

    if { $type != "atomFont" && $type != "globalFont" } {
	ErrorDialog "wrong type $type for ModAtomLabels:CloseUpdate"
	return
    }
    
    set font $atomLabel($type)

    set bright_color [xcModifyColorGet $atomLabel($type.brightID) float RGB]
    set dark_color   [xcModifyColorGet $atomLabel($type.darkID)   float RGB]
    
    if { $type == "globalFont" } {
	set atomLabel(globalFont.brightColor) $bright_color
	set atomLabel(globalFont.darkColor)  $dark_color
        #puts stderr "*** ModAtomLabels:CloseUpdate: xc_setfont .mesa $font $bright_color $dark_color"
	xc_setfont .mesa $font $bright_color $dark_color

    } elseif { $type == "atomFont" } {
	# here
	set id $atomLabel($type.id)
	set atomLabel($id.$type.label)      $atomLabel($type.label)
	set atomLabel($id.$type.font)       $font
	set atomLabel($id.$type.brightColor) $bright_color
	set atomLabel($id.$type.darkColor)  $dark_color

	xc_setatomlabel .mesa $id $atomLabel($type.label) $font $bright_color $dark_color
    }
    
    if { $action == "close" } {
	CancelProc $t
    }
}


# ------------------------------------------------------------------------
# ModAtomLabels:advanced* procs
# ------------------------------------------------------------------------

proc ModAtomLabels:advancedCheckButton {mode} {
    global atomLabel mody

    if { $mode == "default" } {
	#
	# do not display all global atom labels
	#
	if { $atomLabel(globalFont.do_display) } {
	    xc_newvalue .mesa $mody(SET_GLOBALATOMLABEL_DO_DISPLAY) 1
	} else {
	    xc_newvalue .mesa $mody(SET_GLOBALATOMLABEL_DO_DISPLAY) 0	    
	}
    } else {
	#
	# do not display all custom atom labels
	#
	set customLabelIDs [xc_getvalue $mody(GET_ATOMLABEL_ALL_ID)]
	if { $atomLabel(atomFont.do_display) } {
	    set display 1
	} else {
	    set display 0
	}
	foreach id $customLabelIDs {
	    # here checking for do_display 
	    xc_newvalue .mesa $mody(SET_ATOMLABEL_DO_DISPLAY) $id $display
	}
    }
}

proc ModAtomLabels:advancedSelectAtom  {action} {
    global atomLabel select

    set select(done) 0
    set selw [WidgetName]
    if { $action == "select" } {
	PreSel [WidgetName] .mesa "Select an Atom" \
	    "For selecting an atom click on particular atom" AtomSel 1
    } else {
	# action == deselect
	PreSel [WidgetName] .mesa "(De-)Select an Atom" \
	    "For (de-)selecting an atom click on particular atom" AtomSel 1
    }

    tkwait window $selw
    if { ! $select(done) } {
	return
    }

    if { $action == "select" } {
	lappend atomLabel(atomIDs.do_not_display) $select(Sqn1)
    } else {
	# action == deselect
	set ind [lsearch $atomLabel(atomIDs.do_not_display) $select(Sqn1)]
	if { $ind > -1 } {
	    set atomLabel(atomIDs.do_not_display) \
			      [lreplace $atomLabel(atomIDs.do_not_display) $ind $ind]
	}
    }
}

proc ModAtomLabels:advancedSelectType  {action} {
    global atomLabel ptable

    ptable . -command ptableSelectElement
    if { $ptable(result) == "" } { 
	return 
    }
    
    if { $action == "select" } {
	lappend atomLabel(atomTypes.do_not_display) $ptable(result)	
    } else {	
	# action == deselect	
	set ind [lsearch $atomLabel(atomTypes.do_not_display) $ptable(result)]
	if { $ind > -1 } {
	    set atomLabel(atomTypes.do_not_display) \
			      [lreplace $atomLabel(atomTypes.do_not_display) $ind $ind]
	}
    }
}

proc ModAtomLabels:advancedCloseUpdate {t action} {
    global atomLabel mody

    # build an do_not_display ID-list, we have: 
    #
    #    atomLabel(atomIDs.do_not_display)
    #    atomLabel(atomTypes.do_not_display)

    if { ! [info exists atomLabel(atomIDs.do_not_display)] } {
	set atomLabel(atomIDs.do_not_display) ""
    }
    if { ! [info exists atomLabel(atomTypes.do_not_display)] } {
	set atomLabel(atomTypes.do_not_display) ""
    }

    # IDs-part from atomIDs ...
    set IDs $atomLabel(atomIDs.do_not_display) 

    # IDs-part from atomTypes ...
    set natList ""
    foreach type $atomLabel(atomTypes.do_not_display) {
	set _nat [Aname2Nat $type]
	if { ! [string is integer $_nat] } {
	    ErrorDialog "$type is not atomic symbol"
	} else {		
	    lappend natList $_nat
	}
    }
    if { $natList != "" } {
	set natoms [xc_getvalue $mody(GET_NATOMS)]
	for {set i 1} {$i <= $natoms} {incr i} {
	    set ind [lsearch $natList [xc_getvalue $mody(GET_NAT) $i]]
	    if { $ind > -1 } {
		lappend IDs $i
	    }
	}
    }

    if { ! [info exists atomLabel(atomIDs.do_no_display.previous)] } {
	set atomLabel(atomIDs.do_no_display.previous) ""
    }
    
    # make visible all the labels that were disabled in the last call
    # to ModAtomLabels:advancedCloseUpdate    
    
    foreach id $atomLabel(atomIDs.do_no_display.previous) {
	xc_newvalue .mesa $mody(SET_DO_NOT_DISPLAY_ATOMLABEL) $id 0
    }
    
    # disable all curently specified labels

    foreach id $IDs {			    
	xc_newvalue .mesa $mody(SET_DO_NOT_DISPLAY_ATOMLABEL) $id 1
    }
    
    # current IDs list become previous-list

    set atomLabel(atomIDs.do_no_display.previous) $IDs

    if { $action == "close" } {
	CancelProc $t
    }
}


proc ModAtomLabels:advancedClearCustomLabels {} {
    set button [tk_messageBox -message "All custom atomic labels are about to be cleared. All your custom atomic-label editings will be lost.\n\nDo you really want to do it?" \
		    -type yesno -icon question]
    if { $button == "yes" } {
	xc_clearatomlabel .mesa all
    }
}
