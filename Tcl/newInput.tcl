#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/newInput.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#-------------------------------------------------------------
# VARIABLES used for INPUT
#-------------------------------------------------------------
# job_title ... title of job
# A     ....... "A" lattice vector
# B     ....... "B" lattice vector
# C     ....... "C" lattice vector
# ALFA  ....... "BC" angle(degrees) - between "B" & "C"
# BETA  ....... "AC" angle(degrees) - between "A" & "C"
# GAMMA    .... "AB" angle(degrees) - between "A" & "B"
# NATR  ....... number of non-equivalent atoms
# Nat   ....... atomic number (an array)
# X,Y,Z ....... atomic coordinates (an array)
# IFSO  ....... settings for the origin for the crystal reference frame
# IX,IY,IZ .... coordinates of non-standard shift of origin 
#----------------------------------------------------------------------


# ---------------------------------------------------
# procedure's block
# ---------------------------------------------------
proc new_file {} {
    # $button is used for dialogs
    global XCState xcMisc system

    if { ! $system(c95_exist) } { return }

    set xcMisc(titlefile) Unsaved
    wm title . "XCrySDen: $xcMisc(titlefile)"

    # define a state variable XCState as c95
    set XCState(state) c95_newinput
    xcUpdateState; # XCState(state) has changed, so -> xcUpdateState
    xcAdvGeomState reset
    ChooseSpecies
}


proc ChooseSpecies {} {
    global XCState

    # this procedure creates window: Choose
    xcToplevel .cho "XInput: Choose" "SeqTable" . 100 50
    AlwaysOnTopON . .cho

    # TOP & BOTTOM & CLOSE frames; all widget within .cho will be displayed 
    # in one of those
    frame .cho.top -relief raised -bd 1
    frame .cho.bot -relief raised -bd 1
    frame .cho.close -relief raised -bd 1
    pack .cho.top -side top -fill both -padx 1 -pady 1
    pack .cho.bot .cho.close -side top -fill both

    # in the TOP there will be label, entry & message
    label .cho.top.label -text "Title:"
    entry .cho.top.entry -relief sunken -textvariable job_title
    label .cho.top.lab -text "Title is just a comment;\nnot used otherwise"\
	    -anchor center
    
    pack .cho.top.lab -side bottom -fill both -padx 5 -pady 5 -expand 1
    pack .cho.top.label -side left -fill x -padx 3 -pady 10 -expand 1
    pack .cho.top.entry -side left -fill x -padx 3 -pady 10 -expand 1
    
    #set focut to .CHO.TOP.ENTRY
    focus .cho.top.entry

    # in the BOT there will go label & five buttons (Molecule,Polymer,
    #                                           Slab,Crystal,External)

    label .cho.bot.label -text "Choose one of the following:" 
    button .cho.bot.button1 -text "Molecule" -command [list Molecule .cho]
    button .cho.bot.button2 -text "Polymer" -command [list Polymer .cho]
    button .cho.bot.button3 -text "Slab" -command [list Slab .cho] 
    button .cho.bot.button4 -text "Crystal" -command [list Crystal .cho]

    button .cho.close.but -text "Cancel" \
	    -command [list NewInputCancelProc .cho]

    pack .cho.bot.label -side top -fill both -expand true -padx 10 -pady 15
    pack .cho.bot.button1 .cho.bot.button2 .cho.bot.button3 \
	    .cho.bot.button4 \
	    -side top -fill x -expand 1 -padx 25 -pady 10 -ipady 2
    pack .cho.close.but -expand 1 -padx 25 -pady 5
    return
}


proc Molecule {tplw} {
    # do we need all those variables to be global?????
    global type_group type_group1 species group_list groupsel inp

    # $tplw is the name of toplevel widget from whom we arrive in this toplevel
    AlwaysOnTopOFF 
    destroy $tplw

    set w .mol
    set type_group point
    set type_group1 Point
    set species molecule 
    # now "load" a group_list variable; it's done by point_group PROC
    point_group
    set inp(IFHR) 0
    #now we are able to call geom_sym_input proc
    geom_sym_input $w $tplw
}


proc Polymer {tplw} {
    global type_group type_group1 species group_list groupsel inp
    
    # $tplw is the name of toplevel widget from whom we arrive in this toplevel
    AlwaysOnTopOFF 
    destroy $tplw

    set w .poly
    set type_group line
    set type_group1 Line
    set species polymer 
    
    # now "load" a group_list variable; it's done by line_group PROC
    line_group 
    set inp(IFHR) 0
    #now we are able to call geom_sym_input proc
    geom_sym_input $w $tplw
}

proc Slab {tplw} {
    # do we need all those variables to be global?????
    global type_group type_group1 species group_list groupsel inp
    
    # $tplw is the name of toplevel widget from whom we arrive in this toplevel
    AlwaysOnTopOFF 
    destroy $tplw

    set w .slab
    set type_group plane
    set type_group1 Plane
    set species slab
    
    # now "load" a plane_list variable; it's done by plane_group PROC
    plane_group 
    set inp(IFHR) 0
    #now we are able to call geom_sym_input proc
    crys_slab_sym $w $tplw
}


proc Crystal {tplw} {
    # do we need all those variables to be global?????
    global type_group type_group1 species group_list groupsel inp
    
    # $tplw is the name of toplevel widget from whom we arrive in this toplevel
    AlwaysOnTopOFF 
    destroy $tplw

    set w .crys
    set type_group space
    set type_group1 Space
    set species crystal 
    # now "load" a group_list variable; it's done by space_group PROC
    space_group 
    set inp(IFLAG) 1
    set inp(IFHR)  0
    #now we are able to call geom_sym_input proc
    crys_slab_sym $w tplw
}


proc crys_slab_sym {w tplw} {
    global done okay type_group type_group1 species group_list \
	    groupsel n_groupsel \
	    crdatom 

    set done {}
    set okay {}

    xcToplevel $w "Crystal95's Input: $species" "Input: $species" . 100 50
    AlwaysOnTopON . $w

    frame $w.frame1 -relief raised -bd 1
    frame $w.frame2 -relief raised -bd 1
    frame $w.frame3
    
    pack $w.frame1 $w.frame2 $w.frame3 -side top -fill both -expand 1    
    
    #----------------------------
    # FRAME #1
    label $w.frame1.lab -text "Symmetry input for $species"
    pack $w.frame1.lab -expand 1 -fill both -padx 10 -pady 10

    #---------------------------
    # widgets in FRAME #2
    frame $w.frame2.frm
    label $w.frame2.frm.label -text "$type_group1 group of $species:"
    entry $w.frame2.frm.entry -relief sunken -textvariable groupsel
    button $w.frame2.button -text "Select a $type_group group" \
	    -command [list CheckGroup button]
    pack $w.frame2.frm -side top
    pack $w.frame2.button -side bottom -padx 10 -pady 7
    pack $w.frame2.frm.label $w.frame2.frm.entry -side left \
	    -padx 1 -pady 7 -expand 1
    
    # now we made focus & binding for entry
    focus $w.frame2.frm.entry 
    bind $w.frame2.frm.entry <Return> { CheckGroup entry }

    #---------------------------------------------------------
    # widgets in FRAME #3
    # --------------------------------------------------------
    button $w.frame3.but_prev -text "<< Previous"  \
	    -command "CancelProc $w; ChooseSpecies"
    button $w.frame3.but_next -text "Next >>" -command gotoGeom
    pack $w.frame3.but_prev $w.frame3.but_next \
	    -side left -anchor c -padx 10 -pady 10 -expand 1

    
    tkwait variable okay
    set tplw $w
    AlwaysOnTopOFF 
    destroy $w
    set w ${w}1
    geom_sym_input $w $tplw
}

proc CheckGroup {{act {}}} {
    # this proc is called by CRYS_SLAB_SYM
    global done okay type_group type_group1 species group_list \
	    groupsel n_groupsel crdatom inp var 

    set done "" 
    #puts stdout "CG.1> $act"
    if { $act == "button" } {
	#puts stdout "CG.2> $act"
 	group_sel .sel_group $type_group $group_list
    } else {
	group_sel_ok $group_list
    }

    if { $done == 1 } {
	# $rhombohexa var. is for rhombo. groups
	# rhombo. groups are all that are R-centered
	
	# this need to be modified T.K_MOD
	if { [lindex $groupsel 0] == "R" } {
	    # query for rhombo./hexa. axes
	    xcToplevel .ifhr "IFHR" "IFHR" . 120 70 1
	    catch { grab .ifhr }
	    
	    frame .ifhr.frame -relief raised -bd 2
	    label .ifhr.frame.lab -text "Select the type of cell for\n\
		    Rhombohedral Group"
	    pack .ifhr.frame .ifhr.frame.lab -side top -fill both \
		    -padx 1 -pady 1
	    # from IFHR we must assign VAR:
	    switch -exact -- $inp(IFHR) {
		0 { set var "hexagonal cell" }
		1 { set var "rhombohedral cell" }
	    }
	    
	    RadioButtons .ifhr.frame var top "hexagonal cell" \
		    "rhombohedral cell"
	    
	    button .ifhr.but -text OK -command {set ok 1}
	    pack .ifhr.but -pady 7 -expand 1
	    
	    # when OK button is pressed toplevel can be destroyed
	    tkwait variable ok
	    destroy .ifhr
	    
	    # now we must assing a value for $IFHR, which is suitable for 
	    # crystal program
	    switch -glob -- $var {
		hexa*  { set inp(IFHR) 0 }
		rhomb* { set inp(IFHR) 1 }
	    }
	}	
    }
}


proc gotoGeom {{w {}}} {
    global okay done type_group type_group1 species group_list groupsel \
	    crdatom XCState

    group_sel_ok $group_list
    if { $done == 1 } {
	if { [xcIsActive c95] && [xcIsActive openinput] } {
	    if [winfo exists $w] { 
		AlwaysOnTopOFF 
		destroy $w 
	    }
	    set okay 1
	} else {
	    set okay 1
	}
    }
}

    
##############################################################################
proc geom_sym_input {w tplw} {
    # do we need all those variables to be global?????
    global err okay done type_group type_group1 species class \
	    group_list groupsel n_groupsel inp crdatom XCState

    # $done is used with bind ENTRY <Return> {set done 1}
    # this means when you press return in an entry textvariable is set !!!!!
    set done 0
    if ![info exists class] {set class {}}
    if ![info exists class] {set crdatom 0}

    xcToplevel $w "Crystal95's Input: $species" "Input" . 100 50
    AlwaysOnTopON . $w

    frame $w.frame1 -relief raised -bd 1
    frame $w.frame2 -relief raised -bd 1
    frame $w.frame3 -relief raised -bd 1
    frame $w.frame4 -relief raised -bd 1    
    frame $w.frame5 -bd 1
    frame $w.frame6 -bd 1
    
    # now we will pack this frames
    # $species IS IMPORTANT, BECAUSE THE DISPLAY IS DIFFERENT 
    # FOR EACH SPECIES

    if { $species == "molecule" } {
	# there is no FRAME #3 for MOLECULES
	pack $w.frame1 $w.frame2 $w.frame4 $w.frame5 \
		-side top -fill both -expand 1
	pack $w.frame6 -padx 40 \
		-side top -fill both -expand 1
    } elseif { $species == "polymer" } {
	#pack everything
	pack $w.frame1 $w.frame2 $w.frame3 $w.frame4 $w.frame5 $w.frame6 \
		-side top -fill both -expand 1
	pack $w.frame6 -padx 40 \
		-side top -fill both -expand 1
    } else {
	pack $w.frame1  $w.frame3 $w.frame4 $w.frame5 $w.frame6 \
		-side top -fill both -expand 1
	pack $w.frame6 -padx 40 \
		-side top -fill both -expand 1
    }

    #---------------------------
    # widget in FRAME #1
    label $w.frame1.lab -text "Geometry input for $species"
    pack $w.frame1.lab -expand 1 -fill both -padx 10 -pady 10
    
    #---------------------------
    # widgets in FRAME #2; ONLY FOR MOLECULE & POLYMER
    if {$species == "molecule" || $species == "polymer"} {
	frame $w.frame2.frm
	label $w.frame2.frm.label -text "$type_group1 group of $species:"
	entry $w.frame2.frm.entry -relief sunken -textvariable groupsel
	button $w.frame2.button -text "Select a $type_group group" \
		-command [list group_sel .sel_group $type_group $group_list]
	pack $w.frame2.frm -side top
	pack $w.frame2.button -side bottom -padx 10 -pady 7
	pack $w.frame2.frm.label $w.frame2.frm.entry -side left \
		-padx 1 -pady 7 -expand 1
    
	# now we made focus & binding for entry
	focus $w.frame2.frm.entry 
	bind $w.frame2.frm.entry <Return> { group_sel_ok $group_list }
    }
    
    # IF----IF----IF----IF----IF----IF----IF----IF
    # FRAME #3 olny for $species != "molecule"
    #-------------------------------------------------------------
    # widgets in FRAME #3
    #
	
    if { $species != "molecule" } {
	
	# only for CRYSTALS; about ORIGIN SHIFT (var. inp(IFSO))
	# default value for IFSO is set to 0
	set inp(IFSO) 0
	if { $species == "crystal" } {
	    button $w.frame3.button -text "Set the ORIGIN of the \n\
		    crystal reference frame" -command set_origin 
	    pack $w.frame3.button -side top -ipadx 3 -ipady 3 -pady 7 \
		    -expand 1
	}
	
	label $w.frame3.lab -text \
		"Minimal set of lattice parameters"
	pack $w.frame3.lab -side top -expand 1 -fill both -padx 10 -pady 5
	
	# we must read MINIMAL SET OF LATTICE PARAMETERS
	if ![info exists inp(A)] {set inp(A) {}}
	if ![info exists inp(B)] {set inp(B) {}}
	if ![info exists inp(C)] {set inp(C) {}}
	if ![info exists inp(ALFA)] {set inp(ALFA) {}}
	if ![info exists inp(BETA)] {set inp(BETA) {}}
	if ![info exists inp(GAMMA)] {set inp(GAMMA) {}}
	read_param $w.frame3  
	frame $w.frame3.frm 
    }
    
    #---------------------------------------------------------
    # widgets in FRAME #4    

    # do we need this !!!
    if ![info exists inp(NATR)] {set inp(NATR) {}}

    label $w.frame4.lab -anchor center -text \
	    "Specify atomic numbers and coordinates\nof non-equivalent atoms"
    button $w.frame4.but -text "Atomic numbers/coordinates" \
	    -command atom_num_coord
    pack $w.frame4.lab -side top -padx 10 -pady 7 -fill both -expand true
    pack $w.frame4.but -side top -padx 10 -pady 7 \
	    -ipadx 3 -ipady 3 -expand true
    
    #---------------------------------------------------------
    # widgets in frame #5
    button $w.frame5.but_view -text "View $species" \
	    -command [list contrl_var $w View]    
    # is this OK!!!!!!!!!!!!
    button $w.frame5.but_prev -text "<< Previous"  \
	    -command [list geom_sym_prev $w $tplw]
    #button $w.frame5.but_next -text "Next >>" \
	    #-command [list contrl_var $w Next]
    puts stdout "GEOM_SYM_INPUT: \$w: $w"
    
    pack $w.frame5.but_view $w.frame5.but_prev \
	    -side left -anchor c -padx 10 -pady 10 -expand 1
    
    #----------------------------------------------------------
    # widgets in frame #6
    # button Done is displayed only if XCState == c95_openinput
    if { [xcIsActive c95] && [xcIsActive openinput] } {
	button $w.frame6.but1 -text "Done" -command [list contrl_var $w Modify]
	pack $w.frame6.but1 -side left -anchor c \
		-padx 10 -pady 10 -expand 1
    } else {
	button $w.frame6.but2 -text "Cancel" \
		-command [list NewInputCancelProc $w]
	pack $w.frame6.but2 -side left -anchor c \
		-padx 10 -pady 10 -expand 1
    }
            
    #__________________________________________________________
    # END OF GEOMETRY SPECIFIERS BLOCK
    #----------------------------------------------------------
    
    tkwait variable okay
    CancelProc $w
    # maybe XCState(toplevel) is not empty
    #if [info exists XCState(toplevel)] {
    #	 if { [winfo exists $XCState(toplevel)] || $XCState(toplevel) != "" } {
    #	     AlwaysOnTopON .
    #	 }
    #}
}


proc geom_sym_prev {w tplw} {
    global done type_group type_group1 species group_list groupsel

    if [winfo exists $w] { CancelProc $w }

    if { $tplw == ".cho" } {
	ChooseSpecies
    } elseif  { $tplw == ".crys" || $tplw == ".slab" || \
	    ( $tplw == ".opflmod" && ( $species == "slab" ||\
	    $species == "crystl" ) ) } {
	crys_slab_sym $tplw $tplw; # THIS IS OK; 
    } 
}


proc set_origin {} {
    # leave global VAR there, It's needed !!!!!!!!!!!!
    global var button inp

    xcToplevel .ifso "Origin Settings" "Origin Settings" . 120 70 1
    catch { grab .ifso }

    frame .ifso.frame -relief raised -bd 2
    pack .ifso.frame -side top -fill both -padx 1 -pady 1

    # from IFSO we must assign VAR:
    switch -exact -- $inp(IFSO) {
	0 { set var "Origin derived from the symbol (IFSO=0)" }
	1 { set var "Standard shift of origin (IFSO=1)" }
	2 { set var "Non-standard shift of origin (IFSO=2)" }
    }

    RadioButtons .ifso.frame var top \
	"Origin derived from the symbol (IFSO=0)" \
	"Standard shift of origin (IFSO=1)" \
	"Non-standard shift of origin (IFSO=2)" 

    button .ifso.but -text OK -command {set ok 1}
    pack .ifso.but -pady 7 -expand 1
    
    # when OK button is pressed toplevel can be destroyed
    tkwait variable ok
    destroy .ifso


    # now we must assing a value for $IFSO, which is suitable for 
    # crystal program
    switch -glob -- $var {
	Origin*       { set inp(IFSO) 0 }
	Standard*     { set inp(IFSO) 1 }
	Non-standard* { set inp(IFSO) 2 }
    }
    
    # if inp(IFSO)>=1 --> specify coordinates of non-standard shift
    if { $inp(IFSO) == 2 } { orig_coord }
    
}


proc orig_coord {} {
    global inp ok err geo

    xcToplevel .orig_coord "Origin Shift" "Origin Shift" . 120 70 1
    catch { grab .orig_coord }
    # this procedure create toplevel, where coord. of origin-shift 
    # will be specified
    
    frame .orig_coord.frame -relief raised -bd 2
    label .orig_coord.frame.lbl -text "Specify coordinates of Origin-Shift"
    pack .orig_coord.frame -fill both -expand 1
    pack .orig_coord.frame.lbl -side top -padx 7 -pady 10 

    Entries .orig_coord.frame "Fract-X*24: Fract-Y*24: Fract-Z*24:"\
	    "inp(IX) inp(IY) inp(IZ)" 6
    # DO WE NEED $FOC HERE !!!!! -TEST-TEST
    focus .orig_coord.frame.frame.entry1
    button .orig_coord.but -text OK -command [list check_var {{inp(IX) int}\
	    {inp(IY) int} {inp(IZ) int}} {.orig_coord.frame.frame.entry1 \
	    .orig_coord.frame.frame.entry2 .orig_coord.frame.frame.entry3}]
    pack .orig_coord.but -side top -expand 1 -pady 7
    
    # will this work
    tkwait variable ok
    destroy .orig_coord
}


proc n_atoms {} {
    global inp crdatom
    # this PROC display a TOPLEVEL where a number of non-eq. atoms is specefied!
    # NATR -number of non. eq. atoms

    xcToplevel .natom "N. of Atoms" "N. of Atoms" . 120 70 1
    catch { grab .natom }

    frame .natom.top -relief raised -bd 1
    frame .natom.bot
    pack .natom.top .natom.bot -side top -expand 1 -fill both
    
    label .natom.top.lab -text "Number of non-equivalent atoms:"
    entry .natom.top.entry -relief sunken -width 4 -textvariable inp(NATR)
    pack .natom.top.lab -side top -padx 10 -pady 7 -expand 1
    pack .natom.top.entry -side bottom   -pady 7 -expand 1
    bind .natom <Return> [list check_var \
	    {{inp(NATR) posint}} .natom.top.entry]
    focus .natom.top.entry
    
    button .natom.bot.but -text OK -command [list check_var \
	    {{inp(NATR) posint}} .natom.top.entry] 
    pack .natom.bot.but -side top -pady 7 -expand 1
    bind .natom.bot.but <Return> [list check_var \
	    {{inp(NATR) posint}} .natom.top.entry]

    tkwait variable ok
    destroy .natom
}

    
proc atom_num_coord {} {
    global inp species ok crdatom
    # NATR - num. of non-equivalent atoms
    
    # n_atom is PROC, where $inp(NATR) is specefied
    n_atoms

    # just in any case
    set ok {}
    # if $inp(NATR) == {} -> first specify inp(NATR)
    # WARNING - what if inp(NATR) is not a number
    if { $inp(NATR) == "" } {
	dialog .atom1 WARNING "You must first specify a number of \
		non-equivalent atoms" warning 0 OK	
    } else {
	xcToplevel .atom "Atom Coordinates" "Atom Coordinates" . 120 70 1
	catch { grab .atom }

	# bottom frame where OK button will be
	set Frm [frame .atom.frame]
	pack .atom.frame -side bottom -expand true -fill both 
	# and one frame where canvas&scrollbar will be!!
	frame .atom.f -relief sunken -bd 2
	pack .atom.f -side top -expand true -fill both 
	
	canvas .atom.f.canv -yscrollcommand [list .atom.f.yscroll set]
	scrollbar .atom.f.yscroll -orient vertical -command \
		[list .atom.f.canv yview]
	pack .atom.f.yscroll -side right -fill y
	pack .atom.f.canv -side left -fill both -expand true
	
	# create FRAME to hold every LABELS & ENTRIES
	set f [frame .atom.f.canv.f -bd 0]
	.atom.f.canv create window 0 0 -anchor nw -window $f

	for {set i 1} {$i <= $inp(NATR)} {incr i 1} {	    
	    frame $f.fr$i -relief groove -bd 2
	    pack $f.fr$i -padx 5 -pady 5
	    label $f.fr${i}.label$i -text "Atom N.: $i" 
	    pack $f.fr${i}.label$i -anchor w -padx 7 -pady 7
	    frame $f.fr${i}.frm$i
	    frame $f.fr${i}.frm1_$i
	    frame $f.fr${i}.frm2_$i	    
	    frame $f.fr${i}.frame$i 
	    pack $f.fr${i}.frm1_$i $f.fr${i}.frm2_$i \
		    -in $f.fr${i}.frm$i -side left -anchor w
	    pack $f.fr${i}.frm$i $f.fr${i}.frame$i -side top -anchor w
	    Entries $f.fr${i}.frm1_$i {{Atomic Number:}} inp(NAT,$i) 5
	    button $f.fr${i}.frm2_$i.b -text "Periodic Table" \
		    -command [list ptable .atom \
		    -command ptableSelectAtomicNumber -variable inp(NAT,$i)]
	    pack $f.fr${i}.frm2_$i.b -side left -padx 5
	    # make a varlist & foclist for PROC CHECK_VAR
	    lappend varlist [list inp(NAT,$i) nat]
	    lappend foclist $f.fr${i}.frm1_$i.frame.entry1
	    switch -exact -- $species {
		molecule { 
		    Entries $f.fr${i}.frame$i \
			    "X: Y: Z:" " inp(X,$i) inp(Y,$i) inp(Z,$i)" 8 
		    # make a varlist for PROC CHECK_VAR
		    lappend varlist inp(X,$i) inp(Y,$i) inp(Z,$i)
		}
		polymer { 
		    Entries $f.fr${i}.frame$i \
			    "Fractional-X: Y: Z:" \
			    " inp(X,$i) inp(Y,$i) inp(Z,$i)" 8 
		    # make a varlist for PROC CHECK_VAR
		    lappend varlist [list inp(X,$i) fract] inp(Y,$i) inp(Z,$i)
		}
		slab { 
		    Entries $f.fr${i}.frame$i "Fractional-X: \
			    Fractional-Y: Z:" " inp(X,$i) inp(Y,$i) inp(Z,$i)" 8 
		    # make a varlist for PROC CHECK_VAR
		    lappend varlist [list inp(X,$i) fract] \
			    [list inp(Y,$i) fract] inp(Z,$i)
		} 
		crystal { 
		    Entries $f.fr${i}.frame$i "Fractional-X: \
			    Fractional-Y: Fractional-Z:" \
			    " inp(X,$i) inp(Y,$i) inp(Z,$i)" 8 
		    # make a varlist for PROC CHECK_VAR
		    lappend varlist [list inp(X,$i) fract] \
			    [list inp(Y,$i) fract] [list inp(Z,$i) fract]
		}
	    }
	    # make foclist for PROC CHECK_VAR
	    lappend foclist $f.fr${i}.frame$i.frame.entry1 \
		    $f.fr${i}.frame$i.frame.entry2 \
		    $f.fr${i}.frame$i.frame.entry3    
	}
	
	puts stdout "FOCLIST: $foclist\n\n"
	puts stdout "VARLIST: $varlist"
	set child [lindex [pack slaves $f] 0]
	
	# set the focus to first entry that upper FOR-LOOP create
	focus $f.fr1.frm1_1.frame.entry1

	tkwait visibility $child
	set width [winfo width $f]
	set height [winfo height $f]
	if { $inp(NATR) < 5 } {
	    .atom.f.canv config -width $width -height $height 
	} else {
	    .atom.f.canv config -width $width -height \
		    [expr $height / $inp(NATR) * 4] \
		    -scrollregion "0 0 $width $height"
	}

	button $Frm.butok -text OK -command \
		[list check_var $varlist $foclist]

	button $Frm.butcan -text "Cancel" -command "destroy .atom"

	bind $Frm.butok <Return> [list check_var $varlist $foclist]
	pack $Frm.butok $Frm.butcan -side left \
		-expand 1 -padx 10 -pady 10
	
	tkwait variable ok

	#===========================================
	# THIS IS FOR CHEACKING IN THE GEM_SYM_INPUT
	set crdatom 1
	#===========================================
	destroy .atom
    }
}


proc contrl_var {topl action} {
    global crdatom err okay groupsel group_list class species inp XCState
    # topl -name of toplevel

    # PROC wait for variable okay; then exit (this means everything is OK
    # check entries in FRAME #2 --> $groupsel
    ########
    
    set err 0
    group_sel_ok $group_list

    puts stdout "CLASS, species> $species $class"
    flush stdout
    if { $err } {
	focus $topl.frame2.frm.entry
    } 
    if { $err == 0 } {
	# if there is no error in specifying $groupsel, we could proceed
	# check entries in FRAME #3 --> Min. lattice parameter list
	set n 1
	if {$species != "molecule"} {
	    foreach elem $class {
		foreach elem1 $elem { 
		    lappend varlist inp($elem1) 
		}
		for {set i 1} {$i <= [llength $elem]} {incr i} {
		    lappend foclist $topl.frame3.$n.frame.entry$i
		}
		incr n
	    }
	    # check_var return OK variable; this is not good this time, so:
	    puts stdout "CONTRL_VAR:: $varlist $foclist"
	    check_var $varlist $foclist	    
	}
    } 
    if { $err == 0 } {
	###############
	# CHECKING FRAME #4
	if ![info exists crdatom] {
	    tk_dialog [WidgetName] ERROR "ERROR: You forget to specify atomic \
		    coordinates" error 0 OK
	    set err 1
	    return 0
	}
	    
	if { $crdatom != "1" } {
	    dialog .crdatm ERROR "ERROR: Your specification of atoms \
		    coordinates is bad. Try again !!!" error 0 OK
	    set err 1
	}
    }
    
    #============================================
    if { $err == 0 } { 
	# THIS IS PROBABLY OBSOLETE; CHECK THIS
	# if there were no ERROR($err=0) --> set ok 1
	set okay 1 
	puts stdout $action
	if { $action == "View" } {
	    puts stdout "GOING TO CalStru"
	    #append XCState(state) "_render"; CalStru takes care of that
	    CalStru
	}
    }
    #============================================ 
    
    #puts stdout "RESULT\n\
    #	     ============="
    #puts stdout "GROUPSEL: $groupsel"
    #puts stdout "inp(IFSO): $inp(IFSO)"
    #puts stdout "inp(IX),inp(IY),inp(IZ): $inp(IX) $inp(IY) $inp(IZ)"
    #puts stdout "A: $inp(A)   B: $inp(B)   C: $inp(C)\n\
    #inp(ALFA): $inp(ALFA)    inp(BETA): $inp(BETA)    inp(GAMMA): $inp(GAMMA)"
    #puts stdout "inp(NATR): $inp(NATR)"
    #for {set i 1} {$i <= $inp(NATR)} {incr i} {
    #	 puts stdout "$Nat($i)\n \
    #		 $X($i)     $Y($i)     $Z($i)"  
    #}
}


proc NewInputCancelProc {w} {
    global XCState

    CancelProc $w
    if { [info exists XCState] } {
	unset XCState
    }
    xcUpdateState
}


#proc InputBuilder {} {
#    global button
#
#    set button {}
#    label   .xfile_top -relief raised -text \
#	     "XCrySDen\n\nA X-based Interactive Input Builder for Crystal95" \
#	     -font "-adobe-helvetica-bold-o-normal--25-180-100-100-p-138-iso8859-1"
#    button .xfile_new -text "New File" -command new_file 
#    button .xfile_open -text "Open File" -command OpenCyrstalInput
# 
#    button .xfile_exit -text "Exit" -command exit_pr
#
#    if { $button == 0 } {exit} 
#
#    pack .xfile_top -side top -fill x -ipadx 3m -ipady 10m -anchor c 
#    pack .xfile_new .xfile_open .xfile_exit  \
#	     -side left -anchor c -padx 20m -pady 3m -expand 1
#}
