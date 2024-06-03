#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/undoAdvGeom.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc UndoMenuMotion {menu y} {
    global undoMenu

    set n [$menu index $y]

    # it may easily happen that $n is not an integer
    if [catch {expr abs($n)}] {
	return
    }

    for {set i 0} {$i < $n} {incr i} {	
	$menu entryconfigure $i \
		-foreground $undoMenu(active_fg) \
		-background $undoMenu(active_bg)
    }
    set m [$menu yposition end]
    for {set i $n} {$i < $m} {incr i} {
	$menu entryconfigure $i \
		-foreground $undoMenu(default_fg) \
		-background $undoMenu(default_bg)
    }
}


proc Undo {menu redo index} {
    global undoAdvGeom

    set n $index

    set ld [expr $n - $undoAdvGeom(start_index)]    

    set undoAdvGeom(redo_list) \
	    [GetUndoRedoList [lrange $undoAdvGeom(list) 0 $ld] \
	    $undoAdvGeom(redo_list)]
    $undoAdvGeom(menu) entryconfigure "Redo ..." -state normal

    set undoAdvGeom(list) [lrange $undoAdvGeom(list) [expr $ld + 1] end]
    if { $undoAdvGeom(list) == {} } {
	$undoAdvGeom(menu) entryconfigure "Undo ..." -state disabled
    }
    CreateUndoRedoMenu $undoAdvGeom(list) $menu $redo Undo \
	    undoAdvGeom(start_index) undoAdvGeom(current_index)
    CreateUndoRedoMenu $undoAdvGeom(redo_list) $redo $menu Redo \
	    undoAdvGeom(redo_start_index) undoAdvGeom(redo_current_index)

    # insert the real UNDO command here !!!
    xcAdvGeomState undo [expr $ld + 1]
    CalStru
    ###
}


proc GetUndoRedoList {listnew listold} {
    #
    # we must revert the order of newlist
    #

    set n [llength $listnew]
    set redo_list {}
    for {set i $n} {$i > 0} {incr i -1} {
	set ii [expr $i - 1]
	lappend redo_list [lindex $listnew $ii]
    }
    return [concat $redo_list $listold]
}


proc Redo {redo menu index} {
    global undoAdvGeom

    set n $index

    set ld [expr $n - $undoAdvGeom(redo_start_index)]

    set undoAdvGeom(list) \
	    [GetUndoRedoList [lrange $undoAdvGeom(redo_list) 0 $ld] \
	    $undoAdvGeom(list)]
    $undoAdvGeom(menu) entryconfigure "Undo ..." -state normal

    set undoAdvGeom(redo_list) \
	    [lrange $undoAdvGeom(redo_list) [expr $ld + 1] end]
    if { $undoAdvGeom(redo_list) == {} } {
	$undoAdvGeom(menu) entryconfigure "Redo ..." -state disabled
    }
    CreateUndoRedoMenu $undoAdvGeom(list) $menu $redo Undo \
	    undoAdvGeom(start_index) undoAdvGeom(current_index)
    CreateUndoRedoMenu $undoAdvGeom(redo_list) $redo $menu Redo \
	    undoAdvGeom(redo_start_index) undoAdvGeom(redo_current_index)

    # insert the real REDO command here !!!
    xcAdvGeomState new [expr $ld + 1]
    CalStru

}


proc GenCommUndoRedo item {
    global undoAdvGeom undoMenu

    $undoAdvGeom(menu) entryconfigure "Undo ..." -state normal

    ##########################
    #generate command HERE !!!
    #^^^^^^^^^^^^^^^^^^^^^^^^^
    set item [list $item]
    set undo $undoMenu(undo)
    set redo $undoMenu(redo)
    #-----------
    set undoAdvGeom(list) [concat "$item" $undoAdvGeom(list)]
    CreateUndoRedoMenu $undoAdvGeom(list) $undo $redo Undo \
	    undoAdvGeom(start_index) undoAdvGeom(current_index)
    
    # delete redo-list
    $undoAdvGeom(menu) entryconfigure "Redo ..." -state disabled
    set undoAdvGeom(redo_list) {}
    $redo delete 0 end
}


proc CreateUndoRedoMenu {list menu1 menu2 command start_index current_index} {
    global undoMenu
    upvar #0 $start_index start $current_index current
    # first delete the whole list, then rebuilt
    $menu1 delete 0 end

    set start $current
    foreach item $list {
	$menu1 add command -label $item \
		-command [list $command $menu1 $menu2 $current] \
		-foreground $undoMenu(default_fg) \
		-background $undoMenu(default_bg) \
		-activeforeground $undoMenu(active_fg) \
		-activebackground $undoMenu(active_bg)
	
	incr current
    }
}

proc ClearUndoRedoBuffer {} {
    global undoAdvGeom
    
    set undoAdvGeom(start_index)        0
    set undoAdvGeom(current_index)      0
    set undoAdvGeom(list)               {}
    set undoAdvGeom(redo_start_index)   0
    set undoAdvGeom(redo_current_index) 0
    set undoAdvGeom(redo_list)          {}
}

#set m [menubutton .advgeom -text AdvGeom -menu .advgeom.ur]
#
#set menu [menu .advgeom.ur -tearoff 0]
#
#$menu add cascade -label Undo -menu $menu.undo
#$menu add cascade -label Redo -menu $menu.redo
#$menu entryconfig Undo -state disabled
#$menu entryconfig Redo -state disabled
#set undoAdvGeom(menu) $menu
#
#
#set undo  [menu $menu.undo -tearoff 0]
#set redo  [menu $menu.redo -tearoff 0]
#button .b -text "Generate Command" -command [list GenComm $undo $redo]
#pack $m .b -side left
#
#
#bind $undo <Motion> [list UndoMenuMotion $undo @%y]
#bind $redo <Motion> [list UndoMenuMotion $redo @%y]
#
#flush stdout





