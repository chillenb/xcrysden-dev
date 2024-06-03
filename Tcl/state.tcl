#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/state.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

##############################################################################
# DISABLE WIDGETS:
# ---------------
# disable all children of parent $w if possible
#proc xcDisableAll {wlist} {
#    
#    foreach w $wlist {
#	 if ![winfo exists $w] continue
#	 set children [winfo children $w]
#	 if { $children == "" } {
#	     #puts stdout $w 
#	     catch [list $w configure -state disabled]
#	 } else {
#	     foreach child $children {
#		 xcDisableAll $child
#	     }
#	 }
#    }	
#}
proc xcDisableEntry {args} {
    global xcColors
    foreach entry $args {
	$entry config -bd 0 -fg $xcColors(disabled_fg) -state disabled 
    }	
}

proc xcEnableEntry  {args} {
    global xcColors
    foreach entry $args {
	$entry config -state normal -bd 2 -fg $xcColors(enabled_fg)
    }	
}

proc xcAll:wlist {wlist} {
    global XCstate
    
    if { [lindex $wlist 0] == "-disabledfg" } {
	set XCstate(-disabledfg) 1
	return [lrange $wlist 1 end]
    } else {
	set XCstate(-disabledfg) 0
	return $wlist
    }
} 

proc xcDisableAll {args} {
    set wlist [xcAll:wlist $args]
    _xcDisableAll $wlist
}
proc _xcDisableAll {wlist} {
    global xcColors XCstate    
    foreach w $wlist {
	if { ![winfo exists $w] } continue
	set children [winfo children $w]
	if { $children != "" } {
	    foreach child $children {
		_xcDisableAll $child
	    }
	}
	catch [list $w configure -state disabled]
	if { $XCstate(-disabledfg) } {
	    catch [list $w configure -fg $xcColors(disabled_fg)]
	}
    }
}


proc xcDisableMenuOne {menu list} {
    foreach {w e} $list {
	lappend wlist $w
	lappend elist $e
    }
    xcDisableOne $wlist
    xcDisableMenuentryOne $menu $elist
}


# disable only widgets in $wlist
proc xcDisableOne {wlist} {
    foreach wid $wlist {
	catch [list $wid configure -state disabled]
    }
}

proc xcDisableMenuentryOne {menu elist} {
    foreach entry $elist {
	catch {$menu entryconfigure $entry -state disabled}
    }
}
    

##############################################################################
# ENABLE WIDGETS:
# --------------
# enable all children of parent $w if possible
proc xcEnableAll {args} {
    _xcEnableAll [xcAll:wlist $args]
}
proc _xcEnableAll {wlist} {
    global xcColors XCstate    
    
    foreach w $wlist {
	if { ![winfo exists $w] } { continue }
	set children [winfo children $w]
	if { $children != "" } {
	    foreach child $children {
		_xcEnableAll $child
	    }
	}
	catch [list $w configure -state normal]
	if { $XCstate(-disabledfg) } {
	    catch [list $w configure -fg $xcColors(enabled_fg)]
	}
    }
}

proc xcEnableMenuOne {menu list} {
    foreach {w e} $list {
	lappend wlist $w
	lappend elist $e
    }
    xcEnableOne $wlist
    xcEnableMenuentryOne $menu $elist
}


# enable only widgets in $wlist
proc xcEnableOne {wlist} {
    foreach wid $wlist {
	catch [list $wid configure -state normal]
    }
}

proc xcEnableMenuentryOne {menu elist} {
    foreach entry $elist {
	catch {$menu entryconfigure $entry -state normal}
    }
}

proc xcConfigAll {wlist args} {
    foreach w $wlist {
	if { ![winfo exists $w] } continue
	set children [winfo children $w]
	if { $children != "" } {
	    foreach child $children {
		xcConfigAll $child $args
		catch {eval $child configure $args}
	    }
	}
    }
}

proc xcIsActive {word} {
    global XCState

    #--------
    # states:
    #--------
    # c95         .... c95 case (newinput, openinput, ..
    # newinput    .... new input will be/has been made
    # openinput   .... existeing input is opened
    # render      .... structure is rendered
    # properties  .... crytal fort.9 is opened
    # wien        .... when dealing with WIEN2k
    # multislab   .... when dealing with WIEN2k multi-slab
    # external    .... when reading structre from some external file format
    # external34  .... when reading structre from some external file format 
    #                  via FORT.34; !!! *** probably temporal *** !!!
    if { ![info exists XCState(state)] } { return 0 }
    return [string match *${word}* $XCState(state)]
}


############################################################################
proc xcUpdateState {{menu .menu} {menu2 .menu}} {
    global XCState dispmode periodic radio sInfo system xsfAnim xcMisc

    # disable all entries but exit&print-setup-entry in file-menu 
    foreach entry { 
	{New CRYSTAL Input} 
	{Open Structure ...} 
	{Open PWscf ...} 
	{Open CRYSTAL ...} 
	{Open WIEN2k ...} 
	Close 
	{Save XSF Structure} 
	{Save Current State and Structure}
	{Save Current State}
	{Save CRYSTAL Input} 
	{Save WIEN2k Struct File} 
	Print* 
	{XCrySDen Examples ...}
    } {
	$menu.vmfile${menu2} entryconfig $entry -state disabled	
    }

    # disable all menus
    xcDisableMenuOne $menu [list \
				$menu.vmfile "File" \
				$menu.vmdis  "Display" \
				$menu.vmmod  "Modify" \
				$menu.vmadvg "AdvGeom" \
				$menu.vmpro  "Properties" \
				$menu.vmdat  "Tools"]

    xcDisableOne {.mea.f2.sty .mea.f2.shd}
    
    # when we start the program XCState array is not set
    if { ![array exists XCState] } {
	# enable File menu
	xcEnableMenuOne $menu [list $menu.vmfile File]
	
	# disable all widgets in CONTROL-MAIN FRAME & MEASURE-MAIN FRAME
	xcDisableAll .ctrl .mea

	# in file menu enable new*/open* entries
	foreach entry { 
	    {New CRYSTAL Input} 
	    {Open Structure ...} 
	    {Open PWscf ...} 
	    {Open CRYSTAL ...}
	    {Open WIEN2k ...} 
	    {XCrySDen Examples ...}
	} {
	    $menu.vmfile${menu2} entryconfig $entry -state normal
	}
    }

    if { ![xcIsActive render] } {
	# Enable just File menu
	xcEnableMenuOne $menu [list $menu.vmfile File]

	# disable all widgets in CONTROL-MAIN FRAME & MEASURE-MAIN FRAME
	xcDisableAll .ctrl .mea
    } else {	
	xcEnableMenuOne $menu [list \
				   $menu.vmfile "File" \
				   $menu.vmdis  "Display" \
				   $menu.vmmod  "Modify" \
				   $menu.vmdat  "Tools"]

	# for "render state" enable in File menu:
	foreach entry {
            Close
            {Save XSF Structure} 
	    {Save Current State and Structure}
	    {Save Current State}
            Print*
        } {
	    $menu.vmfile${menu2} entryconfig $entry -state normal
	}
	
	# if structure is CRYSTAL, disable all cell options
	if { $periodic(dim) < 3 } {
	    foreach entry { 
		{Primitive Cell Mode} 
		{Conventional Cell Mode} 				
		{Hexagonal/Rhombohedral ...} 
		{Wigner-Seitz Cells}
	    } {
		$menu.vmdis${menu2} entryconfig $entry -state disabled
	    }
	    $menu.vmmod${menu2} entryconfig {Wigner-Seitz Cell Settings} \
		-state disabled
	} else {
	    foreach entry { 
		{Primitive Cell Mode} 
		{Conventional Cell Mode} 
		{Wigner-Seitz Cells}
	    } {
		$menu.vmdis${menu2} entryconfig $entry -state normal
	    }
	    $menu.vmmod${menu2} entryconfig {Wigner-Seitz Cell Settings} \
		-state normal
	    #if { $periodic(igroup) < 7 } { 
	    #	 $menu.vmdis${menu2} entryconfig {Hexagonal/Rhombohedral ...} \
	    #		 -state disabled
	    #} else {
	    #	 $menu.vmdis${menu2} entryconfig {Hexagonal/Rhombohedral ...} \
	    #		 -state normal	    
	    #}
	    if { [array exists sInfo] } {
		if { !$sInfo(lprimvec) } {
		   #$menu.vmdis${menu2} entryconfig {Hexagonal/Rhombohedral ...} \
		   #	     -state disabled
		    $menu.vmdis${menu2} entryconfig {Primitive Cell Mode} \
			    -state disabled
		    set radio(cellmode) conv
		}
		if { !$sInfo(lconvvec) } {
		    #$menu.vmdis${menu2} entryconfig {Hexagonal/Rhombohedral ...} \
		   # 	    -state disabled		    
		    $menu.vmdis${menu2} entryconfig {Conventional Cell Mode} \
			-state disabled
		    set radio(cellmode) prim
		}
	    }  
	}
	
	if { $periodic(dim) == 0 } {
	    foreach entry { 
		{Crystal Cells} {Wigner-Seitz Cells} 
		{Crystal Cells As ...} {Primitive Cell Mode} 
		{Conventional Cell Mode} {Hexagonal/Rhombohedral ...} 
		{Unit of Repetition ...} 
	    } {
		$menu.vmdis${menu2} entryconfig $entry -state disabled
	    }
	    $menu.vmmod${menu2} entryconfig {Wigner-Seitz Cell Settings} \
		-state disabled
	} else {
	    foreach entry { 
		{Crystal Cells}
		{Crystal Cells As ...}
		{Unit of Repetition ...} 
	    } {
		$menu.vmdis${menu2} entryconfig $entry -state normal
	    }
	    #$menu.vmmod${menu2} entryconfig {Wigner-Seitz Cell Settings} \
	    #	-state normal
	}	
	if { ![info exists xsfAnim(nstep)] } {
	    $menu.vmmod${menu2} entryconfig {Animation Controls} \
		-state disabled
	} else {
	    if { $xsfAnim(nstep) < 1 } {
		$menu.vmmod${menu2} entryconfig {Animation Controls} \
		    -state disabled
	    } else {
		$menu.vmmod${menu2} entryconfig {Animation Controls} \
		    -state normal
	    }
	}
	
	if { $dispmode(style) == "2D" } {
	    xcDisableOne { .mea.f2.sty .mea.f2.shd }
	} else {
	    xcEnableOne { .mea.f2.sty .mea.f2.shd }
	}

        # enable all widgets in CONTROL-MAIN FRAME & MEASURE-MAIN FRAME
	xcEnableAll .ctrl .mea
	# disable PLANE button (temporarily)
	#xcDisableOne .mea.f1.pln

        reloadButtonQueryState
    }
    
    if { [xcIsActive render] && [xcIsActive c95] && ![xcIsActive properties]} {
	xcEnableMenuOne $menu [list $menu.vmadvg "AdvGeom"]
    }
    
    if { [xcIsActive c95] || [xcIsActive newinput] || \
	     [xcIsActive openinput] } {
	# enable the following entries in File menu
	foreach entry { 
	    Close 
	    {Save XSF Structure} 
	    {Save Current State and Structure}
	    {Save Current State}
	    {Save CRYSTAL Input} 
	} {
	    $menu.vmfile${menu2} entryconfig $entry -state normal
	}
    }
    
    if { [xcIsActive properties] && [xcIsActive render] } {
	xcEnableMenuOne $menu [list $menu.vmpro "Properties"]
	foreach entry { 
	    Close 
	    {Save XSF Structure} 
	    {Save Current State and Structure}
	    {Save Current State}
	    {Save CRYSTAL Input} 
	    {Save WIEN2k Struct File} 
	} {
	    $menu.vmfile${menu2} entryconfig $entry -state normal
	}
    }
    
    if { $periodic(dim) == 0 } {
	$menu.vmdat${menu2} entryconfig {k-path Selection} -state disabled
    }
    if { ![array exists sInfo] } {
	$menu.vmdat${menu2} entryconfig {Data Grid} -state disabled
    } elseif { !$sInfo(ldatagrid2D) && !$sInfo(ldatagrid3D) } {
	$menu.vmdat${menu2} entryconfig {Data Grid} -state disabled
    }
    
    if { $dispmode(style) == "2D" } {
	$menu.vmdis${menu2} entryconfig {Wigner-Seitz Cells} -state disabled
    } else {
	if { $periodic(dim) == 3 } {
	    $menu.vmdis${menu2} entryconfig {Wigner-Seitz Cells} -state normal
	}
    }
    
    if { [xcIsActive render] && [xcIsActive c95] && \
	     ([xcIsActive openinput] || [xcIsActive newinput]) && \
	     $periodic(dim) == 2 } {
	foreach entry { 
	    {Atoms & Cell Manipulation ...} 
	    {Cut a SLAB} 
	    {Cut a non-Periodical Structure ...} 
	} {
	    $menu.vmadvg${menu2} entryconfig $entry -state normal
	}
	$menu.vmadvg${menu2}.3 entryconfig {Create a Multi-Slab} -state normal
    } else {
	$menu.vmadvg${menu2}.3 entryconfig {Create a Multi-Slab} -state disabled
    }
    
    if { [xcIsActive multislab] } {
	xcEnableMenuOne $menu [list $menu.vmadvg "AdvGeom"]
	set end [$menu.vmadvg${menu2} index end]
	for {set i 0} {$i <= $end} {incr i} {
	    catch {$menu.vmadvg${menu2} entryconfig $i -state disabled}
	}
	$menu.vmadvg${menu2} entryconfig {Multi-Slab ...} -state normal
	$menu.vmadvg${menu2}.3 entryconfig {Change Multi-Slab Vacuum Thickness} \
	    -state normal
    } else {
	$menu.vmadvg${menu2}.3 entryconfig {Change Multi-Slab Vacuum Thickness} \
	    -state disabled
    }
    
    if { [xcIsActive render] && $periodic(dim) == 3 } {
	$menu.vmfile${menu2} entryconfig {Save WIEN2k Struct File} \
	    -state normal
    }
    
    if { $periodic(dim) == 0 } {
	$menu.vmmod${menu2} entryconfig {Number of Units Drawn} -state disabled
	.mea.f2.unit config -state disabled
	.mea.f2.asym config -state disabled
    } else {
	$menu.vmmod${menu2} entryconfig {Number of Units Drawn} -state normal
	.mea.f2.unit config -state active
	.mea.f2.asym config -state active
    }
    
    if { $system(c95_exist) == 0 } {
	foreach entry { 
	    {New CRYSTAL Input} 
	    {Open CRYSTAL ...} 
	} {
	    $menu.vmfile${menu2} entryconfig $entry -state disabled
	}
    }
        
    global xcMisc
    if { ! [info exists xcMisc(babel)] } {
	catch { $menu.vmfile${menu2}.opstr entryconfig "Gaussian Z-Matrix File" -state disabled }
    }

    # t.k.: this is temporal
    if { ![info exists xcMisc(0.2.x_hexa/rhombo_debug)] } {
	$menu.vmdis${menu2} entryconfig "Hexagonal/Rhombohedral ..." -state disabled
    }

    if { ![info exists xcMisc(movie_encoder)] && ![info exists xcMisc(gif_encoder)]} {
	$menu.vmdat${menu2} entryconfig {Movie Maker} -state disabled
    } else {
	# prevent simultaneous AnimationControl's and MovieMaker's animation capability
	
	global gifAnim
	if { [info exists gifAnim(anim_ctrl_widgets)] } {
	    if { $gifAnim(anim_ctrl_widgets) == 1 } {
		# AnimationControl's "movie options" are On, disable MovieMaker menu
		$menu.vmdat${menu2} entryconfig {Movie Maker} -state disabled
	    } else {
		if { [info exists xcMisc(movie_encoder)] || [info exists xcMisc(gif_encoder)]} {
		    $menu.vmdat${menu2} entryconfig {Movie Maker} -state normal		
		}
	    }
	}

	global movieMaker xsfAnim
	if { [info exists movieMaker(tplw)] } {	   
	    
	    if { [winfo exists $movieMaker(tplw)] } {
		if { [info exists xsfAnim(anim_ctrl_button)] } {
		    if { [winfo exists $xsfAnim(anim_ctrl_button)] } {
			$xsfAnim(anim_ctrl_button) config -state disabled
		    }
		}
	    } else {

		if { [info exists xsfAnim(anim_ctrl_button)] } {
		    if { [winfo exists $xsfAnim(anim_ctrl_button)] } {		    
			$xsfAnim(anim_ctrl_button) config -state normal
		    }
		}
	    }
	}
    }

    # stereo button should be always disabled if stereo is not supported
    global stereo_visual
    if { $stereo_visual == "false" } {
	.mea.f2.stereo configure -state disable
    }
}


############################################################################
# proc append word to XCState variable
proc xcAppendState {word} {
    global XCState

    set app 1
    # maybe XCState do not exists or XCstate(state) is empty
    if ![info exists XCState(state)] { 
	set XCState(state) $word 
	return
    }
    if { $XCState(state) == "" } { 
	set XCState(state) $word 
	return
    }
    
    # if word is already in XCState(state), do not append
    set statelist [split $XCState(state) _]
    foreach item $statelist {
	if { $item == $word } { set app 0 }
    }

    if $app {
	append XCState(state) _$word
    }

    # we will return 1 if appending was successful
    return $app
}


############################################################################
# proc delete word form XCState variable
proc xcDeleteState {word} {
    global XCState

    # if XCState does not exist, there is nothing to delete; return silently
    if ![info exists XCState(state)] { return }

    set statelist [split $XCState(state) _]
    
    set XCState(state) ""
    foreach item $statelist {
	if { $item != $word } {
	    append XCState(state) "${item}_"
	}
    }

    if { $XCState(state) != "" } {
	string trimright $XCState(state) _
    }
}



#############################################################################
# TRACE utility procedure
proc xcTrace {name1 name2 op} {
    global species groupsel XCTrace nxdir nydir nzdir prop \
	    isoControl isosurf radio pDen

    # this is for Print command in File-menu
    if { $name1 == "species" } {
	set name [capitalize $species]
	.menu.vmfile.menu entryconfig Print* -label "Print $name to File"
    } elseif { $name1 == "groupsel" } {
	# in MODIFY/CHANGE toplevel there is 
	# "Type of Cell for Rhombohedral Groups" button, which is activ 
	# only if $groupsel in RHOMBOHEDRAL ONE
	if { [winfo exists $XCTrace(RHOMBO_TYPE_BUTTON)] } {
	    if { [lindex $groupsel 0] == "R" } {
		$XCTrace(RHOMBO_TYPE_BUTTON) config -state normal
	    } else {
		$XCTrace(RHOMBO_TYPE_BUTTON) config -state disabled
	    }
	}
    } elseif { $name1 == "nxdir" } {
	# look in isoControl.tcl file
	if [info exists XCTrace(scX)] {
	    if [winfo exists $XCTrace(scX)] {
		$XCTrace(scX) config -to $nxdir
	    }
	}
	if [info exists XCTrace(2DscX)] {
	    if [winfo exists $XCTrace(2DscX)] {
		$XCTrace(2DscX) config -to $nxdir
	    }
	}
    } elseif { $name1 == "nydir" } {
	# look in isoControl.tcl file
	if [info exists XCTrace(scY)] {
	    if [winfo exists $XCTrace(scY)] {
		$XCTrace(scY) config -to $nydir
	    }
	}
	if [info exists XCTrace(2DscY)] {
	    if [winfo exists $XCTrace(2DscY)] {
		$XCTrace(2DscY) config -to $nydir
	    }
	}
    } elseif { $name1 == "nzdir" } {
	# look in isoControl.tcl file
	if [info exists XCTrace(scZ)] {
	    if [winfo exists $XCTrace(scZ)] {
		$XCTrace(scZ) config -to $nzdir
	    }
	}
	if [info exists XCTrace(2DscZ)] {
	    if [winfo exists $XCTrace(2DscZ)] {
		$XCTrace(2DscZ) config -to $nzdir
	    }
	}
    } elseif { $name1 == "prop" } {
	if { $name2 == "doss_criteria" } {
	    # enable or disable the scale for specifying number
	    if { $prop(doss_criteria) == "band-interval criteria" } {
		DOSS_Init_BandIntv
	    } else {
		DOSS_Init_EnerIntv
	    }
	}
    } elseif { $name1 == "isosurf" } {
	if { $name2 == "3Dinterpl_degree" } {
	    # determine the number of slides for ISO_PLANE123
	    # xc_iso grid returns {nx, ny, nz}
	    if ![info exists isosurf(3Dinterpl_degree_old)] {
		set isosurf(3Dinterpl_degree_old) $isosurf(3Dinterpl_degree)
	    }
	    set n     $isosurf(3Dinterpl_degree)
	    set n_old $isosurf(3Dinterpl_degree_old)
	    xcDebug "XCTrace::::3Dinterpl_degree: n=$n; n_old=$n_old"
	    set r(1) [expr double($n) / double($n_old)]
	    set r(2) [expr double($n) / double($n_old)]
	    set r(3) [expr double($n) / double($n_old)]
	    
	    set isoControl(1,nslide) \
		    [expr ([lindex [xc_iso grid] 2] - 1) * $n + 1] 
	    set isoControl(2,nslide) \
		    [expr ([lindex [xc_iso grid] 1] - 1) * $n + 1]
	    set isoControl(3,nslide) \
		    [expr ([lindex [xc_iso grid] 0] - 1) * $n + 1]
	    xcDebug "XCTrace::::3Dinterpl_degree: 1,2,3,plane: $isoControl(1,nslide) $isoControl(2,nslide) $isoControl(3,nslide)"

	    if ![info exists isoControl(current_slide)] {
		set isoControl(current_slide) 1
	    }
	    foreach i {1 2 3} {
		if ![info exists isoControl($i,current_slide)] {
		    set isoControl($i,current_slide) 1
		} else {
		    if { $r($i) > 1.0 } { set r($i) 1.0 };  # just in any case
		    set isoControl($i,current_slide) [expr round(\
			    ($isoControl($i,current_slide) - 1) * $r($i)) + 1]
		    # just in any case
		    if { $isoControl($i,current_slide) > \
			    $isoControl($i,nslide) } { 
			set isoControl($i,current_slide) \
				$isoControl($i,nslide) 
		    }
		    if { $isoControl($i,current_slide) < 1 } {
			set isoControl($i,current_slide) 1
		    }
		}
		set isoControl($i,current_text_slide) "Current slide:  $isoControl($i,current_slide) / $isoControl($i,nslide)"
		xcDebug "XCTrace::::3Dinterpl_degree: text_slide $isoControl($i,current_text_slide)"
	    }
	}
	set isosurf(3Dinterpl_degree_old) $isosurf(3Dinterpl_degree)
    } elseif { $name1 == "radio" } {
	if { $name2 == ".mesa,bg" } {
            global radio mesa_bg	
	    .mesa.unmap config -bg $radio($name2)
	    set mesa_bg(current) $radio($name2)
	}    
    } elseif { $name1 == "pDen" } {
	for {set i 0} {$i < $pDen(nsurface)} {incr i} {
	    if { $name2 == "$i,monocolor" } {
		$pDen($i,colbut) config -bg [rgb_f2h $pDen($i,monocolor)]
	    }
	}
    } elseif { $name1 == "colSh" } {
	global colSh
	if { $name2 == "slabrange_min" } {
	    if { [winfo exists colSh(scMax)] } {
		$colSh(scMax) config -from $colSh(slabrange_min)
	    }
	} elseif { $name2 == "slabrange_max" } {
	    $colSh(scMin) config -to   $colSh(slabrange_max)
	} elseif { $name2 == "slab_fractional" } {
	    if { $colSh(slab_fractional) } {
		# from absolute --> fractional
		$colSh(scMin) config -from 0
		$colSh(scMax) config -to   1
		ColorScheme:fromAbsToFrac
	    } else {
		# from fractional --> absolute
		$colSh(scMin) config -from $colSh(slab_absrange_min_from)
		$colSh(scMax) config -to   $colSh(slab_absrange_max_to)
		ColorScheme:fromFracToAbs
	    }
	}
    } elseif { $name1 == "toglEps" } {
	global toglEps
	if { $name2 == "EPStype" } {
	    if { $toglEps(EPStype) == "BITMAP" } {
		xcDisableAll $toglEps(frame2)
	    } else { 
		# VECTORIAL
		xcEnableAll $toglEps(frame2)
	    }
	}
    }
}
	

# delete trace on "var" if it exists
proc xcTraceDelete {var} {
    
    set tr [trace vinfo $var]
    if { $tr != "" } { trace vdelete $var }
}


###############################################################################
proc xcAdvGeomState {type {num 1}} {
    global advGeomState AdvGeom

    if { $type == "reset" } {
	if [array exists AdvGeom] { unset AdvGeom }
	set advGeomState(count) 0
    } elseif { $type == "new" } {
	if ![info exists advGeomState(count)] { set advGeomState(count) 0 }
	incr advGeomState(count) $num
	xcDebug "\nxcAdvGeomState:: count == $advGeomState(count)\n"
    } elseif { $type == "undo" } {
	# just decrease count by 1; leave AvdGeom($count,*) alive
	incr advGeomState(count) -$num
    } elseif { $type == "current" } {
	return $advGeomState(count)
    } elseif { $type == "delete" } {
	# decrease count by 1 and unset AdvGeom($count,*)
	for {set i 0} {$i < $num} {incr i} {
	    set alist [array names AdvGeom $advGeomState(count),*]
	    foreach name $alist {
		unset AdvGeom($name)
	    }
	    incr advGeomState(count) -1
	}
    }

    return $advGeomState(count)
}
