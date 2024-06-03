#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/saveState.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

# saveState ?state_and_structure? --> save structure and state
# saveState state --> sane only state
proc saveState {{what state_and_structure}} {
    global saveState

    # get the filename to write the script into

    saveStruct:_queryFilename $what

    # open the channel for writing the script into
    
    set wID [open $saveState(file) w]
    set saveState(wID) $wID

    # write the top info-head
    
    saveStruct:_writeInfoHead $wID $what

    # write the Structure
    if { $what == "state_and_structure" } {
        saveStruct:_writeStructure $wID
    }

    # write the State
    
    saveStruct:_writeState $wID

    # close the channel
    
    flush $wID
    close $wID    
}

proc saveStruct:_writeState {wID} {

    puts $wID {

# ======================================================================== #
#                                                                          #
#                        STATE-PART OF THE FILE                            #
#                                                                          #
# ======================================================================== #

    }

    # --------------------------------------------------
    # definition of xcMisc & printSetup
    # --------------------------------------------------
    global xcMisc
    saveState:header "definition of xcMisc array"
    puts $wID "array set xcMisc [list [array get xcMisc]]"

    global printSetup
    if { [array exists printSetup] } {
	saveState:header "definition of printSetup"     
	puts $wID "array set printSetup [list [array get printSetup]]"
    }

    # --------------------------------------------------
    # display watch-cursor
    # --------------------------------------------------

    saveState:header "display \"waiting\" toplevel and watch-cursor"
    puts $wID "set wait_window \[DisplayUpdateWidget \"Reconstructing\" \"Reconstructing the structure and display parameters. Please wait\"\]"
    puts $wID "SetWatchCursor"
    puts $wID "set xcCursor(dont_update) 1"

    # --------------------------------------------------
    # remember the size of the main window
    # --------------------------------------------------

    saveState:header "size of the main window"
    set w [winfo width  .]
    set h [winfo height .]
    puts $wID "wm geometry . ${w}x${h}"

    # --------------------------------------------------
    # create all fonts
    # --------------------------------------------------

    saveState:beginHeader "create fonts"
    foreach font [font names] {
	puts $wID "saveState:fontCreate $font [font actual $font]"
    }    
    saveState:endHeader

    # --------------------------------------------------
    # display-mode
    # --------------------------------------------------

    global light mode2D mode3D dispmode style3D viewer translationStep rotstep
    saveState:beginHeader "take care of display-mode"
    puts $wID "set translationStep $translationStep"
    puts $wID "set rotstep $rotstep"
    puts $wID "set light $light" 
    puts $wID "Lighting $light"
    puts $wID "array set mode2D   [list [array get mode2D]]"
    puts $wID "array set mode3D   [list [array get mode3D]]"
    puts $wID "array set dispmode [list [array get dispmode]]"
    puts $wID "saveState:displayMode"

    puts $wID "set style3D(draw)  $style3D(draw); Style3D draw $style3D(draw)"
    puts $wID "set style3D(shade) $style3D(shade); Style3D shade $style3D(shade)"
    puts $wID "set viewer(rot_zoom_button_mode) $viewer(rot_zoom_button_mode); Viewer:rotZoomButtonMode"
    saveState:endHeader

    # ------------------------------------------------------------------------
    # Take care of MODIFY menu parameters
    # ------------------------------------------------------------------------

    # N.B: loading of script is very slow if (nxdir,nydir,nzdir) are
    # specified that soon, consider shifting just befor the xc_rotationmatrix

    # t.k. TESTING: let's how it works if this is shifted down !!!
    
    # # This item must be the first one:
    # # Number of Units Drawn (must be after DISPLAY-menu checkbuttons)
    # global nxdir nydir nzdir
    # saveState:header "Number of Units Drawn"
    # puts $wID "set nxdir $nxdir"
    # puts $wID "set nydir $nydir"
    # puts $wID "set nzdir $nzdir"
  
    saveState:beginHeader "Various colors"

    # Atomic-Color
    global mody    
    saveState:_foreachAtomicType $mody(L_ATCOL_ONE) $mody(D_ATCOL_ONE)

    # Unibond-Color
    puts $wID "xc_newvalue .mesa $mody(L_UNIBONDCOLOR) [xc_getvalue $mody(D_UNIBONDCOLOR)] 1.0"

    # Crystal Cell's Color
    global cellcol
    puts $wID "xc_newvalue .mesa $mody(L_FRAMECOL) [xc_getvalue $mody(D_FRAMECOL)]"
    saveState:endHeader

    # balls/spacefill based on covalent or van der Waals ?

    global radio
    #if { $radio(ball) == "Balls based on covalent radii" } {
    #	puts $wID "xc_newvalue .mesa $mody(L_BALL_COV)"
    #} elseif { $radio(ball) == "Balls based on Van der Waals radii" } {
    #	puts $wID "xc_newvalue .mesa $mody(L_BALL_VDW)"
    #}
    if { $radio(space) == "SpaceFill based on covalent radii" } {
	puts $wID "xc_newvalue .mesa $mody(L_SPACE_COV)"
    } elseif { $radio(space) == "SpaceFill based on Van der Waals radii" } {
	puts $wID "xc_newvalue .mesa $mody(L_SPACE_VDW)"
    }

    # chemical connectivity factor, ball-factor, ball/stick ratio, point-radius, ...
    saveState:header "Various parameters"
    # on Thu 17 Aug, this has been skipped from below loop L_ATRAD_SCALE      D_ATRAD_SCALE
    foreach {newValueID getValueID} {
	L_COV_SCALE        D_COVF
	L_BALLF            D_BALLF
	L_RODF             D_RODF
	L_WFLINEWIDTH      D_WFLINEWIDTH
	L_PLLINEWIDTH      D_PLLINEWIDTH
	L_FRAMELINEWIDTH   D_FRAMELINEWIDTH
	L_OUTLINEWIDTH     D_OUTLINEWIDTH
	L_WF3DLINEWIDTH    D_WF3DLINEWIDTH
	L_PLRADIUS         D_PLRADIUS
	L_FRAMERODF        D_FRAMERODF
	L_TESSELLATION     D_TESSELLATION
	L_PERSPECTIVEBACK   D_PERSPECTIVEBACK 
	L_PERSPECTIVEFOVY  D_PERSPECTIVEFOVY
	L_PERSPECTIVEFRONT  D_PERSPECTIVEFRONT
	L_BACKGROUND       D_BACKGROUND
    } {
	puts $wID "xc_newvalue .mesa $mody($newValueID) [xc_getvalue $mody($getValueID)]"
    }

    # Lights
    global glLight
    if { [info exists glLight(nlights)] } {
	saveState:beginHeader "OpenGL Lights"
	puts $wID "array set glLight [list [array get glLight]]"
	for {set i 0} {$i < $glLight(nlights)} {incr i} {
	    puts $wID "glLight:update .mesa $i"
	}
	saveState:endHeader
    }

    # Materials/Fog/Antialias
    global glModParam
    if { [array exists glModParam] } {
	saveState:beginHeader "Materials/Fog/Antialias"
	puts $wID "array set glModParam [list [array get glModParam]]"
	puts $wID "glModParam:Materials:Update  .mesa"
	puts $wID "glModParam:DepthCuing:Update .mesa"
	puts $wID "glModParam:AntiAlias:Update  .mesa"
	saveState:endHeader
    }

    # Forces ...
    global forceVec
    if { [array exists forceVec] } {
	saveState:beginHeader "Forces"
	puts $wID "array set forceVec [list [array get forceVec]]"
	puts $wID "xc_forces .mesa scalefunction  $forceVec(scalefunction)"
	puts $wID "xc_forces .mesa threshold      $forceVec(threshold)"
	puts $wID "xc_forces .mesa lengthfactor   $forceVec(lengthfactor)"
	puts $wID "forceUpdate .mesa"
	saveState:endHeader
    }

    # H-bonds
    global Hbonds
    if { [array exists Hbonds] } {
	if { [info exists Hbonds(colorID)] } {
	    set _hbonds_color_id $Hbonds(colorID)
	    unset Hbonds(colorID)
	}

	saveState:beginHeader "H-bonds"
	puts $wID "array set Hbonds [list [array get Hbonds]]"
	puts $wID "HbondsSetting:update .mesa"	
	saveState:endHeader

	if { [info exists _hbonds_color_id] } {
	    set Hbonds(colorID) $_hbonds_color_id
	}
    }

    # Molecular Surfaces ...
    global pDen
    if { [info exists pDen(nsurface)] } {
	saveState:beginHeader "Molecular-Surface"
    	puts $wID "array set pDen [list [array get pDen]]"
	for {set i 0} {$i < $pDen(nsurface)} {incr i} {
	    set pDen($i,type_old) unknown
	    puts $wID "set pDen($i,ident) \[xc_molsurfreg .mesa\]"
	    puts $wID "PseudoDen:_xc_molsurf $i"
	}
	saveState:endHeader
    }

    # ------------------------------------------------------------------------
    # DISPLAY-menu checkbuttons
    # ------------------------------------------------------------------------

    global check
    saveState:header "Various displays (i.e. checkbuttons of DISPLAY menu)"
    puts $wID "array set check [list [array get check]]"
    foreach cmd {
	CrdSist AtomLabels CrysFrames Unibond {forceVectors .mesa} {Hbonds .mesa}
	WignerSeitz Perspective
    } {
	puts $wID $cmd
    }

    # ------------------------------------------------------------------------
    # DISPLAY-menu radiobuttons
    # ------------------------------------------------------------------------

    global periodic
    saveState:header "Various displays (i.e. radiobuttons DISPLAY menu)"
    puts $wID "array set radio [list [array get radio]]"
    if { $periodic(dim) > 0 } {
	foreach cmd {
	    {CellMode 1}
	    DispFramesAs
	} {
	    puts $wID $cmd
	}
    }

    global radio
    if { $radio(ball) == "Balls based on covalent radii" } {
	puts $wID "xc_newvalue .mesa $mody(L_BALL_COV)"
    } elseif { $radio(ball) == "Balls based on Van der Waals radii" } {
	puts $wID "xc_newvalue .mesa $mody(L_BALL_VDW)"
    }

    # ------------------------------------------------------------------------
    # This is from MODIFY menu, but should presumably be here ?!?!?
    #
    # Wigner-Seitz Cell (BEWARE: for Wigner-Setz we should know the
    # cell-mode [prim,conv])
    # ------------------------------------------------------------------------

    global ws wsp wsc ws_lfpos ws_npos
    if { $periodic(dim) > 2 } {
	saveState:beginHeader "Wigner-Seitz Cell"
	puts $wID "set check(wigner) 1"
	puts $wID "WignerSeitz"
	foreach array {ws wsp wsc ws_lfpos ws_npos} {
	    if { [array exists $array] } {
		puts $wID "array set $array [list [array get $array]]"
	    }
	}
	puts $wID "SetWignerSeitzInit; update; SetWignerSeitz_OK test; SetWignerSeitz_Cancel .wgnset"
	saveState:endHeader
    }

    # t.k: testing

    # Atomic Radii
    saveState:beginHeader "Atomic radii"
    saveState:_foreachAtomicType $mody(L_ATRAD_ONE) $mody(D_ATRAD_ONE)
    saveState:_foreachAtomicType $mody(L_RCOV_ONE)  $mody(D_RCOV_ONE)
    saveState:endHeader

    # t.k./

    # ------------------------------------------------------------------------
    # TOOLS-menu
    # ------------------------------------------------------------------------

    xcDebug -stderr "saving colSh ..."
    # color scheme
    global colSh
    if { [array exists colSh] } {
	saveState:header "load appropriate color-scheme"
	puts $wID "array set colSh [list [array get colSh]]" 
	puts $wID "ColorSchemeUpdate .mesa"
    }

    # --------------------------------------------------
    # take care of 2D|3D datagrid
    # --------------------------------------------------

    xcDebug -stderr "saving datagrid ..."
    global sInfo DataGrid DG isoControl prop isosurf
    if { ! [info exists prop(datagridDim)] } {
	set prop(datagridDim) 0
    }
    if { $sInfo(ldatagrid2D) || $sInfo(ldatagrid3D) \
	     || $prop(datagridDim) == 2 || $prop(datagridDim) == 3 } {
	if { [info exists DataGrid(launch_command)] } {
	    # the scalar field is rendered	    
	    saveState:beginHeader "scalar field (controur or isosurface) settings"
	    puts $wID "array set DG [list [array get DG]]"
	    puts $wID "DataGridOK"
	    foreach array {isoControl prop isosurf} {
		puts $wID "array set $array [list [array get $array]]"
	    }
	    if { $DataGrid(dim) == "2D" } {
		puts $wID "UpdatePropertyPlane"
	    } else {
		puts $wID "UpdateIsosurf"
	    }
	    if { [info exists DataGrid(first_time)] } {
		puts $wID "array set DataGrid [list [array get DataGrid]]"
	    }
	    saveState:endHeader
	}
    }

    # Isosurfaces ...
    global openGL
    if { [array exists openGL] } {
	saveState:beginHeader "Isosurface colors/transparency ..."
	
	# colors
	foreach what {
	    {isosurf pos front}
	    {isosurf pos back}
	    {isosurf neg front}
	    {isosurf neg back}
	    {isosurf one front}
	    {isosurf one back}
	} {
	    set com [concat xc_getGLparam material -what $what]
	    set ambient   [eval $com -get ambient]
	    set diffuse   [eval $com -get diffuse]
	    set specular  [eval $com -get specular]
	    set emission  [eval $com -get emission]
	    set shininess [eval $com -get shininess]
	    puts $wID "xc_setGLparam material -what $what \\"
	    puts $wID "   -shininess [list $shininess] \\"
	    puts $wID "   -specular  [list $specular]  \\"
	    puts $wID "   -ambient   [list $ambient]  \\"
	    puts $wID "   -diffuse   [list $diffuse]  \\"
	    puts $wID "   -emission  [list $emission]"
	}
	
	# blendfunc
	set bldf      [xc_getGLparam blendfunc -what isosurf]
	set src_blend [lindex $bldf 0]
	set dst_blend [lindex $bldf 1]
	puts $wID "xc_setGLparam blendfunc -what isosurf -sfunc $src_blend -dfunc $dst_blend"	

	saveState:endHeader
    }

    # t.k. EXPERIMENTAL: testing for faster speed !!!
    global nxdir nydir nzdir
    saveState:header "Number of Units Drawn"
    puts $wID "set nxdir $nxdir"
    puts $wID "set nydir $nydir"
    puts $wID "set nzdir $nzdir"

    # --------------------------------------------------
    # PseudoDensity & Anti-Aliasing & Depth Cuing (they are the most
    # time consuming, so we do it at end
    # --------------------------------------------------

    xcDebug -stderr "saving anti-alias ..."

    saveState:header "Anti-Aliasing & Depth-Cuing & PseudoDensity (these are time consuming)"
    puts $wID "DepthCuing; PseudoDensity; AntiAlias"    

    # update the display
    if { $periodic(dim) > 0 } {
	saveState:header "update the display"
	puts $wID "CellMode 1"
    }

    #------------------------------------------------------------------------
    # Atomic-Labels/Fonts
    #------------------------------------------------------------------------
    global atomLabel
    saveState:beginHeader "Atomic-Labels/Fonts"
    puts $wID "array set atomLabel [list [array get atomLabel]]"
    puts $wID "set t \[ModAtomLabels\]"
    puts $wID "xc_setfont .mesa \
        [list $atomLabel(globalFont)] \
        [list $atomLabel(globalFont.brightColor)] \
        [list $atomLabel(globalFont.darkColor)]"
    foreach name [array names atomLabel *.atomFont.label] {
	#regsub ?switches? exp string subSpec varName
	regsub -- {\.atomFont\.label$} $name {} id
	
	puts $wID "xc_setatomlabel .mesa $id \
          [list $atomLabel($id.atomFont.label)] \
          [list $atomLabel($id.atomFont.font)] \
          [list $atomLabel($id.atomFont.brightColor)] \
          [list $atomLabel($id.atomFont.darkColor)]"
    }
    puts $wID "ModAtomLabels:advancedCheckButton default"
    puts $wID "ModAtomLabels:advancedCheckButton custom"
    puts $wID "ModAtomLabels:advancedCloseUpdate dummy update"
    puts $wID "CancelProc \$t"
    saveState:endHeader

    # --------------------------------------------------
    # rotation matrix and zooming factor, and translation displacements
    # --------------------------------------------------

    xcDebug -stderr "saving rotations/translations ..."

    saveState:header "rotation matrix and zooming factor, and translation displacements"
    puts $wID "xc_rotationmatrix set [xc_rotationmatrix get]"
    puts $wID "xc_translparam    set [xc_translparam get]"
    puts $wID ""
    puts $wID "# this is used to force the update of display"
    puts $wID "cry_toglzoom .mesa 0.0"

    # --------------------------------------------------
    # reset cursor and delete the "waiting" window
    # --------------------------------------------------
    saveState:header "reset cursor"
    puts $wID "set xcCursor(dont_update) 0"
    puts $wID "ResetCursor"
    puts $wID "destroy \$wait_window"
}



proc saveStruct:_queryFilename {what} {
    global system saveState

    if { ![xcIsActive render] } { return -code return }

    set filetypes {
	{{XCrySDen Files}         {.xcrysden}    }
	{{Gzipped XCrySDen Files} {.xcrysden.gz} }
	{{All Files}              *              }
    }

    if { $what == "state" } {
        set title "Save Current State"
    } else {
        set title "Save Current State and Structure"
    }
    set saveState(file) [tk_getSaveFile -initialdir $system(PWD) \
			     -title $title \
			     -defaultextension .xcrysden \
			     -filetypes $filetypes]
    
    # maybe Cancel button was pressed
    
    if { $saveState(file) == "" } { 	
	return -code return
    }

    # open the channel for writing the script into
    set saveState(wID) [open $saveState(file) w]

    return $saveState(wID)
}

proc saveStruct:_writeStructure {wID} {
    global system mody
    
    # --------------------------------------------------
    # non-XSF to XSF transformations (i.e. Take care of PDB & XYZ files) 
    # --------------------------------------------------
    
    set FORMAT_NONE  0
    set FORMAT_XSF   1
    set FORMAT_XYZ   2
    set FORMAT_PDB   3
    set format [xc_getvalue $mody(D_CURRENTFILEFORMAT)]
    if { $format == $FORMAT_NONE } {
	return -code return
    }

    set xsf_tmp_file [xcTempFile xc_xsf]
    xc_writeXSF $xsf_tmp_file        
    
    # ------------------------------
    # Take care of CRYSTAL98 Properties
    # ------------------------------

    if { [xcIsActive c95] && [xcIsActive properties] } {
	global prop
	# we are dealing with CRYSTAL98 properties. If we are
	# rendering the isosurface/contours, then this requires some
	# special treatment.
	if { $prop(datagridDim) == 2 || $prop(datagridDim) == 3 } {
	    saveState:_saveDatagrid $xsf_tmp_file
	} 
    }

    # --------------------------------------------------
    # the XSF file is stored in $xsf_tmp_file, write the XSF part of
    # *.xcrysden script
    # --------------------------------------------------

    puts $wID {
        
# ======================================================================== #
#                                                                          #
#                      STRUCTURE-PART OF THE FILE                          #
#                                                                          #
# ======================================================================== #

    }

    saveState:beginHeader "XSF structure data"    
    set xsfFile "\$system(SCRDIR)/xc_xsf.\$system(PID)"

    puts $wID "# Store the content of the XSF file in the xsfStructure variable"
    puts $wID ""
    puts $wID "set xsfStructure \{"
    set rID [open $xsf_tmp_file r]
    fcopy $rID $wID
    close $rID

    # check for scalar field
    global sInfo prop DataGrid

    if { ! [info exists prop(datagridDim)] } {
	set prop(datagridDim) 0
    }
	
    if { ! [xcIsActive c95] } {
	# NOTE: for CRYSTALxx program the datagrids are handled separately

	if { $sInfo(ldatagrid2D) || $sInfo(ldatagrid3D) \
		 || $prop(datagridDim) == 2 || $prop(datagridDim) == 3 } {
	    if { [info exists DataGrid(launch_command)] } {
		
		# the scalar field is rendered	    
		
		set sfile  [xcTempFile xc_tmp]
		set ident1 "UNKNOWN"
		set ident2 "DataGrid_Generated_by_XCrySDen"
		set dg     "DATAGRID${DataGrid(dim)}D"
		
		xc_iso save $sfile $ident1
		
		puts $wID "BEGIN_BLOCK_$dg"
		puts $wID "$ident2"
		puts $wID [ReadFile -nonewline $sfile]
		puts $wID "\nEND_BLOCK_$dg"	    
	    }
	}
    }

    puts $wID "\}"
    puts $wID "# END of defintion of xsfStructure variable"
    puts $wID ""
    puts $wID "WriteFile $xsfFile \$xsfStructure w"
    puts $wID ""
    puts $wID "::scripting::open --xsf $xsfFile"

    saveState:endHeader
   
    return $wID
}


proc saveState:_foreachAtomicType {newValueID getValueID} {
    global mody saveState
    foreach atom [AtomNames] {
	set nat [Aname2Nat $atom]
	puts $saveState(wID) \
	    "xc_newvalue .mesa $newValueID $nat [xc_getvalue $getValueID $nat]"
    }
}


proc saveState:header {text} {
    global saveState

    set wID $saveState(wID)
    puts $wID "\n\# ------------------------------------------------------------------------\n\# $text\n\# ------------------------------------------------------------------------\n"
    set saveState(headerText) $text
}

proc saveState:beginHeader {text} {
    global saveState

    set wID $saveState(wID)
    puts $wID "\n\# ------------------------------------------------------------------------\n\# BEGIN: $text\n\# ------------------------------------------------------------------------\n"
    set saveState(headerText) $text
}
proc saveState:endHeader {} {
    global saveState

    set wID $saveState(wID)
    puts $wID "\n\# ------------------------------------------------------------------------\n\# END: $saveState(headerText)\n\# ------------------------------------------------------------------------\n\n"
}

proc saveState:displayMode {} {
    global light mode2D mode3D dispmode

    if { $light == "On" } {
	if { $dispmode(mode3D) == "Preset" } {
	    DisplayOver3D $dispmode(mode3D_name)
	} else {
	    # $dispmode(mode3D) == "Logic" --> revert to preset mode
	    Mode3D Preset	    
	}
    } else {
	Display2D $dispmode(mode2D_name)
    }
}


proc saveState:fontCreate {font args} {
    set fonts [font names]
    if { [lsearch $fonts $font] < 0 } {
	# font does not yet exists
	eval font create $font $args
    }
}



proc saveState:_saveDatagrid {xsf_tmp_file} {
    global xcMisc system isoControl DataGrid DG
    
    set gridFile [xcTempFile xc_xsf_datagrid]

    set ident UNKNOWN
    xc_iso save $gridFile $ident
    
    set ident DataGrid_Generated_By_XCrySDen
    
    # set DG correspondingly
    set DG(radio)             0
    set DG(n_subblock,0)      1
    set DG(cb0,0)             1
    set DG(envar0,0)          1.0
    set DataGrid(first_time)  exists
    # replace this ...
    if { $isoControl(datagridDim) == 3 } {
	# 3D
	set dg                       DATAGRID3D
	set DataGrid(launch_command) IsoControl
	set DataGrid(dim)            3
	set DG(type,0)               3D
    } else {
	#2D
	set dg                       DATAGRID2D
	set DataGrid(launch_command) IsoControl2D
	set DataGrid(dim)            2
	set DG(type,0)               2D
     }
    
    # take care about the newlines
    set content   [ReadFile -nonewline $xsf_tmp_file]
    set datablock 0
    set out       {}
    foreach line [split $content \n] {
	if {[string match BEGIN_BLOCK_* $line]} {
	    set datablock 1
	}
	if {!$datablock} {
	    append out "$line\n"
	}
	if {[string match END_BLOCK_* $line]} {
	    set datablock 0
	}
    }
    append out "BEGIN_BLOCK_$dg\n"
    append out "$ident\n"
    append out [ReadFile -nonewline $gridFile]
    append out "\nEND_BLOCK_$dg"
    WriteFile $xsf_tmp_file $out w
}

proc saveStruct:_writeInfoHead {wID what} {
    global saveState env

    set date     "Xxx Xxx XX XX:XX:XX CEST 20XX"
    set user     "unknown"
    set hostname "localhost"
    set domain   "local"
    catch {set date     [clock format [clock seconds]]}
    catch {set user     $env(USER)}
    catch {set hostname $env(HOSTNAME)}
    catch {set domain   [exec dnsdomainname]}

    if { $what != "state" } {
        puts $wID [subst {
# ------------------------------------------------------------------------
# This is the script saved via the "File-->Save Current State and
# Structure" menu. This feature saves the currently displayed
# structure. When it will be loaded again it will recall all the
# display parameters, such as orientation, zoom, perspective, ...
#
# Execute this script as: xcrysden --script $saveState(file)
#
# File created by $user on $hostname.$domain
# File creation date: $date    
# ------------------------------------------------------------------------
        }]
    } else {
        puts $wID [subst {
# ------------------------------------------------------------------------
# This is a script saved via the "File-->Save Current State" menu. It
# is useful for the user written scripts as to get the desired "state"
# of display, but as standalone this script is useless.
#
# File created by $user on $hostname.$domain
# File creation date: $date    
# ------------------------------------------------------------------------
        }]
    }
}

