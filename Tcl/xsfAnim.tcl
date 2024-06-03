#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xsfAnim.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#
# tis proc is called form menu File->...->(Open XSF File)
#
proc xsfAnimOpenMenu {can} {
    global fileselect system

    set fileselect(path) [tk_getOpenFile -defaultextension .xsf \
                              -filetypes { 
                                  {{All Files}                  {.*}  }
                                  {{Animated XSF Files}         {.axsf}}
                                  {{XSF Files}                  {.xsf}}
                                  {{GZipped Animated XSF Files} {.axsf.gz}}
                                  {{GZipped XSF Files}          {.xsf.gz}}
                              } -initialdir $system(PWD) -title "Open AXSF File"]
    if { $fileselect(path) == "" } {
        return
    }
    reloadRegister $fileselect(path) XSF
    xsfOpen $fileselect(path) $can
    #xsfAnimInit $fileselect(path) $can
}


proc xsfAnimInit {filedir can {update 0}} {
    global system geng xcMisc periodic radio sInfo periodic xsfAnim

    set FileName [file tail $filedir]

    ##################################
    # everything must happen in SCRDIR
    cd $system(SCRDIR)

    set xsfAnim(nstep) [xsfParseAnim $filedir]
    if { $xsfAnim(nstep) < 0 } {
	# an Error has occured
	tk_dialog [WidgetName] "ERROR" "ERROR: An error occured when parsing AXSF file $FileName !!!" error 0 OK
	return 0
    } elseif { $xsfAnim(nstep) == 0 } {
	tk_dialog [WidgetName] "NOTIFICATION" "NOTIFICATION: Specified file $FileName seems to be an ordinary XSF file !!!" warning 0 OK
	set xsfAnim(not_anim) 1
	xsfOpen $filedir $can
	return 1
    } 
    
    set xsfAnim(anim_step) 1
    set xsfAnim(current)   1
    set xsfAnim(label_text_current) "Current slide: $xsfAnim(current)/$xsfAnim(nstep)"
    
    xsfOpenParsed $system(SCRDIR)/$FileName.1 $can $update

    set xsfAnim(filename) $FileName
    xsfAnimWid $can
}


proc xsfAnimWid can {
    global xsfAnim gifAnim xcMisc
    
    if { ! [info exists xsfAnim(nstep)] } {
	return
    }
    if { $xsfAnim(nstep) < 1 } {
	return
    }
    
    # set the daley between slides (in msec)
    set xsfAnim(delay) 50
    set c [lindex [split $can .] end]
    if { [winfo exists .anim$c] } {
	return
    }

    set xsfAnim(label_text_current) "Current slide: $xsfAnim(current)/$xsfAnim(nstep)"
    set gifAnim(anim_ctrl_widgets)  0

    # from here on is just the widget managing
    set t [xcToplevel .anim$c "Animation Control" "Animation" . -0 0 1]
    set xsfAnim(tplw) $t
    bind $t <Destroy> {xsfAnimWidClose %W}

    # this is for "Hide" button
    global unmapWin
    xcRegisterUnmapWindow $t $unmapWin(frame,main) anim \
	    -textimage {Animate unmap}
    bind $t <Unmap> [list xcUnmapWindow unmap %W $t \
	    $unmapWin(frame,main) anim]
    bind $t <Map>   [list xcUnmapWindow map %W $t \
	    $unmapWin(frame,main) anim]

    # container frames
    set f1 [frame $t.f1 -relief ridge -bd 2]
    set f2 [frame $t.f2 -relief ridge -bd 2]
    set f3 [frame $t.f3]
    pack $f1 -side top -expand 1 -fill x \
	    -padx 2m -pady 2m -ipadx 2m -ipady 2m 
    pack $f2 -side top -expand 1 -fill x \
	    -padx 2m -ipadx 2m -ipady 1m 
    pack $f3 -side top


    # --------------------------------------------------
    # xsfAnim widgets
    # --------------------------------------------------

    set e1 [OneEntries $f1 \
	    {"Delay between slides (in msec):" "Animation step:"} \
	    {xsfAnim(delay) xsfAnim(anim_step)} 31 15 5 -side top]    
    set slide_label [label $f1.l1 -textvariable xsfAnim(label_text_current) -anchor c]
    set lfont [ModifyFont [$slide_label cget -font] $slide_label -size 18 -underline 1]
    $slide_label configure -font $lfont
    set f_play [frame $f1.1]

    set first [button $f_play.1st -image first    -anchor center \
	    -command [list xsfAnimAction $can first]]
    set backw [button $f_play.bck -image backward -anchor center \
	    -command [list xsfAnimAction $can backward]]
    set previ [button $f_play.prv -image previous -anchor center \
	    -command [list xsfAnimAction $can previous]]
    set stop [button $f_play.sto -image stop -anchor center \
	    -command [list xsfAnimAction $can stop]]
    set next  [button $f_play.nxt -image next     -anchor center \
	    -command [list xsfAnimAction $can next]]
    set forw  [button $f_play.frw -image forward  -anchor center \
	    -command [list xsfAnimAction $can forward]]
    set last  [button $f_play.lst -image last     -anchor center \
	    -command [list xsfAnimAction $can last]]

    pack $slide_label $f_play -side top -expand 1 -pady 1m
    pack $first $backw $previ $stop $next $forw $last -side left

    foreach {wid text} {
	1st  "First image"
	bck  "Play backward"
	prv  "Previous image"
	sto  "Stop playing"
	nxt  "Next image"
	frw  "Play forward"
	lst  "Last image"
    } {
	set path $f_play.$wid		
	DynamicHelp::register $path balloon $text
    }

    # --------------------------------------------------
    # gifAnim widgets
    # --------------------------------------------------
    
    set gifAnim(make_gifAnim) 0
    
    set gifAnim(button_text) "Animated GIF/MPEG/AVI >>"
    set f_animgif [frame $f2.1]
    set gifb [button $f2.gifb \
		  -textvariable gifAnim(button_text) \
		  -command [list gifAnimWidPack $f_animgif $can]]
    pack $gifb -side top -padx 10 -pady 3 -fill x -expand 1

    set xsfAnim(anim_ctrl_button) $gifb

    if { ![info exists xcMisc(gif_encoder)] && ![info exists xcMisc(movie_encoder)] } {
	$gifb config -state disabled
    }

    # make a movie-control windgets
    gifAnim_controlWidgets $f_animgif $can $slide_label

    
    set gifAnim(make_text)    "Start Recording Animation"
    set gifAnim(make_gifAnim) 0
    set makeGIF [button $f_animgif.make -textvariable gifAnim(make_text) \
		     -command [list gifAnimMake $can $slide_label $f_animgif.make]]
    pack $makeGIF -side top -padx 5 -pady 5

    #
    # CLOSE/HIDE widgets
    #
    set hide  [button $f3.hide  -text "Hide" \
	    -default active -command [list HideWin $t anim]]
    set close [button $f3.close -text "Close" \
	    -command [list destroy $t]]
    pack $hide $close -side left -ipadx 1m -ipady 1m -padx 2m -pady 2m

    xcUpdateState
}

proc xsfAnimWidClose {t} {
    global gifAnim
    if { [winfo exists $t] } {
	destroy $t
    }
    set gifAnim(anim_ctrl_widgets) 0
    xcUpdateState
}

proc xsfAnimAction {can what} {
    global xsfAnim
    
    if { ![info exists xsfAnim(stop)] } {
	set xsfAnim(stop) 0
    }
    if { ![info exists xsfAnim(anim_step)] } {
	set xsfAnim(anim_step) 1
    }
    if { $xsfAnim(anim_step) < 1 } {
	set xsfAnim(anim_step) 1
    }

    if { ![number xsfAnim(delay)] } {
	tk_dialog [WidgetName] "ERROR" "ERROR: Please specify a number for delay !!!" error 0 OK
	return 0
    }

    switch -exact -- $what {
	stop   {
	    set xsfAnim(stop) [expr {$xsfAnim(stop) ? 0 : 1}]
	}
	first   {		   
	    if { $xsfAnim(current) == 1 } {
		return
	    }
	    set xsfAnim(current) 1	
	    xsfAnimRenderCurrent $can
	}
	backward {
	    while { $xsfAnim(current) > 1 && ! $xsfAnim(stop)} {
		incr xsfAnim(current) -$xsfAnim(anim_step)
		if { $xsfAnim(current) < 1 } {
		    set xsfAnim(current) 1
		}
		xsfAnimRenderCurrent $can
		update
		set xsfAnim(label_text_current) \
			"Current slide: $xsfAnim(current)/$xsfAnim(nstep)"
		after $xsfAnim(delay)
	    }
	    set xsfAnim(stop) 0
	    return
	}
	previous {
	    if { $xsfAnim(current) > 1 } {
		incr xsfAnim(current) -$xsfAnim(anim_step)
		if { $xsfAnim(current) < 1 } {
		    set xsfAnim(current) 1
		}
		xsfAnimRenderCurrent $can
	    } else {
		return
	    }
	}
	next     {
	    if { $xsfAnim(current) < $xsfAnim(nstep) } {
		incr xsfAnim(current) $xsfAnim(anim_step)
		if { $xsfAnim(current) > $xsfAnim(nstep) } {
		    set xsfAnim(current) $xsfAnim(nstep)
		}
		xsfAnimRenderCurrent $can
	    } else {
		return
	    }
	}
	forward  {
	    while { $xsfAnim(current) < $xsfAnim(nstep) && ! $xsfAnim(stop) } {
		incr xsfAnim(current) $xsfAnim(anim_step)
		if { $xsfAnim(current) > $xsfAnim(nstep) } {
		    set xsfAnim(current) $xsfAnim(nstep)
		}
		xsfAnimRenderCurrent $can
		update
		set xsfAnim(label_text_current) \
			"Current slide: $xsfAnim(current)/$xsfAnim(nstep)"
		after $xsfAnim(delay)
	    }
	    set xsfAnim(stop) 0
	    return
	}
	last     {
	    if { $xsfAnim(current) == $xsfAnim(nstep) } {
		return
	    }
	    set xsfAnim(current) $xsfAnim(nstep)
	    xsfAnimRenderCurrent $can
	}
    }
    set xsfAnim(label_text_current) "Current slide: $xsfAnim(current)/$xsfAnim(nstep)"    
}

proc xsfAnimRenderCurrent can {
    global xsfAnim gifAnim periodic nxdir nydir nzdir geng system \
	working_XSF_file

    SetWatchCursor
    set working_XSF_file $xsfAnim(filename).$xsfAnim(current)

    if { $periodic(dim) > 0 } {
	# new: tk Mon Jan 27 15:27:23 CET 2003
	CellMode 1
	#GenGeom $geng(M1_PRIM) $geng(M2_CELL) 11 1 $nxdir $nydir $nzdir \
	#	xc_struc.$system(PID)	
	#UpdateStruct $can xc_struc.$system(PID)
    } else {
	UpdateStruct $can $system(SCRDIR)/$working_XSF_file
    }
    ResetCursor    

    if { $gifAnim(make_gifAnim) == 1 } {
	# 
	# Animated-GIF option:  print this snapshot to GIF file
	gifAnimPrintCurrent $can
    }
    return
}


# sample AXSF file (fixed-cell):
#------------------
# ANIMSTEPS 18
# DIM-GROUP
# 3 1
# PRIMVEC
#   5.8859563944     0.0000000000     0.0000000000
#   0.0000000000     5.8859563944     0.0000000000
#   0.0000000000     0.0000000000    23.3615963678
# PRIMCOORD 1
#32 1
# 47    3.8888048966     4.1163007841    19.5974532965
# 47    3.8888048966     1.7696556103    19.5974532965
# 47    1.9971514978     4.1163007841    19.5974532965
# 47    1.9971514978     1.7696556103    19.5974532965
#  6    2.9429781972     3.5961067857    19.2368415845
#  6    2.9429781972     2.2898496088    19.2368415845
#  ..............

# sample AXSF file (variable-cell):
#------------------
# ANIMSTEPS 18
# DIM-GROUP
# 3 1
# PRIMVEC 1
#   5.8859563944     0.0000000000     0.0000000000
#   0.0000000000     5.8859563944     0.0000000000
#   0.0000000000     0.0000000000    23.3615963678
# PRIMCOORD 1
#32 1
# 47    3.8888048966     4.1163007841    19.5974532965
# 47    3.8888048966     1.7696556103    19.5974532965
# 47    1.9971514978     4.1163007841    19.5974532965
# 47    1.9971514978     1.7696556103    19.5974532965
#  6    2.9429781972     3.5961067857    19.2368415845
#  6    2.9429781972     2.2898496088    19.2368415845
#  ..............

######################################################################
# xsfParseAnim --
#     Converts AXSF file to sequence of XSF files. If the AXSF tail of 
#     filename is "name", then the generated sequence of XSF files will
#     be named "name.n", where "n" is the sequent number starting from 
#     1 and up. The files are generated in `pwd`.
#     
# Arguments:
#     file      file is the the FULL-PATH/FILENAME
#
# Results:
#     Returns the number of XSF files genetared. It return 0 if
#     ANIMSTEP keyword is not found and -1 if parsing error occurs
proc xsfParseAnim file {   

    # default type of animated XSF file
    set axsf_type fixed-cell

    set filetail [file tail $file]

    set file [gunzipXSF $file]

    set content [ReadFile -nonewline $file]
    set ConList [split $content \n]
    unset content

    # read NUMBER of ANIMATION STEPS    
    set AnimSteps [lindex [lindex $ConList \
	    [lsearch -glob $ConList "*ANIMSTEP*"]] 1]
    if { $AnimSteps < 0 } {
	return 0
    }

    # read HEADER of FILE
    set ind  [lsearch -glob $ConList "*DIM-GROUP*"]
    if { $ind >= 0 } {	
	lappend HList [lrange $ConList $ind [expr $ind + 1]]
    }

    set ind  [lsearch -regexp $ConList "^ *MOLECULE*"]
    if { $ind >= 0 } {	
	lappend HList [linedex $ConList $ind]
    }
    set ind  [lsearch -regexp $ConList "^ *POLYMER*"]
    if { $ind >= 0 } {	
	lappend HList [lindex $ConList $ind]
    }
    set ind  [lsearch -regexp $ConList "^ *SLAB*"]
    if { $ind >= 0 } {	
	lappend HList [lindex $ConList $ind]
    }
    set ind  [lsearch -regexp $ConList "^ *CRYSTAL*"]
    if { $ind >= 0 } {	
	lappend HList [lindex $ConList $ind]
    }

    set ind [lsearch -glob $ConList "*PRIMVEC*"]
    if { $ind >= 0 } {
	# check if it is variable- or fixed-cell animated XSF
	if { [llength [lindex $ConList $ind]] == 1 } {	    
	    # fixed-cell XSF
	    set axsf_type "fixed-cell"
	    lappend HList [lrange $ConList $ind [expr $ind + 3]]
	} else {
	    # variable-cell XSF
	    set axsf_type "variable-cell"	    
	}
    }
    set ind [lsearch -glob $ConList "*CONVVEC*"]
    if { $ind >= 0 } {
	if { [llength [lindex $ConList $ind]] == 1 } {
	    # fixed-cell XSF
	    lappend HList [lrange $ConList $ind [expr $ind + 3]]
	}
    }
    # transform HList for writing to XSF file
    if { [info exists HList] } {
	foreach record $HList {
	    foreach elem $record {
		append Header ${elem}\n
	    }
	}
	unset HList
    } else {
	set Header {}
    }

    # check if there is a PRIMCOORD record and if it is composed like:
    # PRIMCOORD <n>, where n is a number
    set ind [lsearch -glob $ConList "*PRIMCOORD*"]
    if { $ind >= 0 && [llength [lindex $ConList $ind]] > 1 } {
	# parse PRIMCOORD records and write sequence of XSF files
	for {set i 1} {$i <= $AnimSteps} {incr i} {
	    set FileName $filetail.$i

	    append FileContent $Header

	    #
	    # for variable-cell XSF we should write the vectors !!!
	    #
	    if { $axsf_type == "variable-cell" } {
		# PRIMVEC
		set ind [lsearch -glob $ConList "*PRIMVEC*$i*"]
		if { $ind < 0 } {
		    # mandatory keyword: ERROR parsing AXSF file
		    return -1
		}
		append FileContent [_xsfAnimParse_getVec $ConList $ind]
		
		# CONVVEC
		set ind [lsearch -glob $ConList "*CONVVEC*$i*"]
		if { $ind >= 0 } {
		    # optional keyword
		    append FileContent [_xsfAnimParse_getVec $ConList $ind]
		}
	    }

	    # PRIMCOORD
	    set ind [lsearch -glob $ConList "*PRIMCOORD*$i*"]	    

	    if { $ind < 0 } {
		# mandatory keyword: ERROR parsing AXSF file
		return -1
	    }
	    incr ind
	    set Natm  [lindex [lindex $ConList $ind] 0]
	    set AList [lrange $ConList $ind [expr $ind + $Natm]]

	    append FileContent "PRIMCOORD\n"
	    foreach elem $AList {
		append atoms ${elem}\n
	    }
	    append FileContent $atoms
	    WriteFile $FileName $FileContent w	    
	    unset FileContent atoms
	}
    } else {
	# must be ATOMS <n> type of AXSF file
	set ind [lsearch -glob $ConList "*ATOMS*"]
	if { $ind >= 0 && [llength [lindex $ConList $ind]] > 1 } {
	    # parse ATOMS records and write sequence of XSF files
	    for {set i 1} {$i <= $AnimSteps} {incr i} {
		set FileName $filetail.$i
		set ind [lsearch -glob $ConList "*ATOMS*$i*"]
		if { $ind < 0 } {
		    # ERROR parsing AXSF file
		    return -1
		}
		# ATOMS section does not specify how many atoms there is,
		# hence search for next ATOMS <n> string !!!
		if { $i < $AnimSteps } {
		    set ii [expr $i + 1]
		    set ind1  [lsearch -glob $ConList "*ATOMS*$ii*"]	    
		    set AList [lrange $ConList \
				   [expr $ind + 1] [expr $ind1 - 1]]
		} else {
		    set AList [lrange $ConList [expr $ind + 1] end]
		}		    
		
		#append FileContent $Header
		append FileContent ATOMS\n
		foreach elem $AList {
		    append atoms ${elem}\n
		}
		if { [info exists atoms] } {
		    append FileContent $atoms
		    WriteFile $FileName $FileContent w
		    #xcDebug -stderr "ANIMSTEP: $i"
		    #xcDebug -stderr $FileContent
		    unset FileContent atoms
		}
	    }
	}
    }

    return $AnimSteps
}


#
# utility library: returns the whole PRIMVEC or CONVVEC section
#
proc _xsfAnimParse_getVec {ConList ind} {
    append _vec [lindex [lindex $ConList $ind] 0]\n
    incr ind
    append _vec [lindex $ConList $ind]\n
    incr ind
    append _vec [lindex $ConList $ind]\n
    incr ind
    append _vec [lindex $ConList $ind]\n
    return $_vec
}
