#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/gifAnim.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc gifAnimWidPack {frame can} {
    global gifAnim

    if { $gifAnim(button_text) == "Animated GIF/MPEG/AVI >>" } {
	pack $frame -side top -padx 0 -pady 2
	set gifAnim(button_text) "<< Animated GIF/MPEG/AVI"
	set gifAnim(anim_ctrl_widgets) 1

	if { [info exists gifAnim(radioframes,$frame)] } {
	    # dirty hack
	    set i 0
	    foreach rf $gifAnim(radioframes,$frame) {
		set gifAnim(radio_frame,$i) $rf
		incr i
	    }
	}

	_mencoder_temporary_files
    } else {
	pack forget $frame
	set gifAnim(button_text) "Animated GIF/MPEG/AVI >>"
	set gifAnim(anim_ctrl_widgets) 0
    }

    xcUpdateState 
}


proc gifAnimMake {can {play_frame ""} {record_button ""} {outfile {}}} {
    global gifAnim xcMisc mesa_bg system
    
    if { $gifAnim(make_gifAnim) == 0 } {
	
	# --------------------------------------------------
	# The "Start Recording Animation" button was pressed
	#
	# This means we start recording. Also change the label of the 
	# button to "Stop Recording and Save"
	# --------------------------------------------------

	set gifAnim(make_gifAnim) 1
	set gifAnim(filelist)     {}  
	set gifAnim(nframe)       0
	
	if { $play_frame != "" } {
	    set gifAnim(make_text)        "Stop Recording and Save"
	    set gifAnim(play_frame_bg)    [$play_frame cget -bg]
	    set gifAnim(record_button_bg) [$record_button cget -bg]	
	    $play_frame    configure -background \#ff4444
	    $record_button configure -background \#ff4444
	    $record_button configure -relief sunken
	}
	if { [info exists gifAnim(radio_frame,1)] } {
	    if { [winfo exists $gifAnim(radio_frame,1)] } {
		# disable the radiobuttons	
		xcDisableAll  $gifAnim(radio_frame,0) $gifAnim(radio_frame,1) $gifAnim(radio_frame,2) $gifAnim(radio_frame,3) $gifAnim(radio_frame,4)
	    }
	}
    } else {


	# --------------------------------------------------
	# The "Stop Recording and Save" button was pressed
	#
	# This means we save the Animated-GIF/MPEG file. Also change the 
	# label of the  button to "Start Recording Animated-GIF/MPEG"
	# --------------------------------------------------

	set gifAnim(make_gifAnim) 0
	
	if { $play_frame != "" } {
	    set gifAnim(make_text)   "Start Recording Animation"
	    $play_frame    configure -background $gifAnim(play_frame_bg)
	    $record_button configure -background $gifAnim(record_button_bg)
	    $record_button configure -relief raised
	    # enable the radiobuttons	
	}
	if { [info exists gifAnim(radio_frame,1)] } {
	    if { [winfo exists $gifAnim(radio_frame,1)] } {
		xcEnableAll $gifAnim(radio_frame,0) $gifAnim(radio_frame,1) $gifAnim(radio_frame,2) $gifAnim(radio_frame,3) $gifAnim(radio_frame,4)
		_mencoder_temporary_files gifAnim movie_format ""	    
	    }
	}

	set first_frame [lindex $gifAnim(filelist) 0]
	set last_frame  [lindex $gifAnim(filelist) end]
	set first {}
	set last  {}
	if { $gifAnim(ntime_first_frame) > 1 } {
	    for {set i 1} {$i < $gifAnim(ntime_first_frame)} {incr i} {
		append first "${first_frame} "
	    }
	}
	if { $gifAnim(ntime_last_frame) > 1 } {
	    for {set i 1} {$i < $gifAnim(ntime_last_frame)} {incr i} {
		append last "${last_frame} "
	    }
	}
	set filelist [format "%s %s %s" $first $gifAnim(filelist) $last]

	if { $gifAnim(movie_format) == "gif" } {

	    # --------------------------------------------------
	    # create Animated GIF
	    # --------------------------------------------------
	    if { $outfile == "" } {
		set outfile [tk_getSaveFile -defaultextension .gif \
				 -filetypes { 
				     {{GIF File}  {.gif}}
				     {{All Files} {.*}}
				 } \
				 -defaultextension .gif \
				 -initialdir  $system(PWD) \
				 -title       "Save Animated GIF File"]	    
		if { $outfile == "" } {
		    return
		}
	    }
	    
	    encode_gif $filelist $outfile
	} else {

	    # --------------------------------------------------
	    # create AVI/MPEG
	    # --------------------------------------------------

	    if { $outfile == "" } {
		set filetypes { 
		    {{AVI File}   {.avi}}
		    {{MPEG File}  {.mpg .mpeg}}
		    {{All Files}  {.*}}
		}
		set defExt .avi
		if { $xcMisc(movie_encoder) == "ppmtompeg" } {
		    set filetypes { 
			{{MPEG File}  {.mpg .mpeg}}
			{{All Files}  {.*}}
		    }
		    set defExt .mpg
		}
		
		set outfile [tk_getSaveFile -defaultextension .mpg \
				 -filetypes $filetypes \
				 -defaultextension $defExt \
				 -initialdir  $system(PWD) \
				 -title       "Save AVI/MPEG File"]
		if { $outfile == "" } {
		    return
		}
	    }
	    
	    encode_movie $filelist $outfile
	}
    }
}

	
proc gifAnimPrintCurrent {can} {
    global gifAnim xcMisc mesa_bg system movieMaker

    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
    } else {
	set dir $system(SCRDIR)
    }
    
    if { [info exists movieMaker(mode)] && [info exists movieMaker(recording)] } {
	# movieMaker's everysnapshot and realtime does not use printing via this proc

	if { $movieMaker(mode) != "manual" && $movieMaker(recording) } {
	    return
	}
    }

    set abs_head [file join $dir anim-[format "%04d" $gifAnim(nframe)]]
    set ext      [string tolower $gifAnim(frame_files_format)]
	
    if { $gifAnim(movie_format) == "gif" } {
	set abs_file $abs_head.gif
	if { [info exists xcMisc(gif_encoder)] } {
	    if { $xcMisc(gif_encoder) == "convert" } {
		# ImageMagick's convert can use whatever format
		set abs_file $abs_head.$ext
	    }    
	}
	append gifAnim(filelist) " $abs_file"
    } else {
	set abs_file $abs_head.$ext
	append gifAnim(filelist) " $abs_file"
    }

    incr gifAnim(nframe)

    # now update the display (maybe $can was obscured)
    scripting::printToFile $abs_file
}


proc gifAnim:editParam {comlin} {
    global gifAnim system
    
    set comlin [subst {
# ------------------------------------------------------------------------
# You may edit the movie encoder command below (don't touch this comment)
# ------------------------------------------------------------------------

$comlin
}]

    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
    } else {
	set dir $system(SCRDIR)
    }	    

    set file [file join $dir encode.sh]
    WriteFile $file $comlin w
    xcEditFile $file foreground

    # drop the comments from the file

    set content {}
    foreach line [split [ReadFile -nonewline $file] \n] {
	if { ! [regexp {^#} [string trim $line]] } {
	    if { $line != {} } {
		append content [format "%s\n" $line]
	    }
	}
    }
    return $content
}

  
proc gifAnim_controlWidgets {parent can {slide_label ""}} {    
    global gifAnim
    
    # groovy frames for checkbuttons+radiobuttons+entries
    
    set gifAnim(radioframes,$parent) ""

    foreach i {
	0 1 2 3 4 5
    } text {
	{Animated-GIF options:}
	{Movie format:}
	{Store temporary image frame files in:}
	{Format of temporary frame files:}
	{Movie options:}
	{Encoding:}
    } {	    
	set f($i) [labelframe $parent.$i -text $text -bd 2]
	pack $f($i) -side top -padx 3 -pady 5 -expand 1 -fill x

	set gifAnim(radio_frame,$i) $f($i)
	lappend gifAnim(radioframes,$parent) $f($i)
    }    

    # checkbuttons

    foreach \
	cb    [list $f(0).cb1 $f(0).cb2 $f(0).cb3] \
	elem  [list gif_global_colormap gif_minimize gif_transp] \
	onvalue "1 1 1" offvalue "0 0 0" \
	text {
	    "use global color-map"
	    "minimize GIF file size"
	    "make transparent background"
	} {
	    checkbutton $cb -text $text -variable gifAnim($elem) \
		-onvalue $onvalue -offvalue $offvalue -anchor w
	    pack $cb -side top -padx 5 -pady 0 -expand 1 -fill x
	}
    
    # radiobuttons

    foreach \
	rb    [list $f(1).r1 $f(1).r2   $f(2).r1 $f(2).r2   $f(3).r1 $f(3).r2 $f(3).r3] \
	pack  {left left   left left   top top top} \
	value {mpeg gif    pwd tmp   PPM PNG JPEG} \
	elem  {
	    movie_format movie_format
	    temp_files_dir temp_files_dir 
	    frame_files_format frame_files_format frame_files_format} \
	text  {
	    "AVI/MPEG"
	    "Animated-GIF" 
	    
	    "current working directory"
	    "scratch directory"
	    
	    "non-compressed PPM"
	    "compressed PNG"
	    "compressed JPEG"
	} {
	    radiobutton $rb -text $text -variable gifAnim($elem) -value $value -anchor w
	    pack $rb -side $pack -padx 5 -pady 0 -fill x -expand 1
	}
    set gifAnim(radiobutton_movie_format_mpeg) $f(1).r1
    set gifAnim(radiobutton_movie_format_gif)  $f(1).r2

    _mencoder_temporary_files gifAnim movie_format ""
    trace add variable gifAnim(movie_format) write _mencoder_temporary_files

    # entries

    FillEntries $f(4) {
	{Repeat first frame No. times:}
	{Repeat last  frame No. times:}
	{Time delay between slides (1/100 sec):}
	{Loop animation number of times (0=forever):}
    } {
	gifAnim(ntime_first_frame)
	gifAnim(ntime_last_frame)
	gifAnim(delay)
	gifAnim(loop)
    } 38 8

    checkbutton $f(5).cb -text "Edit flags or parameter-file before encoding" \
	-variable gifAnim(edit_param) -onvalue 1 -offvalue 0 -anchor w
    pack $f(5).cb -side top -expand 1 -fill both -padx 5 -pady 2
}


proc _mencoder_temporary_files {args} {
    global gifAnim xcMisc
        
    if { ![info exists xcMisc(movie_encoder)] && ![info exists xcMisc(gif_encoder)]} {
	return
    }
    
    if { ! [info exists xcMisc(movie_encoder)] } {
	set gifAnim(movie_format) "gif"
	catch { $gifAnim(radiobutton_movie_format_mpeg) config -state disabled }    
    } elseif { ! [info exists xcMisc(gif_encoder)] } {
	set gifAnim(movie_format) "mpeg"
	catch { $gifAnim(radiobutton_movie_format_gif) config -state disabled }	
    }

    if { $gifAnim(movie_format) == "mpeg" && $xcMisc(movie_encoder) == "mencoder" } {
	puts stderr ::::::::MPEG
	xcDisableAll $gifAnim(radio_frame,0)
	xcDisableOne $gifAnim(radio_frame,3).r1

	set to_png 1
	global movieMaker
	if { [info exists movieMaker(recording)] } {
	    if { $movieMaker(recording) && $movieMaker(mode) != "manual" } {
		# for realtime and everysnapshot movie modes we don't convert to png on the fly
		set to_png 0
	    }
	}
	if { $gifAnim(frame_files_format) == "PPM" && $to_png } {
	    set gifAnim(frame_files_format) PNG
	}
    } else {
	puts stderr ::::::::GIF
	xcEnableAll $gifAnim(radio_frame,0)
	xcEnableOne $gifAnim(radio_frame,3).r1
    }
}
