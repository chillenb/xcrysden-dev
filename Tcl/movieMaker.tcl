#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/movieMaker.tcl                                       
# ------                                                                    #
# Copyright (c) 2008 by Anton Kokalj                                        #
#############################################################################

proc MovieMaker {togl} {
    global movieMaker xcMisc

    if { ![info exists xcMisc(gif_encoder)] && ![info exists xcMisc(movie_encoder)] } {
	# at least one of encoder must exist to create movies
	return
    }

    if { ! [info exists movieMaker(mode)] } {
	set movieMaker(mode) realtime 
    }

    if { ! [info exists movieMaker(fps)] } {
	set movieMaker(fps) 20
    }

    if { ! [info exists movieMaker(count)] } {
	set movieMaker(count) 0
    }

    if { ! [info exists movieMaker(recording_text)] } {
	set movieMaker(recording_text) "Start recording"
    }
   
    set c [lindex [split $togl .] end]
    if { [info exists .movie$c] } {
	return
    }

    # toplevel

    set t [xcToplevel .movie$c "Movie Maker" "Movie" . -0 0 1]
    bind $t <Destroy> {movieMaker_close %W}

    set movieMaker(tplw) $t
    set movieMaker(togl) $togl
    set movieMaker(recording) 0

    # this is for "Hide" button

    global unmapWin
    xcRegisterUnmapWindow $t $unmapWin(frame,main) movie \
	    -textimage {MovieMaker unmap}
    bind $t <Unmap> [list xcUnmapWindow unmap %W $t $unmapWin(frame,main) movie]
    bind $t <Map>   [list xcUnmapWindow map %W $t $unmapWin(frame,main) movie]

    # container frames

    set f1 [labelframe $t.f1 -text "Movie creation mode:"]
    set f2 [labelframe $t.f2 -text "Movie encoding options:"]
    set f3 [frame $t.f3]
    set movieMaker(f1) $f1
    set movieMaker(f2) $f2
    set movieMaker(f3) $f3

    pack $f1 -side top -expand 1 -fill x -padx 2m -pady 2m -ipadx 2m -ipady 2m 
    pack $f2 -side top -expand 1 -fill x -padx 2m -ipadx 5m -ipady 2m 
    pack $f3 -side top
    
    pack [set movieMaker(mode_frame) [frame $f1.f]] -side top -expand 1 -fill x -padx 0 -pady 0 -ipadx 0 -ipady 0

    # mode radiobuttons

    foreach text {
	"real-time capture" "capture every redisplay" "manual"
    } value {
	realtime everysnapshot manual
    } {
	radiobutton $f1.f.$value -text $text -variable movieMaker(mode) -value $value -anchor w
	pack $f1.f.$value -side top -padx 5 -pady 0 -fill x -expand 1	
    }


    set movieMaker(fps_entryframe) [frame $f1.fe]
    set movieMaker(buttonframe)    [frame $f1.bf]
    pack $f1.fe -side top -expand 1 -fill x -padx 5
    pack $f1.bf -side top -expand 1 -fill x -padx 5

    # frequency entry for real-time-capture

    label $f1.fe.fpsl -text "Real-time capture frequency (frames per second):" -anchor w
    entry $f1.fe.fpse -relief sunken -textvariable movieMaker(fps) -width 5
    pack $f1.fe.fpsl -side left -padx 5 
    pack $f1.fe.fpse -side left -fill x -expand 1

    # buttons

    set movieMaker(startstop_button)    [button $f1.bf.startstop -textvariable movieMaker(recording_text) -command movieMaker_recording_]
    set movieMaker(startstop_button_bg) [$f1.bf.startstop cget -background]
    set movieMaker(snapshot_button)     [button $f1.bf.snap -image snapshot -command movieMaker_frame]
    pack $movieMaker(startstop_button)  $movieMaker(snapshot_button) -side left -pady 3 -padx 3

    # Animated GIF/MPEG/AVI control widgets (from gifAnim.tcl)

    gifAnim_controlWidgets $f2 $togl

    #
    movieMaker_mode_
    trace add variable movieMaker(mode) write movieMaker_mode_

    # Close/Hide widgets

    set hide  [button $f3.hide  -text "Hide" -default active -command [list HideWin $t movie]]
    set close [button $f3.close -text "Close" -command [list movieMaker_close $t]]
    pack $hide $close -side left -ipadx 1m -ipady 1m -padx 2m -pady 2m

    xcUpdateState
}

proc movieMaker_frame {} {
    global movieMaker

    if { ! [info exists movieMaker(count)] } {
	set movieMaker(count) 0
    }

    incr movieMaker(count)
    scripting::makeMovie::makeFrame
}

proc movieMaker_end {} {
    global movieMaker
    
    if { $movieMaker(mode) != "manual" } {
	global gifAnim       
	set gifAnim(filelist) [xc_realtimemovie $movieMaker(togl) filelist]
	xc_realtimemovie $movieMaker(togl) end
	if { $movieMaker(mode) == "realtime" } {
	    $movieMaker(togl) configure -time $movieMaker(togl_timedelay)
	}
    }

    scripting::makeMovie::end
    set movieMaker(count) 0
}


proc movieMaker_manual {} {    
    global movieMaker

    set movieMaker(mode) "manual"
    scripting::makeMovie::begin
}


proc movieMaker_realTimeCapture {{ms 50}} {
    global movieMaker gifAnim system

    scripting::makeMovie::begin

    set gifAnim(frame_files_format) PPM

    set movieMaker(mode) "realtime"
    set movieMaker(togl_timedelay) [lindex [$movieMaker(togl) configure -time] end]
    
    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
    } else {
	set dir $system(SCRDIR)
    }

    $movieMaker(togl) configure -time $ms
    puts stderr "xc_realtimemovie $movieMaker(togl) begin realtime $dir"
    xc_realtimemovie $movieMaker(togl) begin realtime $dir
}
	       

proc movieMaker_everySnapshot {} {
    global movieMaker gifAnim system

    scripting::makeMovie::begin

    set gifAnim(frame_files_format) PPM
    set movieMaker(mode) "everysnapshot"
    
    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
    } else {
	set dir $system(SCRDIR)
    }
    xc_realtimemovie $movieMaker(togl) begin everysnapshot $dir
}
	       

proc movieMaker_mode_ {args} {
    global movieMaker gifAnim

    switch -glob -- $movieMaker(mode) {
	everysnap* {
	    xcDisableAll $movieMaker(fps_entryframe) $movieMaker(snapshot_button) $gifAnim(radio_frame,3)
	    movieMaker_destroy_manualWidgets_
	    set gifAnim(frame_files_format) PPM
	}
	realtime* {
	    xcDisableAll $movieMaker(snapshot_button) $gifAnim(radio_frame,3)
	    xcEnableAll  $movieMaker(fps_entryframe)
	    movieMaker_destroy_manualWidgets_
	    set gifAnim(frame_files_format) PPM
	}
	manual {
	    xcDisableAll $movieMaker(fps_entryframe)
	    xcEnableAll  $gifAnim(radio_frame,3)
	    _mencoder_temporary_files
	    if { $movieMaker(recording) } {
		xcEnableAll $movieMaker(snapshot_button)
		if { ! [winfo exists $movieMaker(togl).snap] } {
		    button $movieMaker(togl).snap -image snapshot -command movieMaker_frame
		    pack $movieMaker(togl).snap -side left -expand 0 -fill none -padx 0 -pady 0 -ipadx 0 -ipady 0 -anchor nw
		}

		if { ! [winfo exists $movieMaker(buttonframe).frame] } {
		    pack [label $movieMaker(buttonframe).frame -textvariable movieMaker(count)] -side right -padx 2
		    pack [label $movieMaker(buttonframe).label -text "No. of recorded frames: "] -side right -padx 2
		}
	    } else {
		xcDisableAll $movieMaker(snapshot_button)
	    	movieMaker_destroy_manualWidgets_
	    }
	}
    }
}


proc movieMaker_destroy_manualWidgets_ {} {
    global movieMaker
    if { [winfo exists $movieMaker(togl).snap] } {
	destroy $movieMaker(togl).snap
    }
    if { [winfo exists $movieMaker(buttonframe).frame] } {
	destroy $movieMaker(buttonframe).label $movieMaker(buttonframe).frame
    }
}


proc movieMaker_recording_ {} {
    global movieMaker
    
    if { ! $movieMaker(recording) } {
	
	# start recording

	set movieMaker(recording) 1
	set movieMaker(recording_text) "Stop recording"
	$movieMaker(startstop_button) config -background \#ff4444
	xcDisableAll $movieMaker(mode_frame)

	movieMaker_mode_ 

	switch -glob -- $movieMaker(mode) {
	    everysnap* {
		movieMaker_everySnapshot
	    }
	    realtime* {
		set ms [expr round(1.0 / $movieMaker(fps) * 1000)]
		puts stderr "fps = $movieMaker(fps) ; ms = $ms"
		movieMaker_realTimeCapture $ms
	    }
	    manual {
		movieMaker_manual
	    }
	}
    } else {

	# stop recording

	movieMaker_end

	set movieMaker(recording) 0
	set movieMaker(recording_text) "Start recording"
	$movieMaker(startstop_button) config -background $movieMaker(startstop_button_bg)
	xcEnableAll $movieMaker(mode_frame)
	movieMaker_mode_ 
    }
}
 

proc movieMaker_clear {} {
    global movieMaker movie
    
    set movieMaker(count)     0
    set movieMaker(recording) 0
    
    if { [info exists ::scripting::makeMovie::movie(notificationMovie)] } {
	if { [winfo exists $::scripting::makeMovie::movie(notificationMovie)] } {
	    destroy $::scripting::makeMovie::movie(notificationMovie)
	}    
    }

    xc_realtimemovie $movieMaker(togl) clear
    movieMaker_destroy_manualWidgets_
}


proc movieMaker_close {t} {
    movieMaker_clear

    set movieMaker(tplw) {}
    destroy $t

    xcUpdateState
}
    
