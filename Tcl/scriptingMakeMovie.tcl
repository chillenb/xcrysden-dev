#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/scriptingMakeMovie.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################


# ------------------------------------------------------------------------
#****c* Scripting/scripting::makeMovie
#
# NAME
# scripting::makeMovie
#
# PURPOSE

# This namespace provide the scripting interface for making
# MPEG/Animated-GIF movies. The scripting::makeMovie namespace
# interface is kind of state machine. The "init" command initializes
# the process, and then the movie creation is encapsulated between
# "begin" and "end" commands. Every movie frame is created within the
# "begin" and "end" commands with "makeFrame" command. One can create
# several movies within one scripting file. Make sure the sequence of
# makeMovie calls will have the following order:
#
# init
# begin
#   makeFrame
#   makeFrame
#   ...
# end
#
# init
# begin
#   makeFrame
#   makeFrame
#   ...
# end
#
#
# COMMANDS

# -- scripting::makeMovie::init

# Initializes the movie creation process. It must be called before
# making movie. One passes several configuration options to init call.


#
# -- scripting::makeMovie::begin

# Marks the begining of movie creation.

#
# -- scripting::makeMovie::makeFrame

# Makes one movie frame, that is, saves (i.e. prints to file) the
# content of the currently displayed object.

#
# -- scripting::makeMovie::end

# Finishes the movie creation and encodes the movie.
#
#****
# ------------------------------------------------------------------------

namespace eval scripting::makeMovie {
    variable movie

    set movie(makeMovie) 0
    set movie(movieFile) ""
    set movie(notificationMovie) ""
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::makeMovie::init
#
# NAME
# scripting::makeMovie::init
#
# USAGE
# scripting::makeMovie::init -option value ?-option value? ...
#
# PURPOSE

# This proc initializes the movie creation processes. One can pass several
# configuration options, which determine the technical details of the
# movie.

#
# ARGUMENTS
# args -- various configuration "-option value" pairs
#
# OPTIONS
# ------------------------------------------------------------------------
#  OPTION               ALLOWED-VALUES + DESCRIPTION
# ------------------------------------------------------------------------
#
#  -gif_transp          0|1 --> make oblique|transparent animated-GIF
#
#  -gif_minimize        0|1 --> don't-minimize|minimize animateg-GIF
#
#  -gif_global_colormap 0|1 --> don't-use|use global colormap for animated-GIF
#
#  -movieformat         avi|mpeg|gif --> create AVI|MPEG|Animated-GIF
#
#  -dir                 tmp|pwd --> put temporary (i.e. frame) files to 
#                                   scratch(tmp) or current working
#                                   directory(pwd) 
#
#  -frameformat         PPM|PNG|JPEG --> format of the frame-files 
#
#  -firstframe          positive-integer --> repeat first frame n-times
#
#  -lastframe           positive-integer --> repeat first frame n-times
#
#  -delay               positive-integer --> time dalay between frames 
#                                            in 1/100 sec
#  -loop                repeat in movie animation number of times (0=forever)
#
#  -save_to_file        file --> if specified the movie will be saved to file
#                                otherwise the filename will be queried
# 
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::makeMovie::init \
#    -movieformat mpeg \
#    -dir         tmp \
#    -frameformat PPM \
#    -firstframe  10 \
#    -lastframe   10 \
#    -delay       0
#
#****
# ------------------------------------------------------------------------

proc scripting::makeMovie::init {args} {
    variable movie
    global   gifAnim
    
    if { $movie(makeMovie) } {
	error "makeMovie::init called within makeMovie::begin/makeMovie::end, should be called before makeMovie::begin"
    }

    # load defaults
    
    set gifAnim(gif_transp)          0
    set gifAnim(gif_minimize)        0
    set gifAnim(gif_global_colormap) 0
    set gifAnim(edit_param)          1
    set gifAnim(movie_format)        mpeg
    set gifAnim(temp_files_dir)      tmp
    set gifAnim(frame_files_format)  PPM
    set gifAnim(ntime_first_frame)   1
    set gifAnim(ntime_last_frame)    1
    set gifAnim(delay)               0
    set gifAnim(loop)                0
    
    set movie(movieFile)             ""
    
    set i 0
    foreach option $args {
	incr i
	
	if { $i%2 } {
            set tag $option
	} else {	    
            switch -glob -- $tag {
		"-gif_transp" -         
		"-gif_minimize" -       
		"-gif_global_colormap" -
		"-edit_param" {         
		    set tag [string trimleft $tag -]
		    switch $option {
			1 - on - yes { set gifAnim($tag) 1 }
			0 - off - no { set gifAnim($tag) 0 }
			default {
			    ErrorDialog "wrong value $option for $tag option, should be 0 or 1"
			}
		    }
		} 
		"-movieformat" {
		    set option [string tolower $option]
		    switch $option {
			mpeg - gif { set gifAnim(movie_format) $option }
			avi { set gifAnim(movie_format) mpeg }
			default {
			    ErrorDialog "wrong value $option for $tag option, should be \"mpeg\" or \"gif\""
			}
		    }
		} 

		"-dir" {
		    switch $option {
			tmp - pwd { set gifAnim(temp_files_dir) $option }
			default {
			    ErrorDialog "wrong value $option for $tag option, should be \"tmp\" or \"pwd\""
			}
		    }
		} 

		"-frameformat" { 
		    set option [string toupper $option]
		    switch $option {
			PPM - PNG - JPEG { set gifAnim(frame_files_format) $option }
			default {
			    ErrorDialog "wrong value $option for $tag option, should be \"PPM\" or \"JPEG\""
			}
		    }
		} 
		"-ntime_first_frame" - 
		"-firstframe" {
		    if { [nonnegativeInteger $option] } {
			set gifAnim(ntime_first_frame) $option
		    } else {
			ErrorDialog "expected integer, but got $option for $tag option"
		    }
		}
		"-ntime_last_frame" - 
		"-lastframe" {
		    if { [nonnegativeInteger $option] } {
			set gifAnim(ntime_last_frame) $option
		    } else {
			ErrorDialog "expected integer, but got $option for $tag option"
		    }
		}
		"-delay" -
		"-loop" {		    
		    if { [nonnegativeInteger $option] } {
			set elem [string trimleft $tag -]
			set gifAnim($elem) $option
		    } else {
			ErrorDialog "expected integer, but got $option for $tag option"
		    }
		}
		"-save_to_file" {
		    set movie(movieFile) $option
		}		
	    }
	}
    }
    if { $i%2 } {
	ErrorDialog "scripting::makeMovie::init called with an odd number of arguments !!!"
    }
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::makeMovie::begin
#
# NAME
# scripting::makeMovie::begin
#
# USAGE
# scripting::makeMovie::begin
#
# PURPOSE
#
# This proc marks the beginning of movie creation.
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::makeMovie::begin
#
#****
# ------------------------------------------------------------------------


proc scripting::makeMovie::begin {} {
    variable movie
    global gifAnim

    set movie(makeMovie) 1
    
    if { $gifAnim(movie_format) == "mpeg" } {
	set fmt MPEG
    } else {
	set fmt Animated-GIF
    }
    set movie(notificationMovie) [DisplayUpdateWidget "Recording" "Recording $fmt movie."]
    set gifAnim(make_gifAnim) 0
    gifAnimMake .mesa
}



# ------------------------------------------------------------------------
#****f* Scripting/scripting::makeMovie::makeFrame
#
# NAME
# scripting::makeMovie::makeFrame
#
# USAGE
# scripting::makeMovie::makeFrame
#
# PURPOSE

# This proc makes one movie frame, that is, it flushes (i.e. prints
# to file) the content of the currently displayed object.

#
# WARNINGS

# Note that this proc should be called within scripting::makeMovie::begin 
# and scripting::makeMovie::end calls. The scripting::makeMovie::init
# should be called before "begin; makeFrame; ...; end" sequence.

#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::makeMovie::makeFrame
#
#****
# ------------------------------------------------------------------------

proc scripting::makeMovie::makeFrame {} {
    variable movie

    if { $movie(makeMovie) } {
	gifAnimPrintCurrent .mesa
    } else {
	error "makeMovie::makeFrame called outside makeMovie::begin/makeMovie::end"
    }
}


# ------------------------------------------------------------------------
#****f* Scripting/scripting::makeMovie::end
#
# NAME
# scripting::makeMovie::end
#
# USAGE
# scripting::makeMovie::end
#
# PURPOSE
# This proc finishes the movie creation and encodes the movie.
#
# RETURN VALUE
# Undefined.
#
# EXAMPLE
# scripting::makeMovie::end
#
#****
# ------------------------------------------------------------------------

proc scripting::makeMovie::end {} {
    variable movie
    global gifAnim

    if { [winfo exists $movie(notificationMovie)] } {
	destroy $movie(notificationMovie)
    }
    if { $gifAnim(filelist) == "" } { return }

    if { $movie(makeMovie) } {
	if { $movie(movieFile) == "" } {
	    gifAnimMake .mesa
	} else {
	    gifAnimMake .mesa {} {} $movie(movieFile)
	}
	set movie(makeMovie) 0
    } else {
	error "makeMovie::end called before makeMovie::begin"
    }
}
