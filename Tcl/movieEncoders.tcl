#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/__file__
# ------                                                                    #
# Copyright (c) 2008 by Anton Kokalj                                        #
#############################################################################

proc determine_movie_encoders {} {
    global xcMisc
    
    set gif_encoder_priority_list   {convert gifsicle whirlgif}
    set movie_encoder_priority_list {mencoder ppmtompeg}
    
    if { [info exists xcMisc(gif_encoder)] } {
	set gif_encoder_priority_list [linsert $gif_encoder_priority_list 0 [file tail $xcMisc(gif_encoder)]]
    }
    if { [info exists xcMisc(movie_encoder)] } {
	set movie_encoder_priority_list [linsert $movie_encoder_priority_list 0 [file tail $xcMisc(movie_encoder)]]
    }
    
    foreach enc $gif_encoder_priority_list {
	if { [info exists xcMisc($enc)] } {
	    set xcMisc(gif_encoder) $enc
	    break
	}
    }	 
    
    foreach enc $movie_encoder_priority_list {
	if { [info exists xcMisc($enc)] } {
	    set xcMisc(movie_encoder) $enc
	    break
	}
    }
}


# ------------------------------------------------------------------------
#
# AVI/MPEG movies
#
# ------------------------------------------------------------------------

proc encode_movie {filelist outfile} {
    global gifAnim xcMisc system

    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
	cd $system(PWD)
    } else {
	set dir $system(SCRDIR)
	cd $system(SCRDIR)
    }	    

    set cmdLine [cmdLine_$xcMisc(movie_encoder) $filelist $outfile $dir]

    set cw [DisplayUpdateWidget "Encoding" "Encoding the movie"]
    eval xcCatchExecReturnRedirectStdErr $cmdLine
    destroy $cw
}


proc cmdLine_ppmtompeg {filelist outfile input_dir} {
    global xcMisc
    return "$xcMisc(ppmtompeg) [mpegCreateParamFile $outfile $input_dir [join $filelist\n]]"
}


proc cmdLine_mencoder {filelist outfile input_dir} {
    global gifAnim xcMisc

    if { [regexp {\.ppm[ \n\r,]|\.ppm$} $filelist] } {
	# we need to convert all *.PPM files to *.PNG files
	set cw [DisplayUpdateWidget "Converting" "Converting frame files to PNG format"]
	foreach image $gifAnim(filelist) {
	    set newfile [file rootname $image].png
	    scripting::_printToFile_imageConvert $image $newfile
	}
	set tmplist $filelist
	regsub -all -- {\.ppm} $tmplist .png filelist
	destroy $cw
    }

    set delay [expr {$gifAnim(delay) <= 0.0 ? 1 : $gifAnim(delay)}]
    set fps [expr 100.0 / $delay]
    
    set cmdLine "$xcMisc(mencoder) mf://[join $filelist ,]"

    if { [string toupper [file extension $outfile]] == ".AVI" } {
	# create AVI file
	append cmdLine " -mf fps=$fps -o $outfile -ovc lavc -lavcopts vcodec=mpeg4"
    } else {
	# create MPEG file
	append cmdLine " -of mpeg -mf fps=$fps -o $outfile -ovc lavc -lavcopts vcodec=mpeg2video"
    }
    
    if { $gifAnim(edit_param) } {
	set cmdLine [gifAnim:editParam $cmdLine]
    }
    return $cmdLine    
}


# ------------------------------------------------------------------------
#
# Animated GIF
#
# ------------------------------------------------------------------------

proc encode_gif {filelist outfile} {
    global gifAnim xcMisc mesa_bg system

    # determine the background color

    if { $gifAnim(gif_transp) } {
	if { ! [info exists mesa_bg(current)] } {
	    set mesa_bg(current) \#000000
	}
	if { [string index $mesa_bg(current) 0] == "#" } {
	    set bg_color $mesa_bg(current)
	} else {
	    set bg_color [rgb_f2h $mesa_bg(current)]
	}
	set gifAnim(bg_color) $bg_color
    }

    # construct the execution command-line for encoding

    set cmdLine [cmdLine_$xcMisc(gif_encoder) $filelist $outfile]

    # encode ...

    if { $gifAnim(temp_files_dir) == "pwd" } {
	cd $system(PWD)
    } else {
	cd $system(SCRDIR)
    }
    if { $gifAnim(edit_param) } {
	set cmdLine [gifAnim:editParam $cmdLine]
    }
    
    puts stderr "cmdLine=$cmdLine"
    flush stderr
    
    set cw [DisplayUpdateWidget "Encoding" "Encoding the Animated-GIF movie"]
    if { $xcMisc(gif_encoder) != "gifsicle" } {
	eval xcCatchExecReturnRedirectStdErr $cmdLine
    } else {
	eval exec_gifsicle $cmdLine
    }
    destroy $cw
}


proc cmdLine_convert {filelist outfile} {
    global gifAnim xcMisc 
    
    if { $gifAnim(gif_transp) } {
	append comlin " -dispose Previous"
    }

    append comlin " $filelist -set delay $gifAnim(delay) -loop $gifAnim(loop)"

    if { $gifAnim(gif_global_colormap) && !$gifAnim(gif_transp) } {
	append comlin " +map"
    }
    if { $gifAnim(gif_transp) } {
	append comlin " -transparent $gifAnim(bg_color)"
    }
    if { $gifAnim(gif_minimize) } {
	append comlin " -layers optimize"
    }
    return [format "%s %s %s" $xcMisc(convert) $comlin $outfile]
}


proc cmdLine_gifsicle {filelist outfile} {
    global gifAnim xcMisc 
    
    set flags "--no-warnings --delay $gifAnim(delay) --loopcount=$gifAnim(loop)"
    if { $gifAnim(gif_transp) } {
	append flags " --disposal background --transparent=$gifAnim(bg_color)"
    }
    if { $gifAnim(gif_minimize) } {
	append flags " -O2"
    }
    if { $gifAnim(gif_global_colormap) } {
	# ignored
    }
    return [format "%s %s %s > %s" $xcMisc(gifsicle) $flags $filelist $outfile]
}


proc cmdLine_whirlgif {filelist outfile} {
    global gifAnim xcMisc 

    set flags " -time $gifAnim(delay) -loop $gifAnim(loop)"
    if { $gifAnim(gif_minimize) } {
	append flags " -minimize"
    }
    if { $gifAnim(gif_global_colormap) } {
	append flags " -globalmap"
    }
    if { $gifAnim(gif_transp) } {
	append flags " -disp prev -trans $gifAnim(bg_color)"
    }
    
    return [format "%s %s %s %s %s" $xcMisc(whirlgif) $flags -o $outfile $filelist]
}


proc exec_gifsicle {args} {
    xcDebug -stderr "Executing: $args"
    if { [catch {eval exec $args 2> gifsicle.stderr} errMsg] } {
	# gifsicle returns wrong status of 1 when -transparent option
	# is used
	set msg [ReadFile gifsicle.stderr]
	if { [string match {*Usage: gifsicle*} $msg] } {
	    ErrorDialogInfo "while executing\nexec $args" $msg\n$errMsg
	    uplevel 1 {
		return 1
	    }
	}
    }
    return 0
}