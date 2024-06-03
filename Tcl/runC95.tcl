#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/runC95.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc RunAndGetC95Output {command bg_how input_var} {
    global system

    set out $system(SCRDIR)/xc_output.$system(PID)
    if { ! [RunC95 $command $bg_how $input_var $out ] } {
	return 0
    }

    return [ReadFile $out]
}


proc RunC95 {command bg_how input_var \
	{output_file {}} {is_string {}} {dir {}}} {
    global system prop runC95
    # command     ... one of integrals, scf, scfdir, properties
    #
    # bg_how      ... how to handle background execution:
    #            $bg_how == {} | {message <message_text>} | bg
    #            {}    ... calculation is fast; give no feed-back to the user
    #            message . calculation can be long; display some message to
    #                      the user
    #            bg    ... calculation is probably very long; 
    #                      just submit to the background and forget about it.
    #            "{}" & "message" works via fileevent, whereas bg submit the 
    #            job to background and forget about it.
    #
    # input_var   ... variable where input is stored
    #
    # output_file ... name of output file
    #
    # is_string   ... check if $is_string appear in crystal95's output

    # this is just in any case
    set runC95(output) {}
    #
    # check if the $command is allowed command
    #
    set match 0
    foreach com [list $system(c95_integrals) $system(c95_scf) \
	    $system(c95_scfdir) $system(c95_properties)] {
	if { $com == $command } { set match 1 }
    }
    if { ! $match } {
	tk_dialog [WidgetName] "ERROR" "ERROR: CRYSTAL-$system(c95_version)'s command \"${command}:\" not found, must be one of: $system(c95_integrals), $system(c95_scf), $system(c95_scfdir), $system(c95_properties)" error 0 OK
	return 0
    }

    #
    # if dir is not defined, use default one
    #
    if { $dir == "" } {
	set dir $system(SCRDIR)
    } 

    #
    # what is the output file ???
    #
    if { $output_file == "" } {
	set output_file $system(SCRDIR)/xc_output.$system(PID)
    }
    set runC95(output_file) $output_file

    set inp $dir/xc_tmp.$system(PID)
    WriteFile $inp $input_var

    #
    # if the species is EXTERNAL
    #
    cxxManageExternal

    ####CD DIR#######
    # now go to dir #
    cd $dir
    xcDebug -stderr "Crystal Input:\n--------------\n$input_var"
    xcDebug "PWD = [pwd]"

    if { $bg_how != "bg" } {
	if { [lindex $bg_how 0] == "message" } {
	    ###################################################################
	    # calculation can take some time, it's needed to 
	    # give some feed back to the user.
	    set oldgrab [grab current]
	    set mw [DisplayUpdateWidget "Calculating" [lindex $bg_how 1]]
	    catch { grab $mw }
	}
	
	#
	# change cursor to watch, to indicate that something is going on
	#
	SetWatchCursor
	
	# ================================
	# run CRYSTAL in background mode
	# ================================
	set runC95(command) $command
	set runC95(fileID) [open "| $command < $inp 2> /dev/null" r]
	fconfigure $runC95(fileID) -blocking 0
	fileevent $runC95(fileID) readable C95FileEvent
	
	tkwait variable runC95(done)
	ResetCursor
	xcSwapBuffers

	if { [lindex $bg_how 0] == "message" } {
	    destroy $mw
	    if { $oldgrab != {} } {
		catch { grab $oldgrab }
	    }
	}
	#
	# has an error occured during CRYSTAL execution
	#
	if { $runC95(error) == 1 } {
	    #if [catch {exec $command < $inp > $output_file} err]
	    #C95Error $command $output_file $err
	    cd $system(SCRDIR)
	    unset runC95
	    ResetCursor
	    return 0
	} else {
	    # if CRYSTAL has exited nicely, it does not mean that an
	    # error didn't occure; check for ERROR string (here
	    # CRYSTAL06 and later seems aexception when keyword STOP
	    # is used)
	    set is_str 0
	    set is_err 0
	    set fileID [open $output_file r]
	    foreach line [split [read $fileID] \n] {
		if { [string match "*ERROR \*\*\*\**" $line] } {		
		    # check for CRYSTAL06 line: "... STOP KEYWORD - EXECUTION STOPS"
		    if { ($system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14") \
			     && [string match "*STOP KEYWORD - EXECUTION STOPS*" $line] } {
			# OK, we don't have error
			continue
		    } else {
			close $fileID
			C95Error
			cd $system(SCRDIR)
			unset runC95
			ResetCursor
			return 0
		    }
		}
		if { $is_string != "" } {
		    if [string match "*${is_string}*" $line] {
			set is_str 1
		    }
		}
	    }
	    if { $is_string != "" && $is_str == 0} {
		tk_dialog [WidgetName] "CRYSTAL ERROR" "ERROR: pattern \"$is_string\" not found in CRYSTAL-$system(c95_version)'s output; this may be bug in program - please report to author !!! \nCode: RunC95 -pattern not found-" error 0 OK
		close $fileID
		cd $system(SCRDIR)
		unset runC95
		ResetCursor
		return 0
	    }
	    close $fileID
	    
	    # t.k: zakaj skopira unit 34 ??????
	    if { $dir != $system(SCRDIR) } {
		# copy ftn34 to $system(SCRDIR)
		if [file exists $dir/$prop(file)34] {
		    file copy -force $dir/$prop(file)34 $system(SCRDIR)/$prop(file)34
		}
	    }
	    cd $system(SCRDIR)
	    unset runC95
	    ResetCursor
	    xcDebug "end of RunC95"
	    return 1
	}
    }
}
    

proc C95FileEvent {} {
    global runC95
    
    if { ! [eof $runC95(fileID)] } {
	append runC95(output) [gets $runC95(fileID)]\n
	return
    } else {
	xcDebug "In C95FileEvent"
	set runC95(error) 0
	if { [catch {close $runC95(fileID)}] } {
	    xcDebug "In C95FileEvent: error while closing the process"
	    set runC95(error) 1
	}
	
	#
	# write $output to runC95(output_file)
	#
	set fileID [open $runC95(output_file) w]
	puts $fileID $runC95(output)
	flush $fileID
	close $fileID
	
	#
	# if we have an error display tk_dialog
	#
	if { $runC95(error) } {
	    C95Error 
	}
	
	set runC95(done) 1
    }
}


proc C95Error {} {
    global c95error runC95 system
    
    set button [tk_dialog [WidgetName] ERROR \
		    "ERROR occured while executing CRYSTAL-$system(c95_version) command: $runC95(command)" \
		    error 0 OK Details]

    if { $button == 1} {
	set oldgrab [grab current]
	set t [xcToplevel [WidgetName] "CRYSTAL ERROR" "ERROR" \
		. 100 100 1]
	catch { grab $t }
	
	set fileID [open $runC95(output_file) r]
	set output [read $fileID]
	close $fileID

	set text [DispText $t.f $output 80 35]
	set fric [lindex [$text yview] 1]
	$text yview moveto [expr 1.0 - $fric]
	$text config -state disabled
	
	proc C95ErrorCan {} {
	    global c95error
	    set c95error(done) 1
	}
	
	set f2 [frame $t.f2 -height 10]
	pack $f2 -side bottom -before $t.f -fill x
	set close [button $f2.cl -text "Close" \
		-command C95ErrorCan]
	pack $close -side left -expand 1 -ipadx 2 -ipady 2 -pady 10
	
	tkwait variable c95error(done)
	unset c95error
	
	catch { grab release $t }
	destroy $t	
	if { $oldgrab != "" } {
	    catch { grab $oldgrab }
	}
    }
}


#
# Purpose: driver proc for executing CRYSTAL for Advanced-Geom manipulations 
#
proc RunC95_advGeom {input option error_text} {
    global system

    # cd to SCRATCH dir
    cd $system(SCRDIR)
    set output xc_output.$system(PID)

    cxxManageExternal

    if { [catch {exec $system(c95_integrals) < $input > $output} errMsg] } {

	# printing the "FORTRAN STOP" to stderr cause the error. Check !!!
	
	if { $errMsg != "FORTRAN STOP" } {
	    # a real error occurred
	    ErrorDialog "an error occured for option: $option. Please try Again !!!" $errMsg
	    # delete last AdvGeom state
	    xcAdvGeomState delete
	    ResetCursor
	    return 0
	}
    }
    
    # if CRYSTAL has exited nicely, it does not mean that an error
    # didn't occure; check for ERROR (here CRYSTAL06 and later seems
    # an exception when keyword STOP is used)

    set file [ReadFile $output]
    foreach line [split $file \n] {
	if { [string match "*ERROR \*\*\*\**" $line] } {		
	    # check for CRYSTAL06 line: "... STOP KEYWORD - EXECUTION STOPS"
	    if { ($system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14") \
		     && [string match "*STOP KEYWORD - EXECUTION STOPS*" $line] } {
		# OK, we don't have error
		continue
	    } else {
		ErrorDialogInfo "$error_text Please try Again !!!" $file
		# delete last AdvGeom state
		xcAdvGeomState delete
		ResetCursor
		return 0
	    }
	}
    }    

    return 1
}
