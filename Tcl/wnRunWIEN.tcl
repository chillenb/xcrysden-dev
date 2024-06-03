#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnRunWIEN.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnRun_IsERROR {} {
    global system
    if { [file size [exec $system(BINDIR)/wn_errorfile]] } {
	return 1
    }
    return 0    
}

proc wnRunWIEN {command message outf} {
    global wn runWn system xcMisc

    wnRunWIEN_DeleteOutputFile $outf

    set runWn(command) $command
    set runWn(outf) $outf
    #
    # one user has reported that "pipe" execution doesn't work on SUN;
    # provide "direct" execution as a way out
    #
    if { [info exists xcMisc(WIEN_direct_exe)] } {
	if $xcMisc(WIEN_direct_exe) {
	    cd $wn(dir)
	    ###################################################################
	    # calculation can take some time, it's needed to 
	    # give some feed back to the user.
	    set oldgrab [grab current]
	    set mw [DisplayUpdateWidget "Calculating" $message]
	    update
	    catch { grab $mw }
	    #if { [catch "exec $command" err] == 1 } 
	    catch {exec $command} err
	    if { [wnRun_IsERROR] } {
		destroy $mw
		wnRunWIENError 
		return 0
	    }
	    destroy $mw
	    return 1
	}
    }

    ###################################################################
    # calculation can take some time, it's needed to 
    # give some feed back to the user.
    set oldgrab [grab current]
    set mw [DisplayUpdateWidget "Calculating" $message]
    update
    catch { grab $mw }
    #
    # change cursor to watch, to indicate that something is going on
    #
    SetWatchCursor
	
    set runWn(outf) $outf

    cd $wn(dir)
    set command [concat | $command]

    # ================================
    # run WIENXX in background mode
    # ================================

    set runWn(command) $command
    set runWn(fileID) [open $command r]
    fconfigure $runWn(fileID) -blocking 0
    fileevent $runWn(fileID) readable wnRunWIENEvent
    
    tkwait variable runWn(event_done)
    ResetCursor
    xcSwapBuffers

    destroy $mw
    if { $oldgrab != {} } {
	catch { grab $oldgrab }
    }

    #
    # has an error occured during WIEN execution
    #
    if { $runWn(error) == 1 } {
	#if [catch {exec $command < $inp > $output_file} err]
	#C95Error $command $output_file $err
	cd $system(SCRDIR)
	unset runWn
	ResetCursor
	return 0
    } else {
	# if we have parallel execution, this will concatenate output files
	wnRunWIEN_ReadOutputFile $outf
    }
    unset runWn
    ResetCursor

    
    return 1
}
    

proc wnRunWIENEvent {} {
    global runWn
    
    if { ![eof $runWn(fileID)] } {
	append runWn(output) [gets $runWn(fileID)]\n
	xcDebug "$runWn(output)"
	return
    } else {
	xcDebug "In wnRunWIENEvent"
	set runWn(error) 0
	if { [catch {close $runWn(fileID)}] } {
	    set runWn(error) 1
	}
	
	#
	# if we have an error display tk_dialog
	#
	if { $runWn(error) } {
	    # is it really an error; check the error file !!!
	    if { [wnRun_IsERROR] } {
		wnRunWIENError
	    } else {
		set runWn(error) 0
	    }
	}
	
	set runWn(event_done) 1
    }
}


proc wnRunWIENError {} {
    global runWn
    
    set button [tk_dialog [WidgetName] ERROR "ERROR occure while executing WIEN2k command: $runWn(command)" error 0 OK Details]
    if { $button == 1} {
	set oldgrab [grab current]
	set t [xcToplevel [WidgetName] "WIEN2k ERROR" "ERROR" \
		. 100 100 1]
	catch { grab $t }

	set output [wnRunWIEN_ReadOutputFile $runWn(outf)]
	set text [DispText $t.f $output 80 20]
	set fric [lindex [$text yview] 1]
	$text yview moveto [expr 1.0 - $fric]
	$text config -state disabled
	
	proc wnRunWIENErrorClose {} {
	    global runWn
	    set runWn(error_done) 1
	}
	
	set f2 [frame $t.f2 -height 10]
	pack $f2 -side bottom -before $t.f -fill x
	set close [button $f2.cl -text "Close" \
		-command wnRunWIENErrorClose]
	pack $close -side left -expand 1 -ipadx 2 -ipady 2 -pady 10
	
	tkwait variable runWn(error_done)
	
	catch { grab release $t }
	destroy $t	
	if { $oldgrab != "" } {
	    catch { grab $oldgrab }
	}
    }
}


proc wnRunWEIN_IsParallel {} {
    global wn
    
    if { [info exists wn(parallel)] } {
	if { $wn(parallel) == 1 } {
	    return 1
	}
    }
    return 0
}

proc wnRunWIEN_DeleteOutputFile {outf} {
    global wn

    file delete $outf 
    if { [wnRunWEIN_IsParallel] } {
	file delete [glob -nocomplain ${outf}_*]
    }
}

proc wnRunWIEN_ReadOutputFile {outf} {
    global wn

    if { [wnRunWEIN_IsParallel] } {

	# it was a parallel run, concatenate output files
	
	set nproc [wnRunWIEN_NProc]
	puts stderr "wien-parallel: nproc=$nproc"
	
	if { $nproc > 0 } {
	    # check if all parallel output files exist and are larger > 0
	    set concate 1
	    for {set i 1} {$i<=$nproc} {incr i} {
		if { ! [file exists ${outf}_$i] } {
		    set concate 0
		    break
		} else {
		    if { [file size ${outf}_$i] == 0 } {
			set concate 0
			break
		    }
		    if { [file exists $outf] } {
			if { [file mtime $outf] > [file mtime ${outf}_$i] } {
			    # the $outf is newer than the ${outf}_$i, use $outf
			    set concate 0
			    break
			}
		    }
		}
	    }
	    if { $concate } {
		for {set i 1} {$i<=$nproc} {incr i} {
		    puts stderr "wien-parallel: proc=$i, output-file=${outf}_$i"
		    append output [ReadFile ${outf}_$i]
		}
		if { [info exists output] } {		
		    WriteFile $outf $output
		    return $output
		}
	    }
	}
    }

    puts stderr "wien: output-file=$outf"
    return [ReadFile $outf]    	
}


proc wnRunWIEN_NProc {} {
    global wn
    
    set pwdDir [pwd]

    cd $wn(dir)
       # try with this:
       catch {set result [exec grep -v init .processes | wc]}
    cd $pwdDir

    if { [info exists result] } {
	return [lindex $result 0]
    } else {
	return 0
    }    
}
