#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xsfOpen.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xsfUpdate {filedir {can .mesa}} {
    global select
    if { [info exists select(selection_mode)] && $select(selection_mode) } {
	tk_dialog [WidgetName] "WARNING !!!" "WARNING: structure can not be reloaded while in selection mode. First exit from the selection mode, then reload" warning 0 OK
	return
    }
    # we assume that if .iso or .iso2D windows exist we are in DataGrid mode
    if { [winfo exists .iso2D] || [winfo exists .iso] } {
        tk_dialog [WidgetName] "WARNING !!!" "WARNING: structure can not be reloaded while showing datagrids." warning 0 OK
        reloadButtonDisable
	return
    }    
    xsfOpen_ $filedir $can 1    
}

proc xsfAddBlankLine {filedir} {
    ##########################################
    # add a blank line at the end of XSF file
    set ch [open $filedir a]
    puts $ch ""
    flush $ch
    close $ch
    ##########################################
}

proc gunzipXSF {file} {
    global system

    xcDebug -debug "gunzipXSF: file = $file"

    set name [file tail $file]
    cd $system(SCRDIR)
    # maybe file is already locate in $system(SCRDIR)
    if { [file dirname $file] != $system(SCRDIR) && $file != $name } {
	#file copy -force $file $name
	fileCopy $file $name
	catch {exec -- chmod +rw $name}
    }
    ####################
    set gunzipName $name
    ####################

    if { [string match *.gz $name] } {
	xcDebug -debug "gunzipXSF: string matches"
	catch {exec -- gzip -fd $name}	
	set gunzipName [string trimright $name .gz]
	catch {exec -- chmod +rw $gunzipName}
	if { ![file exists $gunzipName] } {
	    ErrorDialog "error when gunzip-ing file $file"
	    uplevel 1 { return }
	}
    }
    set gunzipName $system(SCRDIR)/$gunzipName
    xcDebug -debug "gunzipXSF: gunziped-file = $gunzipName"

    xsfAddBlankLine $gunzipName
    return $gunzipName
}

#
# this proc is called form menu File->...->(Open XSF File)
#
proc xsfOpenMenu {can {dir {}}} {
    global fileselect system

    if { $dir == {} } {
        set dir $system(PWD)
    }
    set fileselect(path) [tk_getOpenFile -defaultextension .xsf \
                              -filetypes { 
                                  {{All Files}                  {.*}  }
                                  {{XSF Files}                  {.xsf}}
                                  {{GZipped XSF Files}          {.xsf.gz}}
                                  {{Animated XSF Files}         {.axsf}}
                                  {{GZipped Animated XSF Files} {.axsf.gz}}
                              } -initialdir $dir -title "Open XSF File"]
    if { $fileselect(path) == "" } {
        return
    }
    reloadRegister $fileselect(path) XSF
    xsfOpen $fileselect(path) $can
}

# call this proc when the XSF file to open is a result of parsing from other file formats
proc xsfOpenParsed {filedir {can .mesa} {update 0}} {
    xsfOpen_ $filedir $can $update
}

# call this proc only for direct opening of the XSF file (for parsed XSF files, call xsfOpenParsed instead)
proc xsfOpen {filedir {can .mesa}} {

    UpdateWMTitle $filedir
    xsfOpen_ $filedir $can
}

proc xsfOpen_ {filedir can {update 0}} {
    global system geng xcMisc periodic radio sInfo periodic xsfAnim working_XSF_file

    set filedir [gunzipXSF $filedir]

    if { ! [info exists xcMisc(reduce_to)] } {
	set xcMisc(reduce_to) {}
    }
    eval {exec $system(BINDIR)/xsf2xsf $filedir $filedir.raw} $xcMisc(reduce_to)
    set filedir $filedir.raw

    #
    # check if xsf is AXSF (Animation XSF)
    #
    set fID  [open $filedir]
    set l1st [gets $fID]
    close $fID
    if { [string match "*ANIMSTEP*" $l1st] } {
	# yes we have AXSF file
	if { $xsfAnim(not_anim) != 1 } {
            puts stderr "*** xsfAnimInit $filedir .mesa"
	    xsfAnimInit $filedir .mesa $update
	    return
	}
    }

    #
    # the working XSF file
    set working_XSF_file $filedir

    if { $update != 0 } {
	GetDimGroupXSF periodic(dim) periodic(igroup) $filedir
	puts stderr "*** update XSF:  dim = $periodic(dim)"
	if { $periodic(dim) < 3 && $radio(cellmode) == "conv" } {
	    set radio(cellmode) prim
	}
	if { $periodic(dim) > 1 } {
	    CellMode 1
	} else {
	    UpdateStruct .mesa $filedir
	}
	xcUpdateState
        return
    }
    
    #
    # now open the file
    #
    ResetDispModes
    if { [catch {xc_openstr xcr $filedir $can PL}] } {
	ErrorDialog "An Error occured, while reading XSF file $filedir"
	CloseCase
	return
    }

    #
    # set the proper state
    #
    Get_sInfoArray
    DisplayDefaultMode
    xcAppendState render
    xcUpdateState
    set periodic(dim) $sInfo(dim)
    
    #
    # if the structure is periodic update it according to the "state"
    #
    if { $periodic(dim) == 0 } {
        # save $filedir to xc_struc.$system(PID) so that SaveXSFStruct will work
        file copy -force $filedir $system(SCRDIR)/xc_struc.$system(PID)
    } else {
	set geng(M3_ARGUMENT) [GetGengM3Arg ANGS]
	
	# the state should be set as well
	cd $system(SCRDIR)
	GenGeom $geng(M1_INFO) $geng(M2_CELL) $geng(M3_ARGUMENT) \
	    1 1 1 1 xc_gengeom.$system(PID)
	set fileID [open "$system(SCRDIR)/xc_gengeom.$system(PID)" r]
	GetDimGroup periodic(dim) periodic(igroup) $fileID
	close $fileID

	#<new>:  t.k. Mon Jan 27 15:24:44 CET 2003
	CellMode 1
	#/

	Get_sInfoArray
	xcUpdateState
    }

    CrysFrames

    #DEBUG:
    #xc_setatomlabel .mesa 1 "1st atom" "" {1.0 1.0 0.0 } {0.0 0.0 0.5}	    
    #/
    return
}


########################################################################
# check if XSF uses atomic symbols instead of atomic numbers
# if it does convert to atomic numbers ...
proc xsfAtmSym2Nat {filedir} {
    global system

    
    file copy -force $filedir $system(SCRDIR)/tmp.xsf

    # open existing XSF file ...
    set fID     [open $system(SCRDIR)/tmp.xsf r]

    # create a new file ...
    set newFile $filedir
    set newID   [open $newFile w]

    # search for ATOMS, PRIMCOORD, CONVCOORD section and replace ...
    set read_atoms 0
    while { ! [eof $fID] } {	
	gets $fID line
	
	if { [regexp -- "^ *ATOMS" $line] } {
	    set read_atoms 1
	    
	} elseif { [regexp -- "^ *PRIMCOORD|^ *CONVCOORD" $line] } {
	    puts $newID $line
	    gets $fID line
	    set read_atoms 1
	    
	} elseif { $read_atoms == 1 } {
	    set len [llength $line]
	    if { $len != 4 && $len != 7 } {
		set read_atoms 0
	    } else {
		set line [_xsfAtmSym2Nat $line]
	    } 
	}
	puts $newID $line
    }
    close $fID
    close $newID
    return $newFile
}


proc _xsfAtmSym2Nat {line} {
    set f1 [lindex $line 0]
    if { ! [string is integer $f1] } {
	set f1 [AnameExt2Nat $f1]
    }
    return [concat $f1 [lrange $line 1 end]]
}


