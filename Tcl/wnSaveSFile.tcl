#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnSaveSFile.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnSaveSFile {} {
    global system periodic

    if { ![xcIsActive render] && $periodic(dim) != 3 } {
	return
    }

    #if [file exists $system(SCRDIR)/xc_wienstruct.$system(PID)] {
    #	set file $system(SCRDIR)/xc_wienstruct.$system(PID)
    #}
    # t.k: I have disabled above "if", see the effects !!!
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    if { 0 } {
	;
    } else {
	if { [catch {exec $system(FORDIR)/savestruct \
			 $system(SCRDIR)/xc_struc.$system(PID) > \
			 $system(SCRDIR)/xc_wnstr.$system(PID)}] } {
	    tk_dialog [WidgetName] ERROR \
		"ERROR while saving WIEN2k STRUCT FILE" \
		error 0 OK
	    return
	}
    
	set pwd [pwd]
	cd $system(SCRDIR)
	set def [wn_nn_def xc_wnstr $system(PID)]
	WriteFile nn.def $def w
	set tmp $system(SCRDIR)/xc_tmp.$system(PID)
	WriteFile $tmp 2.0 w
	# later $system(FORDIR) should be replaced by WIEN2k's nn-path
	catch {exec $system(FORDIR)/nn nn.def < $tmp} error 
	switch -glob -- $error {
	    {*NN ENDS*} {
		# file was OK
		file copy -force xc_wnstr.$system(PID) xc_wnstr.struct
	    }
	    {*NN created*} {
		# new file was created and is writen to xc_wnstr.struct
		;
	    }
	    default {
		# an error has occured
		tk_dialog [WidgetName] ERROR \
			"ERROR while executing \"nn\" program;\n\
			ERORR Code: $error" \
			error 0 OK
		return
	    }
	}
	cd $pwd
	set file $system(SCRDIR)/xc_wnstr.struct
    }

    set filetypes {
	{{WIEN2k Struct File} {.struct} }
	{{All Files}          *         }
    }
    set sfile [tk_getSaveFile -initialdir $system(PWD) \
	    -title "Save WIEN2k Struct File" \
	    -defaultextension .struct \
	    -filetypes $filetypes \
	    -parent .]    
    
    # maybe Cancel button was pressed
    if { $sfile == {} } { return }
    if [catch {file copy -force $file $sfile}] {
	tk_dialog [WidgetName] ERROR \
	    "ERROR saving file $sfile" error 0 OK
	return
    }
    
    wm title . "XCrySDen: [file tail $sfile]"
    set xcMisc(titlefile) $sfile
    return
}
