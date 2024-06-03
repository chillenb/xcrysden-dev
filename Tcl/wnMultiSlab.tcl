#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnMultiSlab.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnMultiSlab {{change 0}} {
    global system multiSlab periodic

    if { $periodic(dim) != 2 && ![xcIsActive multislab] } {
	# multislab mode is not allowed
	xcDebug "WARNING: Multi-Slab mode not allowed"
	return
    }

    if { $change == 0 } {
	set answer [tk_messageBox -message "After performing Multi-Slab option a WIEN2k struct file will be created. This is usually the last AdvGeom option, since AdvGeom options will be disabled after multi-slab will be created\n\nDo You want to Continue?" \
		-type yesno -icon question]    
	if { $answer == "no" } { return }
    }

    # multislab options:
    # ------------------
    # SYMMINFO                 
    # MULTISLAB_NO_INVERSION   
    # MULTISLAB_INVERSION      
    #set inv_center [exec $system(FORDIR)/multislab SYMMINFO \
    #	    $system(SCRDIR)/xc_struc.$system(PID)]
    ########################################################
    # WARNING: multislab SYMMINFO is wrong at the moment !!!
    #          disablbe explotation of center of inversion
    set inv_center 0
    ########################################################

    set title "Create a Multi-Slab"
    if { $change == "change" } {
	set title "Change Multi-Slab Vacuum Thickness"
    }
    set t [xcToplevel [WidgetName] "Multi-Slab" $title . 100 0]
    set f1 [frame $t.1 -class RaisedFrame]
    set f2 [frame $t.2 -class RaisedFrame]
    pack $f1 $f2 -side top -expand 1 -fill both

    set multiSlab(inv) No
    if { $inv_center == 1 } {
	RadioBut $f1 "Structure has center of inversion !!!\nDo You want to exploit it ?" multiSlab(inv) top left 1 1 Yes No
	#label $f1.1 -text "Structure has center of inversion !!!\nDo You want to exploit it ?" -anchor w
	#radiobutton $f1.r1 -variable inv -text "Yes" -value yes
	#radiobutton $f1.r2 -variable inv -text "No"  -value no
    }
    set entry [Entries $f1 {{Vacuum Thickness:}} multiSlab(vacuum) 10]

    if ![info exists multiSlab(mb_angs/bohr)] {
	set multiSlab(mb_angs/bohr) Angstroms
    }
    set frame [string trimleft [filehead $entry] \{\}]
    set mb [menubutton $frame.mb \
	    -textvariable multiSlab(mb_angs/bohr) \
	    -menu $frame.mb.menu \
	    -indicatoron 1 \
	    -relief raised \
	    -width 9 \
	    -anchor w]
    pack $mb -side left -padx 2
    set menu [menu $mb.menu -relief raised -tearoff 0]
    $menu add command -label "Angstroms" \
	    -command [list set multiSlab(mb_angs/bohr) "Angstroms"]
    $menu add command -label "Bohrs" \
	    -command [list set multiSlab(mb_angs/bohr) "Bohrs"]
    
    set can [button $f2.can -text Cancel -command [list CancelProc $t]]
    set ok  [DefaultButton $f2.ok -text OK \
	    -command [list wnMultiSlabOK $t $entry]]
    pack $can $ok -padx 10 -pady 10 -side left -expand 1
    focus $entry
}

proc wnMultiSlabOK {t entry} {
    global multiSlab system geng periodic nzdir

    if ![check_var {{multiSlab(vacuum) real}} $entry] {
	return
    }
    
    if { $multiSlab(mb_angs/bohr) == "Bohrs" } {
	set multiSlab(vacuum) [Bohr2Angs $multiSlab(vacuum)]
    }
    if { $multiSlab(inv) == 0 } {
	set mode MULTISLAB_NO_INVERSION
    } else {
	set mode MULTISLAB_INVERSION
    }

    WriteFile $system(SCRDIR)/xc_tmp.$system(PID) $multiSlab(vacuum) w

    xcCatchExecReturn $system(FORDIR)/multislab $mode \
	$system(SCRDIR)/xc_struc.$system(PID) < \
	$system(SCRDIR)/xc_tmp.$system(PID) > \
	$system(SCRDIR)/xc_wnstr.$system(PID)
    #
    #if [catch {exec $system(FORDIR)/multislab $mode \
    #	    $system(SCRDIR)/xc_struc.$system(PID) < \
    #	    $system(SCRDIR)/xc_tmp.$system(PID) > \
    #	    $system(SCRDIR)/xc_wnstr.$system(PID)}] {
    #	tk_dialog [WidgetName] ERROR \
    #		"ERROR while executing \"multislab\" program" \
    #		error 0 OK
    #	return
    #}
    
    set pwd [pwd]
    cd $system(SCRDIR)
    set def [wn_nn_def xc_wnstr $system(PID)]
    WriteFile nn.def $def w
    # later $system(FORDIR) should be replaced by WIEN2k's nn-path
    set tmp $system(SCRDIR)/xc_tmp.$system(PID)
    WriteFile $tmp 2.0 w
    catch {exec $system(FORDIR)/nn nn.def < $tmp} error 
    switch -glob -- $error {
	{*NN ENDS*} {
	    # file was OK
	    file copy -force xc_wnstr.$system(PID) xc_wnstr.struct
	}
	{*NN created*} {
	    # new file was created
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

    xcCatchExecReturn $system(FORDIR)/str2xcr $system(SCRDIR)/xc_wnstr
    #if [catch {exec $system(FORDIR)/str2xcr $system(SCRDIR)/xc_wnstr}] {
    #	tk_dialog [WidgetName] ERROR \
    #		"ERROR while executing \"str2xcr\" program" \
    #		error 0 OK
    #	return
    #}
    file rename -force $system(SCRDIR)/xc_wnstr.xcr $system(SCRDIR)/xc_str2xcr.$system(PID)

    CancelProc $t

    #
    # WIEN2k struct file is in BOHRs, thatwhy xc_str2xcr.$$ is in BOHRs
    #
    # so far set periodic(group) -> 1
    set periodic(igroup) 1
    set periodic(dim)    3
    set nzdir            1
    set radio(cellmode)  prim
    set radio(hexamode)  "parapipedal"
    set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
    xcAppendState wien
    xcAppendState multislab
    foreach state {c95 openinput newinput} { xcDeleteState $state }
    xcUpdateState
    CellMode
    UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
    return
}


proc wn_nn_def {filehead pid} {
    return [format "%d,'%s','%s','%s',%d\n%d,'%s','%s','%s',%d\n%d,'%s','%s','%s',%d" \
	    66 $filehead.outputnn  UNKNOWN FORMATTED 0 \
	    20 $filehead.$pid      OLD     FORMATTED 0 \
	    21 $filehead.struct    UNKNOWN FORMATTED 0]
}
