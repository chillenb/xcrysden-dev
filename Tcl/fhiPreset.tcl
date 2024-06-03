#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/fhiPreset.tcl                                    #
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc fhiPreset {type file} {
    global system fhi

    if { $type != "inpini" && $type != "coord" } {
	ErrorDialog "invalid type of fhiPreset !!!"
	return 0
    }
    
    # just for any case
    cd $system(SCRDIR)
    if { $type == "inpini" } {
	##############################
	# execute CONVERTING program #
	##############################
	xcCatchExecReturn $system(BINDIR)/fhi_inpini2ftn34 $file getlist
	#if { [catch {exec $system(BINDIR)/fhi_inpini2ftn34 $file getlist}] } {
	#    ErrorDialog "error while executing \"fhi_inpini2ftn34 $file getlist\" program"
	#    return 0
	#}
    } else {
	# type == coord
	##############################
	# execute CONVERTING program #
	##############################
	xcCatchExecReturn $system(BINDIR)/fhi_coord2xcr $file getlist
	#if { [catch {exec $system(BINDIR)/fhi_coord2xcr $file getlist}] } {
	#    ErrorDialog "error while executing \"fhi_coord2xcr $file getlist\" program"
	#    return 0
	#}
    }
    
    # check if "fhi_species_name.list" was created ?!!
    if { ! [file exists fhi_species_name.list] } {
	ErrorDialog "\"FHIxxMD's\" file is probably bad !!!"
	return 0
    } 

    #
    # pop-up the toplevel with corresponding entries
    #

    set t [xcToplevel [WidgetName] "Atomic Numbers of Species" \
	    "Atomic Numbers" . 0 0 1]

    set f1 [frame $t.1]
    set f2 [frame $t.2]
    pack $f1 $f2 -side top -padx 10 -pady 10 -fill x

    frame $f1.title -class StressText
    label $f1.title.l -text "Specify Atomic Numbers of Species" \
	    -relief groove -bd 2
    $f1.title.l config -font [ModifyFontSize $f1.title.l 18]
    pack $f1.title -side top -expand 1
    pack $f1.title.l -ipadx 10 -ipady 3 -padx 10 -pady 10 -expand 1

    #
    # build up for ScrollEntries
    # 
    set list           [split [ReadFile fhi_species_name.list] \n]
    set fhi(n_species) [lindex $list 0]
    set species_names  [concat [lrange $list 1 end]]
    for {set i 0} {$i < $fhi(n_species)} {incr i} {
	set ii [expr $i + 1]
	set fhi(NAME,$ii) [lindex $species_names $i]
    }
    set labellist  [list "Species:" "Atomic Number:"]    
    set arraylist  [list NAME NAT]
    set arraytype  [list text posint]
    set buttonlist [list 1 \
	    [list {Periodic Table} scroll_ptableSelectNAT $t fhi NAT]]	    

    ScrollEntries \
	    $f1 \
	    $fhi(n_species) \
	    "Set Atomic Number for Species #" \
	    $labellist \
	    $arraylist \
	    $arraytype \
	    15 \
	    fhi \
	    $buttonlist \
	    3

    # OK and cancel button should be cerated as well
    set ok  [button $f2.ok -text OK -command [list fhiPresetOK $t]]    
    set can [button $f2.can -text Cancel -command [list fhiPresetCan $t]]
    pack $ok $can -side left -expand 1 -padx 10

    tkwait window $t

    return $fhi(status)
}
    
#
# maybe the program should found out if OK button was pressed by the
# presence of "fhi_species_nat.list" file !!!
#
proc fhiPresetOK t {
    global varlist foclist fhi system
    
    if ![check_var $varlist $foclist] {return 0}
    #
    # its seems OK; write the "fhi_species_nat.list" file and destroy 
    # the toplevel window
    set out "$fhi(n_species)\n"
    for {set i 1} {$i <= $fhi(n_species)} {incr i} {
	append out "$fhi(NAT,$i) "
    }
    ####################
    # just in any case #
    cd $system(SCRDIR)
    ####################
    WriteFile fhi_species_nat.list $out w

    destroy $t
    set fhi(status) 1
}
	    

proc fhiPresetCan t {
    global fhi system

    ####################
    # just in any case #
    cd $system(SCRDIR)
    ####################
    if { [file exists fhi_species_nat.list] } {
	file delete fhi_species_nat.list
    }
    destroy $t
    set fhi(status) 0
}
