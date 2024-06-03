#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnKPath.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnKPath {dir} {
    global Bz BzOK prop system geng periodic xcMisc

    # now struct file should be read and then BZ toplevel displayed
    # with XCrySDen's k-Path Selection: directory

    XCRYSDEN_Logo $dir
    set xcMisc(titlefile) $dir
    
    if [winfo exists .title] { 
	destroy .title 
    }


    ##########################################
    # in WIEN2k there it is always 3D system #
    ##########################################
    set periodic(dim) 3
    set species crystal

    cd $system(SCRDIR)

    set BzOK(wien_kpath) 1    
    set BzOK(wien_dir) [file tail $dir]

    #
    # check if $dir is a struct file or a directory !!!
    #
    if { [file isfile $dir] } {
	# it's a struct-file
	set filehead [file rootname $dir]
    } else {
	# it's a case directory
	set filehead $dir/[file tail $dir]
    }

    #
    # read the struct file and produce xc_struc.$$
    #
    xcDebug "Executing: $system(FORDIR)/str2xcr $filehead, PWD=[pwd]"
    if { [catch {exec $system(FORDIR)/str2xcr $filehead} errMsg] } {
	ErrorDialog "while executing \"str2xcr\" program." $errMsg
	CloseCase
    }
    # now copy $filehead.xcr to $system(SCRDIR)/xc_str2xcr.$$
    file rename -force ${filehead}.xcr $system(SCRDIR)/xc_str2xcr.$system(PID)

    #
    # WIEN2k struct file is in BOHRs, thatwhy xc_str2xcr.$$ is in BOHRs
    #
    set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
    xcAppendState wien
    set periodic(igroup) 1
    xcDebug "Going to GenGeom"
    GenGeom $geng(M1_PRIM) $geng(M2_CELL) $geng(M3_ARGUMENT) \
	    $periodic(igroup) \
	    1 1 1 $system(SCRDIR)/xc_struc.$system(PID)
    # determine the periodic(igroup)
    set fileID [open "$system(SCRDIR)/xc_struc.$system(PID)" r]
    GetDimGroup periodic(dim) periodic(igroup) $fileID
    close $fileID

    xcDebug "Going to xc_readXSF"
    ################################################
    xc_readXSF $system(SCRDIR)/xc_struc.$system(PID)
    ################################################
    
    #
    # render the k-Path selector window
    #
    set BzOK(done) 0
    Bz_MakeToplevel    

    #
    # OK, either the k-path was selected or a cancel button was pressed
    #
    if { !$BzOK(done) } {
	# Cancel button was pressed; EXIT
	exit_pr -silent
	exit
    }
    
    #----------------------------------
    # Prepare k_list-file:
    #----------------------------------
    # band structure script looks like:
    #
    # BAND
    # a-b-c-d-
    # 3 8 200 1 100 1 0
    # 0 0 0    4 4 0
    # 4 4 0    4 2 -2
    # 4 2 -2    5 3 0    
    #

    set script [split $prop(c95_BAND_script) \n]
    if { [lindex $script 0] == "CONVCELL" } {
	set script [lrange $script 1 end]
    }

    set labels [split [lindex $script 1] -]
    set idv    $BzOK(iss)
    set kp     [expr $prop(NLINE) + 1]
    set nk     $BzOK(nK)
  #  set emin   $BzOK(Emin)
  #  set emax   $BzOK(Emax)
    set emin   -8.0
    set emax    8.0
    for {set i 3} {$i < [expr 2 + $prop(NLINE)]} {incr i} {
	set line [lindex $script $i]
	set ii [expr $i - 3]
	set coor(x,$ii) [lindex $line 0]
	set coor(y,$ii) [lindex $line 1] 
	set coor(z,$ii) [lindex $line 2]
    }
    set line [lindex $script $i]
    set ii [expr $i - 3]
    set coor(x,$ii) [lindex $line 0]
    set coor(y,$ii) [lindex $line 1] 
    set coor(z,$ii) [lindex $line 2]
    incr ii
    set coor(x,$ii) [lindex $line 3]
    set coor(y,$ii) [lindex $line 4] 
    set coor(z,$ii) [lindex $line 5]
    
    # now we are ready to construct k_list file
    
    set    k_list "$kp $nk $idv $emin $emax\n"
 #   set    k_list "$kp $nk $idv \n"
    append k_list "$Bz(rendered)\n"
    set unk 0
    for {set i 0} {$i <= $ii} {incr i} {
	set l [lindex $labels $i]
	if { $l == {} } { 
	    set l X$unk 
	    incr unk
	}
	set l [format "'%-10s'" $l]
	append k_list "$coor(x,$i) $coor(y,$i) $coor(z,$i)     $l\n"
    }

    xcDebug -stderr "XCRYSDEN-to-WIEN k-list file:"
    xcDebug -stderr "-----------------------------"
    xcDebug -stderr $k_list

    # prompt for a name of final k-list file

    WriteFile $system(SCRDIR)/xc_klist.$system(PID) $k_list w
    set filetypes {
	{{WIEN2k K-List File} {.klist} }
	{{All Files}          *        }
    }
    set sfile [tk_getSaveFile -initialdir $system(PWD) \
	    -title "Save WIEN2k K-List File" \
	    -defaultextension .klist \
	    -filetypes $filetypes]    
    
    # maybe Cancel button was pressed
    if { $sfile == {} } { exit_pr -silent; exit }

    #
    # finaly execute a kPath program and save the klist into $sfile
    #

    xcDebug "exec $system(FORDIR)/kPath  \
	    $system(SCRDIR)/xc_struc.$system(PID) \
	    $system(SCRDIR)/xc_klist.$system(PID)"

    cd $system(SCRDIR)
    if { [catch {exec $system(FORDIR)/kPath  \
		     $system(SCRDIR)/xc_struc.$system(PID) \
		     $system(SCRDIR)/xc_klist.$system(PID) > \
		     $sfile} errMsg] } {
	ErrorDialog "while executing \"kPath\" program." $errMsg
    }

    #
    # provide some supporting information to check if the generated k-points are OK
    #
    global wnKP

    #set wnKP(npoi) $Bz($can,nselected)
    #set wnKP(M)    $BzOK(iss)
    #set wnKP(type) $type
    #/

    if { $wnKP(type) == "prim" } {
	set type "RECIPROCAL-PRIMITIVE"
    } else {
	set type "RECIPROCAL-CONVENTIONAL"
    }
	
    set    supportInfo "================================================================\n"
    append supportInfo "    * * * S U P P O R T I N G    I N F O R M A T I O N * * *    \n"
    append supportInfo "================================================================\n"
    append supportInfo "\n"
    #append supportInfo "BEWARE: the conversion of selected k-point coordinates to format\n"
    #append supportInfo "        suitable for the WIEN2k program is buggy.\n"
    #append supportInfo "\n"
    #append supportInfo "        !!! USE AT YOUR OWN RISK !!!\n"
    #append supportInfo "        !!! PLEASE CHECK BELOW DATA FOR CONSISTENCY !!!\n"
    #append supportInfo "\n"
    append supportInfo "Number of selected k-points: $wnKP(npoi)\n\n"
    append supportInfo "**ORIGINAL DATASET**\n"
    append supportInfo "--------------------\n"
    append supportInfo "THE FOLLOWING k-POINT WERE SELECTED:\n"
    append supportInfo "(crystal coordinates with respect to $type vectors)\n\n"
    for {set i 1} {$i <= $wnKP(npoi)} {incr i} {
	append supportInfo [format "   % 10.5f   % 10.5f   %10.5f     %s\n" \
				$wnKP(poi,$i,1) $wnKP(poi,$i,2) $wnKP(poi,$i,3)    $wnKP(label,$i)]
    }
    append supportInfo "\n**REPRODUCED(TRANSFORMED) DATASET for WIEN2k**\n"
    append supportInfo "----------------------------------------------\n"
    append supportInfo "IMPORTANT:\n"
    append supportInfo "Check if the data presented below are consistent with the data from ORIGINAL DATASET !!!\n\n"
    append supportInfo [ReadFile supportInfo-1.kpath]
    append supportInfo "\n"
    append supportInfo "Transformed selected k-points for WIEN2k:\n"
    append supportInfo "----------------------------------------\n"
    append supportInfo [ReadFile supportInfo-2.kpath]
    append supportInfo "\n"
    append supportInfo "The content of the generated WIEN2k k-list file is the following:\n"
    append supportInfo "----------------------------------------------------------------\n"
    append supportInfo [ReadFile $sfile]

    WriteFile $system(PWD)/supportInfo.kpath $supportInfo w
    if { [file exists supportInfo-1.kpath] } { file delete supportInfo-1.kpath }
    if { [file exists supportInfo-2.kpath] } { file delete supportInfo-2.kpath }

    set t [xcDisplayVarText $supportInfo "k-path: Supporting Information"]
    if { [winfo exists $t.f1.t] } {
	catch {$t.f1.t configure -height 30 -width 90}
	catch {$t.f1.t configure -state disabled}
    }
    tkwait window $t
    
    exit_pr -silent
    exit
}
