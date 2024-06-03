#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/wnDensity.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc wnDensity dir {
    global system geng

    set filehead [file tail $dir]

    if { ! [file isdirectory $dir] } {
	# $dir is not a directory; ERROR
	ErrorDialog "couldn't find directory \"$dir\" while trying to render wnDensity"
	return 0
    }

    #cd $system(SCRDIR)
    #exec cp $dir/$filehead.struct $filehead.struct
    #if [catch {exec $system(FORDIR)/str2xcr $filehead}] {
    #	 tk_dialog [WidgetName] ERROR \
    #		 "ERROR while executing \"str2xcr\" program" \
    #		 error 0 OK
    #	 return
    #}
    #
    #set geng(M3_ARGUMENT) [GetGengM3Arg BOHR 95]
    ##set xcMisc(titlefile) $filehead
    #xcCatchExecReturn $system(FORDIR)/gengeom $geng(M1_PRIM) $geng(M2_CELL) \
    #	     $geng(M3_ARGUMENT) 1 1 1 1 \
    #	     $system(SCRDIR)/xc_struc.$system(PID) \
    #	     $system(SCRDIR)/$filehead.xcr

    #
    # *.output5 & *.rho file must exists!!!   
    #
    foreach file [list $dir/$filehead.output5 $dir/$filehead.rho] {
	if { ! [file exists $file] } {
	    # $dir is not a directory; ERROR
	    ErrorDialog "couldn't find file \"$file\" while trying to render wnDensity"
	    return 0
	}
    }
    
    #xcCatchExecReturn $system(awk) -f $system(AWKDIR)/getOV.awk $dir/$filehead.output5
    
    if { [catch {set vec [exec $system(awk) -f $system(AWKDIR)/getOV.awk \
    			      $dir/$filehead.output5]}] } {
    	ErrorDialog "error while executing \"getOV.awk\" program"
    	return
    }
    
    #
    # since wnOpenSFile routine was called before wnDensity, the struct file 
    # already converted to $system(SCRDIR)/xc_struc.$system(PID)
    set geom  [ReadFile $system(SCRDIR)/xc_struc.$system(PID)]
    update
    set fID   [open $dir/$filehead.rho r] 
    set outID [open $system(SCRDIR)/xc_struc.$system(PID) w]

    ################################
    # make XCRYSDEN STRUCTURE FILE #
    ################################
    puts $outID $geom
    puts $outID "BEGIN_BLOCK_DATAGRID2D"
    puts $outID "Density_by_WIEN____FILE_generated_by_XCrySDen"
    puts $outID "DATAGRID_2D_Total_Density"
    #
    # WARNING: WIEN is writing grid as (den(i,j),j=1,ny),i=1,nx, BUT
    #          XCrySDen is reading it as (den(i,j),i=1,nx),j=1,ny
    #  FIXING: interchange vec-X with vec-Y and interchange nx with ny
    #
    set i 0
    foreach line [split [read -nonewline $fID] \n] {
	if { $i == 0 } {
	    puts $outID [concat [lindex $line 1] [lindex $line 0]]
	    puts $outID $vec
	} else {
	    # read the rest of the file and write to datagrid
	    puts $outID $line
	}
	incr i
    }
    puts $outID "END_DATAGRID_2D"
    puts $outID "END_BLOCK_DATAGRID2D"
    flush $outID
    close $outID
    close $fID

    #
    # now execute xc_updatestr and ...
    #
    UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
    Get_sInfoArray
    xc_isodatagrid info    
    xc_isodatagrid 0 0 1.0
    DataGrid2Isosurf
    IsoControl2D
}


proc wnDensity2D_or_3D dir {
    global prop wn

    set wn(filehead) [file tail $dir]
    
    if { ! [file isdirectory $dir] } {
	# $dir is not a directory; ERROR
	ErrorDialog "couldn't find directory \"$dir\" while trying to render wnDensity"
	exit_pr
    }
    set wn(dir) $dir

    # check if *.clmval[up|dn] exists
    #if { [glob -nocomplain -- $dir/$filehead.clm*] == {} } {
    #	 tk_dialog [WidgetName] ERROR \
    #		 "ERROR: no \"$filehead.clm*\" files were found !!!" error 0 OK
    #	 exit_pr
    #	 return 0
    #}

    set wn(lapw5) {x lapw5}
    if { $wn(complex) } { set wn(lapw5) {x lapw5 -c} }

    #
    # this is to prevent inconveniences; otherwise prop is used for CRYSTALXX
    #
    set prop(type_of_run)   RHF

    set wn(done) 0
    SetIsoSurfArray
    #
    # ask user what to do: 2D/3D
    #
    set t [xcToplevel [WidgetName] "wnDensity: 2D/3D plot" \
	    "wnDensity" . 25 50 1]
    catch { grab $t }
    set f [frame $t.f -relief raised -bd 2 -class StressText]
    pack $f -fill both -expand 1
    label $f.l -text "Choose 2D/3D Density Plot" -relief groove -bd 2
    set d2 [button $f.d2 -text "2D Density Plot" \
	    -command {global wn; set wn(plot_dim) 2D}]
    set d3 [button $f.d3 -text "3D Density Plot" \
	    -command {global wn; set wn(plot_dim) 3D}]
    pack $f.l $d2 $d3 -side top -padx 20 -pady 10 -ipadx 10 -ipady 10

    tkwait variable wn(plot_dim)

    catch { grab release $t }
    destroy $t

    SetIsoGridSpace \
	    "wnDensity $wn(plot_dim) Map - Grid Specification" "Grid" \
	    {WIEN} $wn(plot_dim)
    tkwait variable wn(done)
}


proc wnMakeIn5_2D3D {pov {dim 2D}} {
    global wn isosurf system
    
    # 2D ... pov stands for points
    # 3D ... pov stands for origin & vectors a,b,c

    set t [WidgetName]
    wnDensityFlags $t
    tkwait window $t

    set in5_file $wn(filehead).in5
    if $wn(complex) {
	set in5_file ${in5_file}c
    }

    if { $dim == "2D" } {
	# point 1 will be taken for origin
	# point 2 as X-end of plot
	# point 3 as Y-end of plot
	#
	#   0--3   Y/0
	#   |  |    |
	#   1--2   O/1--X/2
	#
	# points are in angstroms!!!
    
	# X-vector
	set v(0,0) [expr [lindex $pov 6] - [lindex $pov 3]]
	set v(0,1) [expr [lindex $pov 7] - [lindex $pov 4]]
	set v(0,2) [expr [lindex $pov 8] - [lindex $pov 5]]
	# Y-vector
	set v(1,0) [expr [lindex $pov 0] - [lindex $pov 3]]
	set v(1,1) [expr [lindex $pov 1] - [lindex $pov 4]]
	set v(1,2) [expr [lindex $pov 2] - [lindex $pov 5]]

    } elseif { $dim == "3D" } {
	# pov == origin
	#        vec1
	#        vec2
	#        vec3
	set o(0) [lindex $pov 0]
	set o(1) [lindex $pov 1]
	set o(2) [lindex $pov 2]

	xcDebug -debug "Pov:: $pov"
	xcDebug -debug "Origin:: $o(0)  $o(1)  $o(2)"

	for {set i 0} {$i < 3} {incr i} {
	    for {set j 0} {$j <3} {incr j} {
		set index [expr 3 + $i*3 + $j]
		set v($i,$j) [lindex $pov $index]
	    }
	    xcDebug -debug "Vector:: $v($i,0)  $v($i,1)  $v($i,2)"
	}
    } else {
	xcDebug -stderr "************\nunknown dim $dim, must be 2D or 3D\n************"
	return
    }
    
    set distX [expr sqrt($v(0,0)*$v(0,0) + $v(0,1)*$v(0,1) + $v(0,2)*$v(0,2))]
    set distY [expr sqrt($v(1,0)*$v(1,0) + $v(1,1)*$v(1,1) + $v(1,2)*$v(1,2))]

    if { $isosurf(res_type) == "points" } {
	# let this be number of points in X-dir
	set wn(nx) $isosurf(resol_poi)
	set res [expr $distX / $wn(nx)]
    } else {
	if { $isosurf(mb_angs/bohr) == "Angstroms" } {
	    set res [expr double($isosurf(resol_ang))]
	} else {
	    set res [Bohr2Angs $isosurf(resol_ang)]
	}
	set wn(nx) [expr round( $distX / $res )]
    }    
    set wn(ny) [expr round( $distY / $res )]

    if { $dim == "3D" } {
	set distZ  [expr sqrt($v(2,0)*$v(2,0) + $v(2,1)*$v(2,1) + $v(2,2)*$v(2,2))]
	set wn(nz) [expr round( $distZ / $res )]
    }

    # find out the crystal lattice type from *.struct file??? This is save !!!
    # WIENXX convention:
    # P     -> primitive 
    # F,B,C -> conventional
    # H     -> hexagonal primitive
    # R     -> hexagonal primitive (that is conventional for R)
    switch -glob -- $wn(lattice_type) {
	P*   { set vtype prim }
	F*   { set vtype conv }
	B*   { set vtype conv }
	CXY* { set vtype conv }
	CYZ* { set vtype conv }
	CXZ* { set vtype conv }
	R*   { set vtype prim }
	H*   { set vtype prim }
    }

    if { $dim == "2D" } {
	set frpoi [xc_isospacesel .mesa fractcoor -vtype $vtype]
	
	set m 1.0e7
	set im [expr round($m)]
	for {set i 0} {$i < 3} {incr i} {
	    set y($i) [expr round([lindex $frpoi [expr 0 + $i]] * $m)]
	    set o($i) [expr round([lindex $frpoi [expr 3 + $i]] * $m)]
	    set x($i) [expr round([lindex $frpoi [expr 6 + $i]] * $m)]
	}
	append in5 "[format %12d%12d%12d%14d $o(0) $o(1) $o(2) $im]\n"
	append in5 "[format %12d%12d%12d%14d $x(0) $x(1) $x(2) $im]\n"
	append in5 "[format %12d%12d%12d%14d $y(0) $y(1) $y(2) $im]\n"
	append in5 "$wn(xnsh) $wn(ynsh) $wn(znsh)\n"
	append in5 "$wn(nx) $wn(ny)\n"
	append in5 [format "%-4s%-4s\n" $wn(den_flag) $wn(spin_flag)]
	append in5 [format "%-4s%-4s%-4s\n" \
		$wn(unit_flag) $wn(norm_flag) $wn(debu_flag)]
	append in5 "NONORTHO\n"

	#########################################################
	# now give user the opportunity to check the *.in5 file #
	#########################################################
	set t [xcToplevel $t "wnDensity: check the $in5_file file" \
		"wnDensity" . 25 50 1]
	set f1 [frame $t.1 -relief raised -bd 2]
	set f2 [frame $t.2 -relief raised -bd 2]
	pack $f1 -side top -expand 1 -fill both
	pack $f2 -side top -expand 1 -fill x
	
	set text [DispText $f1.f $in5 60 15]
	
	set ok [DefaultButton $f2.ok -text OK \
		-command [list wnMakeIn5_2D_OK $t $text $in5_file]]
	pack $ok -padx 5 -pady 10 -expand 1  

    } elseif { $dim == "3D" } {

	# convert from Cartesian to Fractional

	set orig [xc_fractcoor -ctype $vtype -coor [list $o(0) $o(1) $o(2)]]
	set ori(0) [lindex $orig 0]
	set ori(1) [lindex $orig 1]
	set ori(2) [lindex $orig 2]
	for {set i 0} {$i < 3} {incr i} {
	    set lvec [xc_fractcoor -ctype $vtype -coor [list $v($i,0) $v($i,1) $v($i,2)]]
	    for {set j 0} {$j < 3} {incr j} {
		set vec($i,$j) [lindex $lvec $j]
	    }
	}	

	set t [xcToplevel [WidgetName] "Calculating 3D grid of points" \
		"Calculating" . 100 250 1]			
	set f [frame $t.f -relief raised -bd 2 -width 400]
	pack $f -expand 1
	set wn(calc_slice_text) "\
	\n\
	***  WIEN is calculating 3D grid of points  ***\n\
	\n\
	Number of 2D slices to calculate:   $wn(nz)\n\
	Progress:   0 / $wn(nz)\n"
	set m [label $f.m \
		-textvariable wn(calc_slice_text) \
		-justify left \
		-relief groove -bd 2]
	set ff [frame $f.f]
	pack $m $ff -side top -expand 1 -ipadx 10 -ipady 10 -padx 10 -pady 10
	set w [expr round(300 / $wn(nz))]
	for {set i 1} {$i <= $wn(nz)} {incr i} {
	    frame $ff.$i -width $w -height 8m -relief raised -bd 2
	    grid $ff.$i -row 0 -column $i
	}

	file delete $system(SCRDIR)/xc_rho.$system(PID)

	#############################
	# vec(-1,$) is DUMMY & must be (0.0, 0.0, 0.0)!!!!!!!!
	set vec(-1,0) 0.0
	set vec(-1,1) 0.0
	set vec(-1,2) 0.0
	for {set i 0} {$i < $wn(nz)} {incr i} {	
	    ############################
	    # calculate $wn(nz) slices #
	    ############################
	    if { $wn(nz) > 1 } {
		set f [expr $i / [expr $wn(nz) - 1.0]]
	    } else {
		set f 1.0
	    }
	    for {set j -1} {$j < 2} {incr j} {	    
		set jj [expr $j + 1]
		set v($i,$jj,0) \
			[expr $ori(0) + $vec($j,0) + $f * $vec(2,0)]
		set v($i,$jj,1) \
			[expr $ori(1) + $vec($j,1) + $f * $vec(2,1)]
		set v($i,$jj,2) \
			[expr $ori(2) + $vec($j,2) + $f * $vec(2,2)]	    
		xcDebug "Factor:: $f\nEndPoints:: $v($i,$jj,0)   $v($i,$jj,1)   $v($i,$jj,2)\n"
	    }
		
	    set m 1.0e7
	    set im [expr round($m)]
	    for {set j 0} {$j < 3} {incr j} {
		set io($j) [expr round( $v($i,0,$j) * $m)]
		set ix($j) [expr round( $v($i,1,$j) * $m)]
		set iy($j) [expr round( $v($i,2,$j) * $m)]
	    }
	    set in5 {}
	    append in5 "[format %12d%12d%12d%14d $io(0) $io(1) $io(2) $im]\n"
	    append in5 "[format %12d%12d%12d%14d $ix(0) $ix(1) $ix(2) $im]\n"
	    append in5 "[format %12d%12d%12d%14d $iy(0) $iy(1) $iy(2) $im]\n"
	    append in5 "$wn(xnsh) $wn(ynsh) $wn(znsh)\n"
	    append in5 "$wn(nx) $wn(ny)\n"
	    append in5 [format "%-4s%-4s\n" $wn(den_flag) $wn(spin_flag)]
	    append in5 [format "%-4s%-4s%-4s\n" \
		    $wn(unit_flag) $wn(norm_flag) $wn(debu_flag)]
	    append in5 "NONORTHO\n"
	    
	    xcDebug "****************************"
	    xcDebug "Calculating slice::::\n$in5"
	    xcDebug "****************************"
	    if { $i > 0 } {
		file copy -force $in5_file ${in5_file}.$i
	    }
	    WriteFile $wn(dir)/$in5_file $in5 w
	    
	    #
	    # now run WIEN2k
	    #
	    if { ![wnRunWIEN $wn(exe_flag) {WIEN program is calculating the grid of points now. It can take some time, so PLEASE WAIT!!!} \
		       $wn(dir)/$wn(filehead).output5] } {
		# exit, since an error occured
		destroy $t
		exit_pr
		return
	    }

	    # now append the *.rho file to $system(SCRDIR)/xc_rho.$$
	    wnAppendRhoFile

	    # update "progress" window
	    set ii [expr $i + 1]
	    set wn(calc_slice_text) "\
		    \n\
		    ***  WIEN is calculating 3D grid of points  ***\n\
		    \n\
		    Number of 2D slices to calculate:   $wn(nz)\n\
		    Progress:   $ii / $wn(nz)\n"	    
	    $ff.$ii config -bg "#5f5"
	    update    
	}
	# destro progress widget
	destroy $t

	################################
	# make XCRYSDEN STRUCTURE FILE #
	################################
	# since wnOpenSFile routine was called before wnDensity, 
	# the struct file already converted to 
	# $system(SCRDIR)/xc_struc.$system(PID)
	set geom  [ReadFile $system(SCRDIR)/xc_struc.$system(PID)]
	update
	set fID   [open $system(SCRDIR)/xc_rho.$system(PID) r] 
	set outID [open $system(SCRDIR)/xc_struc.$system(PID) w]

	#
	# WARNING: WIEN is writing grid as (den(i,j),j=1,ny),i=1,nx, BUT
	#          XCrySDen is reading it as (den(i,j),i=1,nx),j=1,ny
	#  FIXING: interchange vec-X with vec-Y and interchange nx with ny
	# 	
	puts $outID $geom
	puts $outID "BEGIN_BLOCK_DATAGRID3D"
	puts $outID "Density_by_WIEN____FILE_generated_by_XCrySDen"
	puts $outID "DATAGRID_3D_Total_Density"
	puts $outID "$wn(ny) $wn(nx) $wn(nz)"
	puts $outID "$o(0)  $o(1)  $o(2)"
	puts $outID "$v(1,0)  $v(1,1)  $v(1,2)"
	puts $outID "$v(0,0)  $v(0,1)  $v(0,2)"
	puts $outID "$v(2,0)  $v(2,1)  $v(2,2)"

	foreach line [split [read -nonewline $fID] \n] {
	    # read the rho file and write to struct file
	    puts $outID $line
	}
	puts $outID "END_DATAGRID_3D"
	puts $outID "END_BLOCK_DATAGRID3D"
	flush $outID
	close $outID
	close $fID
	# rho file is no more needed, delete it
	file delete $system(SCRDIR)/xc_rho.$system(PID)

	#
	# now execute xc_updatestr and render isosurface
	#
	UpdateStruct .mesa $system(SCRDIR)/xc_struc.$system(PID)
	Get_sInfoArray
	xc_isodatagrid info    
	xc_isodatagrid 0 0 1.0
	DataGrid2Isosurf
	IsoControl	
    }
}


proc wnAppendRhoFile {} {
    global wn system
    
    set out [ReadFile -nonewline $wn(dir)/$wn(filehead).rho]
    set il 0
    foreach line [split $out \n] {
	if { $il > 0 } {
	    append output "$line\n"
	}
	incr il
    }
    #if { ! [file exists $system(SCRDIR)/xc_rho.$system(PID)] } {
    #	#exec touch $system(SCRDIR)/xc_rho.$system(PID)
    #	WriteFile $system(SCRDIR)/xc_rho.$system(PID) {}
    #}
    WriteFile $system(SCRDIR)/xc_rho.$system(PID) $output a
}
    

proc wnMakeIn5_2D_OK {t textw in5_file} {
    global wn
    # t     ... toplevel
    # textw ... text widget

    set text [$textw get 1.0 end]
    xcDebug "${in5_file}::\n$text"
    WriteFile $wn(dir)/$in5_file $text w
    CancelProc $t

    #
    # now run WIEN2k
    #
    if ![wnRunWIEN $wn(exe_flag) {WIEN program is calculating the grid of points now. It can take some time, so PLEASE WAIT!!!} \
	    $wn(dir)/$wn(filehead).output5] {
	# exit, since an error occured
	exit_pr
	return
    }
    wnDensity $wn(dir)
}


proc wnDensityFlags t {
    global wn

    # some parameters must be queried
    xcToplevel $t "wnDensity - Specify Flags" "wnDensity" . 25 50 1
    
    set f1 [frame $t.f1 -relief raised -bd 2]
    set f2 [frame $t.f2 -relief raised -bd 2]
    pack $f1 $f2 -side top -expand 1 -fill both

    set wn(exe_flag) $wn(lapw5)
    if { $wn(spin_polarized) == 1 } {
	set wn(exe_flag)  [concat $wn(lapw5) -up]
    }
    set wn(xnsh) 3
    set wn(ynsh) 2
    set wn(znsh) 3
    set wn(den_flag)  RHO
    set wn(spin_flag) ADD
    set wn(unit_flag) ATU
    set wn(norm_flag) VAL
    set wn(debu_flag) NODEBUG

    frame $f1.f
    FillEntries $f1.f {"# of X shells:" "# of Y shells:" "# of Z shells:"} \
	    {wn(xnsh) wn(ynsh) wn(znsh)} 14 10

    xcDebug -debug "wn(spin_polarized) = $wn(spin_polarized)"
    if $wn(spin_polarized) { 
	set mt1 [concat $wn(lapw5) -up]  
	set mt2 [concat $wn(lapw5) -dn]
	set mb_exeflag [xcMenuButton $f1 -labeltext "Execution Flag:" \
		-labelwidth 14 \
		-textvariable wn(exe_flag) \
		-menu [list \
		$mt1  {set wn(exe_flag) [concat $wn(lapw5) -up]} \
		$mt2  {set wn(exe_flag) [concat $wn(lapw5) -dn]}]]
    }
    
    set mb_den [xcMenuButton $f1 -labeltext "Density Flag:" -labelwidth 14 \
	    -textvariable wn(den_flag) \
	    -menu { 
	"SCF Density  (RHO)"                      {set wn(den_flag) RHO}
	"Difference Density  (DIFF)"              {set wn(den_flag) DIFF}
	"Superposition of Atomic Density  (OVER)" {set wn(den_flag) OVER}
    }]

    set mb_spin [xcMenuButton $f1 -labeltext "Spin Flag:" -labelwidth 14 \
	    -textvariable wn(spin_flag) \
	    -menu {
	"Total Density  (ADD)"  {set wn(spin_flag) ADD}
	"Spin Density  (SUB)"   {set wn(spin_flag) SUB}
    }]
    if !$wn(spin_polarized) {
	xcDisableAll $mb_spin
    }
	
    set mb_unit [xcMenuButton $f1 -labeltext "Unit Flag:" -labelwidth 14 \
	    -textvariable wn(unit_flag) \
	    -menu {
	"e / a.u.^3  (ATU)"  {set wn(unit_flag) ATU}
	"e / A^3  (ANG)"     {set wn(unit_flag) ANG}
    }]

    set mb_norm [xcMenuButton $f1 -labeltext "Norm. Flag:" -labelwidth 14 \
	    -textvariable wn(norm_flag) \
	    -menu {
	"Valence Density  (VAL)"  {set wn(norm_flag) VAL}
	"Total Density  (TOT)"     {set wn(norm_flag) TOT}
    }]

    set mb_debu [xcMenuButton $f1 -labeltext "Debug Flag:" -labelwidth 14 \
	    -textvariable wn(debu_flag) \
	    -menu {
	"No Debugging  (NODEBUG)"  {set wn(debu_flag) NODEBUG}
	"Debugging  (DEBUG)"       {set wn(debu_flag) DEBUG}
    }]

    grid $f1.f       -column 0 -row 0 -sticky w -padx 10 -pady 5	
    if $wn(spin_polarized) {
	grid $mb_exeflag -column 0 -row 1 -sticky w -padx 10 -pady 5	
    }
    grid $mb_den     -column 0 -row 2 -sticky w -padx 10 -pady 5
    grid $mb_spin    -column 0 -row 3 -sticky w -padx 10 -pady 5
    grid $mb_unit    -column 0 -row 4 -sticky w -padx 10 -pady 5
    grid $mb_norm    -column 0 -row 5 -sticky w -padx 10 -pady 5
    grid $mb_debu    -column 0 -row 6 -sticky w -padx 10 -pady 5     

    # let the user determine also execution type
    set	ok [DefaultButton $f2.ok -text OK \
	    -command [list wnDensityFlagsOK $t]]
    set can [button $f2.can -text Cancel -command exit_pr]

    pack $can $ok -side left -padx 5 -pady 10 -expand 1
}


proc wnDensityFlagsOK t {
    global wn    
    #
    # so far user will have to choose correct flags, 
    # no checking will be performed
    #
    CancelProc $t
    set wn(done1) 1; #this is not needed !!!
}


# tmp proc::
#proc PrintDecIdv list {
#
#    set idv [lindex $list 3]
#
#    set i0 [expr double([lindex $list 0]) / double($idv)]
#    set i1 [expr double([lindex $list 1]) / double($idv)]
#    set i2 [expr double([lindex $list 2]) / double($idv)]
#    xcDebug "$i0 $i1 $i2"
#}


#lappend auto_path "/home/tone/prog/XCrys/XCrySDen0.1"
#button .b -text TONE
#entry  .e -text TONE
#set xcFonts(normal)       [lindex [.b configure -font] end]
#set xcFonts(normal_entry) [lindex [.e configure -font] end]
#set xcFonts(small)        [ModifyFontSize .b 10 \
#	     {-family helvetica -slant r -weight bold}]
#set xcFonts(small_entry)  [ModifyFontSize .e 10 \
#	     {-family helvetica -slant r -weight normal}]
#destroy .b .e
#
#wnDensityFlags [WidgetName]
