#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/read25.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc ReadFTN25 {what fileID gID {ns 1} {ng 1}} {
    global grafdata ftn25 system
    # ns --> N_segment
    # ng --> N_graf

    xcDebug "what: $what"
    # this is used for projected band structure
    if { $what == "BAND" } {
	set nproj $ns
	set ns    1
    }

    xcDebug "ReadFTN25 $what $fileID $gID $ns $ng"
    set is0 $ns
    if { $what == "BAND" } {
	set is0 1
    }

    for {set is $is0} {$is <= $ns} {incr is} {
	# CRYSTAL-95::	
	# read FORMAT: A4,A4,2X,6I3,1X,A36,\n)
	#       start:  0  4  8  10
	# CRYSTAL98/03::
	# read FORMAT: A3,I1,A4,2I4,1P,(3E12.5)
	#       start:  0  4  5   8 16  17
	#              1P, 2E12.5
	#              6I3,1X,A36
	# -%-2DOSS   1 202 0.00000E+00 4.54430E+00-1.84013E-01
        # 0.00000E+00-9.05903E+02
        # 1 36  0  0  0  0            **********

	#################
	gets $fileID Line
	#################

	set str [string range $Line 4 7]
	xcDebug "str:: $str"
	if { $what != $str } {
	    tk_dialog .dialog WARNING \
		    "Instead of reading \"$what\" word, \n\
		    a \"$str\" word has been read!!!\n\
		    PLEASE report to Tone.Kokalj@ijs.si" error 0 OK
	    return
	}

	if { $system(c95_version) == "95" } {
	    for {set i 1} {$i <= 6} {incr i} {
		set f [expr 10 + ($i - 1)*3]
		set l [expr $f + 2]
		set ftn25(I$i) [string range $Line $f $l]
		xcDebug "ftn25(I$i): $ftn25(I$i)"
	    }
	    
	    #read FORMAT: 2I4, 5F10.6,\n

	    #################
	    gets $fileID Line
	    #################

	    set ftn25(NROW) [string range $Line 0 3]
	    set ftn25(NCOL) [string range $Line 4 7]
	    xcDebug "(NROW,NCOL): $ftn25(NROW) $ftn25(NCOL)"

	    for {set i 1} {$i <= 5} {incr i} {
		set f [expr 8 + ($i - 1)*10]
		set l [expr $f + 9]
		set ftn25(X$i) [string range $Line $f $l]
		xcDebug "ftn25(X$i): $ftn25(X$i)"
	    }
	} else {
	    if { $system(c95_version) == "98" || $system(c95_version) == "03" } {
		# read FORMAT: A3,I1,A4,2I4,(3E12.5)
		#              2E12.5
		#              6I3,1X,A36
		# -%-2DOSS   1 202 0.00000E+00 4.54430E+00-1.84013E-01		
		# 0.00000E+00-9.05903E+02
		# 1 36  0  0  0  0            **********
		set ftn25(NROW) [string range $Line 8  11]
		set ftn25(NCOL) [string range $Line 12 15]
		xcDebug "(NROW,NCOL): $ftn25(NROW) $ftn25(NCOL)"
		
		set ftn25(X3)   [string range $Line 16 27]
		set ftn25(X4)   [string range $Line 28 39]
		set ftn25(X5)   [string range $Line 40 51]
		
		#################
		gets $fileID Line
		#################
		
		set ftn25(X1)   [string range $Line  0 11]
		set ftn25(X2)   [string range $Line 12 23]
	    } elseif { $system(c95_version) == "06" || $system(c95_version) == "09" || $system(c95_version) == "14" } {
		# CRYSTAL06 read format: a3,i1,a4,2i5,1p,(3e12.5) ****
		#                        1p,3e12.5
		# -%-0DOSS    1  502 0.00000E+00 4.07550E-03-2.97766E-01
		# 1 36  0  0  0  0            **********
		set ftn25(NROW) [string range $Line 8  12]
		set ftn25(NCOL) [string range $Line 13 17]
		xcDebug "(NROW,NCOL): $ftn25(NROW) $ftn25(NCOL)"
		
		set ftn25(X3)   [string range $Line 18 29]
		set ftn25(X4)   [string range $Line 30 41]
		set ftn25(X5)   [string range $Line 42 54]
		
		#################
		gets $fileID Line
		#################
		
		set ftn25(X1)   [string range $Line  0 11]
		set ftn25(X2)   [string range $Line 12 23]

	    }

	    xcDebug "ftn25 X#1-5: $ftn25(X1) $ftn25(X2) $ftn25(X3) $ftn25(X4) $ftn25(X5)"
	    
	    #################
	    gets $fileID Line
	    #################
	    
	    for {set i 1} {$i <= 6} {incr i} {
		set f [expr ($i - 1) * 3]
		set l [expr $f + 2]
		set ftn25(I$i) [string range $Line $f $l]
		xcDebug "ftn25: I#$i= $ftn25(I$i)"
	    }
	}	    

	set N [expr $ftn25(NCOL) * $ftn25(NROW)]
	if { $what == "DOSS" } { set N $ftn25(NCOL) }
        
	set c 0
	if { $what != "BAND" } {
	    for {set i 1} {$i <= $N} {incr i} {
		# read FORMAT: 6(1PE12.5) several times
		if { $c == (6 * ($c / 6)) } {
		    set i6 0
		    #################
		    gets $fileID Line
		    #################
		}
		set f [expr $i6 * 12]
		set l [expr $f + 11]
		incr c
		incr i6

		set grafdata($gID,$c,$is,$ng) [string range $Line $f $l]
		xcDebug "value: $c $grafdata($gID,$c,$is,$ng)"

		if { $c == 1 } {
		    set ftn25(Ymin) $grafdata($gID,1,$is,$ng)
		    set ftn25(Ymax) $grafdata($gID,1,$is,$ng)
		} else {
		    if { $grafdata($gID,$c,$is,$ng) < $ftn25(Ymin) } { 
			set ftn25(Ymin) $grafdata($gID,$c,$is,$ng)
		    } elseif { $grafdata($gID,$c,$is,$ng) > $ftn25(Ymax) } {
			set ftn25(Ymax) $grafdata($gID,$c,$is,$ng)
		    }
		}
	    }
	} else {
	    # t.k: change for UHF
	    # read BANDS
	    set count 0
	    for {set i 1} {$i <= $ftn25(NCOL)} {incr i} {
		set j0 [expr 1 + $ftn25(NROW) * ($nproj - 1)]
		set jN [expr $ftn25(NROW) * $nproj]
		set ip [expr $grafdata($gID,N_point,$ng) + $i]
		xcDebug "(j0,jN,ip): $j0 , $jN , $ip"
		for {set j $j0} {$j <= $jN } {incr j} {
		    if { $count == (6 * ($count / 6)) } {
			set i6 0
			#################
			gets $fileID Line
			#################
		    }
		    set f [expr $i6 * 12]
		    set l [expr $f + 11]
		    set grafdata($gID,$ip,$j,$ng) [string range $Line $f $l]
		    incr i6
		    incr count

		    xcDebug "band: $j;  Kpoint: $ip; value: $grafdata($gID,$ip,$j,$ng)"		
		}
	    }
	}
    }
}


proc FTN25_MustBeNumber {value returnvalue text} {
    global system prop
    if [catch {expr abs($value)}] {
	set button [tk_dialog [WidgetName] ERROR $text error 0 OK Details]
	if $button {
	    DispC95Output $system(SCRDIR) $prop(file)25 \
		    "CRYSTAL Unit 25" 1
	}
	return $returnvalue
    }
    return $value
}
