#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/bwid.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2014 by Anton Kokalj                                   #
#############################################################################

proc BWIDGraph {file} {
    global properties ftn25 grafdata grafsize system prop graf
    
    # FIX THIS
    if { ![file exists $file] } {
	tk_dialog .err ERROR "ERROR! \nFile \"$file\" does not exist !!! \
		This is bug in brogram, please report to author !!!\n\
		Code: BWIDGraph -file not found-" error 0 OK
	return
    }
    
    set fileID [open $file r]

    #
    # query GrapherID
    #
    set gID [NextGrapherID]
    xcDebug "GrapherID:: -- $gID --"

    set grafdata($gID,N_graf)      1
    set grafdata($gID,N_segment,1) 2; # one is min band E, one is max band E
    set grafdata($gID,N_text)      1
    set grafdata($gID,barshadow,1) 0
    set graf($gID,RorigX,1)        0.0
    set graf($gID,RorigY,1)        0.0
    set graf($gID,RsizeX,1)        1.0
    set graf($gID,RsizeY,1)        1.0

    set done 0
    while { !$done } {
	gets $fileID line
	if { $system(c95_version) == "03" || $system(c95_version) == "06" \
		 || $system(c95_version) == "09" || $system(c95_version) == "14" } {
	    #
	    # CRYSTAL03
	    #
	    # grep line: N. OF SCF CYCLES                8  FERMI ENERGY              -0.511E-01
	    set pattern1 "*N. OF SCF CYCLES*FERMI ENERGY*"
	    set pattern2 "*EFERMI(AU)*"
	} else {
	    set pattern1 "*TOTAL *FERMI ENERGY*"
	    set pattern2 "*EFERMI (A.U.)*"
	}

	if { [string match $pattern1 $line] } {
	    putsFlush stderr "EF line: $line"
	    set grafdata($gID,N_Yline,1)        1
	    set grafdata($gID,Yline1,1)         [lindex $line end]
	    set grafdata($gID,Yline1_text,1)    "E_Fermi"
	    set grafsize($gID,Yline1_fill,1)    #ff0000
	    set grafsize($gID,Yline1_width,1)   1
	    set grafsize($gID,Yline1_stipple,1) @$system(BMPDIR)/dotH.bmp
	}
	# may be Fremi energy was recalculated
	if { [string match $pattern2 $line] } {
	    puts stdout "EF line: $line"
	    set grafdata($gID,Yline1,1)         [lindex $line 5]
	}
	if { [regexp "BAND LIMITS" $line] } {
	    xcDebug "$line"
	    set band1 [lindex $line 4]
	    set band2 [lindex $line 7]
	    # read blank line
	    gets $fileID line
	    set done 1
	    break; # leave this break - otherwise it doesn't work (I don't know why)
	}
	xcDebug $line
    }

    set grafdata($gID,X_title,1)     "BAND number"
    set grafdata($gID,Y_title,1)     "E \[a.u.\]"
    set grafdata($gID,N_point,1)     [expr $band2 - $band1 + 1]
    set grafdata($gID,N_Xline,1)     [expr $grafdata($gID,N_point,1) - 1]
    set grafdata($gID,firstbar,1)    $band1
    set grafdata($gID,lastbar,1)     $band2
    # arbitrarily we set: Xmin to 0 & Xmax to 1
    set grafdata($gID,Xmin,1)        0.0
    set grafdata($gID,Xmax,1)        1.0
    set grafdata($gID,Xoffset,1)     0.2
    set grafdata($gID,Yoffset,1)     0.05
    set grafdata($gID,N_MXtick,1)    $grafdata($gID,N_point,1)
    set grafdata($gID,N_mXtick,1)    0
    set grafsize($gID,Mtick_size,1)  0
    set grafsize($gID,mtick_size,1)  0

    if { $grafdata($gID,N_point,1) < 1 } {
	tk_dialog .dialog WARNING "WARNING: Number of BANDs to plot is 0!!!" \
		error 0 OK
	return
    }
    
    set ng 1
    if { $prop(type_of_run) == "UHF" } {
	set ng 2
    }
    set grafdata($gID,barn,1) $ng

    set l [expr 1.0 / $grafdata($gID,N_point,1)]
    puts stdout "Npoint: $grafdata($gID,N_point,1);      lll: $l"
    
    for {set ig 1} {$ig <= $ng} {incr ig} {
	if { $ig == 2 } {
	    set done 0
	    while { !$done } {
		gets $fileID line
		if { [regexp "BAND LIMITS" $line] } {
		    set done 1
		    break
		}
	    }
	    # read blank line
	    gets $fileID line
	}
	for {set i 1} {$i <= $grafdata($gID,N_point,1)} {incr i} {
	    gets $fileID line
	    # t.k.:
	    while { [llength $line] == 0 } {
		gets $fileID line		
	    }
	    #/
	    if { $ig == 1 } {
		set grafdata($gID,Xtick${i}_text,1) [lindex $line 1]
		set grafdata($gID,Xtick$i,1) [expr ($i - 1) * $l + $l / 2]
		set grafdata($gID,X$i,1)     [expr ($i - 1) * $l + $l / 2]
		set grafsize($gID,bar${i}_fill,1)     #ff0000
		set grafsize($gID,bar${i}_outline,1)  #ff0000
		set grafsize($gID,bar${i}_shadow,1)   #aaaaaa
	
		if { $i < $grafdata($gID,N_point,1) } {	    
		    set grafdata($gID,Xline$i,1)          [expr $i * $l]
		    set grafsize($gID,Xline${i}_fill,1)   #000
		    set grafsize($gID,Xline${i}_width,1)  1
		    set grafsize($gID,Xline${i}_stipple,1) \
			    @$system(BMPDIR)/dotV.bmp
		}
	    } else {
		set grafsize($gID,bar${i}_fill,2)     #0000ff
		set grafsize($gID,bar${i}_outline,2)  #0000ff
		set grafsize($gID,bar${i}_shadow,2)   #999999
	    }
	    
	    # Y values
	    if { [llength $line] >= 5 } {
		set grafdata($gID,$i,1,$ig)           [lindex $line 3]
		set grafdata($gID,$i,2,$ig)           [lindex $line 5]
		puts stdout "VALUES: $grafdata($gID,$i,1,$ig)    $grafdata($gID,$i,2,$ig)"
	    } else {
		set grafdata($gID,$i,1,$ig) +999
		set grafdata($gID,$i,2,$ig) -999
	    }

	    if { $i == 1 && $ig == 1 } {
		set Ymin $grafdata($gID,$i,1,$ig)
		set Ymax $grafdata($gID,$i,2,$ig)
	    } else {
		if { $grafdata($gID,$i,1,$ig) < $Ymin } {
		    set Ymin $grafdata($gID,$i,1,$ig)
		}
		if { $grafdata($gID,$i,2,$ig) > $Ymax } {
		    set Ymax $grafdata($gID,$i,2,$ig)
		}
	    }
	    	    
	    set grafsize($gID,bar${i}_stipple,$ig)  gray50
	    set grafsize($gID,bar${i}_width,$ig)    1
	
	}
    }

    set yrange                       [expr $Ymax - $Ymin]
    set grafdata($gID,Ymin,1) \
	    [expr $Ymin - $yrange * $grafdata($gID,Yoffset,1)]
    set grafdata($gID,Ymax,1) \
	    [expr $Ymax + $yrange * $grafdata($gID,Yoffset,1)]
		
    puts stdout "YMin YMax:: $grafdata($gID,Ymin,1) $grafdata($gID,Ymax,1)"
	 
    set grafdata($gID,N_MYtick,1)    6
    set grafdata($gID,N_mYtick,1)    2
    
    set dy [expr ($grafdata($gID,Ymax,1) - $grafdata($gID,Ymin,1)) / \
	    ($grafdata($gID,N_MYtick,1) - 1)]
    for {set i 1} {$i <= $grafdata($gID,N_MYtick,1)} {incr i} {
	set grafdata($gID,Ytick$i,1)        [expr $Ymin + ($i - 1) * $dy]
	set grafdata($gID,Ytick${i}_text,1) \
		[TickFormat $grafdata($gID,Ytick$i,1)]
	puts stdout "YTICK:: $grafdata($gID,Ytick$i,1)"
    }
    set grafsize($gID,canW) 600
    set grafsize($gID,canH) 450	

    close $fileID
}




