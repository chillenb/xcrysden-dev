#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/bar.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

#----------------------
# BARGraph is called as:
#----------------------
#	grapher_BARGraph $file {$args {}}
#	Grapher BARGraph

#############################
# format of bargarph file is:
#----------------------------
# nbar
# 1 Ymin1 Ymin2
# ....
#############################

proc grapher_ReadBARFile file {
    global grafdata grafsize gID system

    set nl 0
    foreach lin [split [ReadFile $file] \n] {
	incr nl
	set line($nl) $lin
    }
    
    set nbar [lindex $line(1) 0]

    set grafdata($gID,N_point,1) $nbar
    set l [expr 1.0 / $grafdata($gID,N_point,1)]

    set i 0
    for {set ii 2} {$ii <= [expr $nbar + 1]} {incr ii} {
	incr i
	# X values
	set grafdata($gID,Xtick${i}_text,1) [lindex $line($ii) 0]
	set grafdata($gID,Xtick$i,1) [expr ($i - 1) * $l + $l / 2]
	set grafdata($gID,X$i,1)     [expr ($i - 1) * $l + $l / 2]
	set grafsize($gID,bar${i}_fill,1)     #ff0000
	set grafsize($gID,bar${i}_outline,1)  #ff0000
	set grafsize($gID,bar${i}_shadow,1)   #aaaaaa

	# Y values	
	set grafdata($gID,$i,1,1) [lindex $line($ii) 1]
	set grafdata($gID,$i,2,1) [lindex $line($ii) 2]
	xcDebug -debug "$grafdata($gID,$i,1,1) $grafdata($gID,$i,2,1)"
	if { $i < $grafdata($gID,N_point,1) } {	    
	    set grafdata($gID,Xline$i,1)          [expr $i * $l]
	    set grafsize($gID,Xline${i}_fill,1)   #000
	    set grafsize($gID,Xline${i}_width,1)  1
	    set grafsize($gID,Xline${i}_stipple,1) \
		    @$system(BMPDIR)/dotV.bmp
	    set grafsize($gID,Xline1_textstipple,1) @$system(BMPDIR)/dotH.bmp
	}

	if { $i == 1 } {
	    set Ymin $grafdata($gID,$i,1,1)
	    set Ymax $grafdata($gID,$i,2,1)
	} else {
	    if { $grafdata($gID,$i,1,1) < $Ymin } {
		set Ymin $grafdata($gID,$i,1,1)
		xcDebug -debug "Ymin1: $i,$Ymin,$grafdata($gID,$i,1,1)"
	    }
	    if { $grafdata($gID,$i,2,1) > $Ymax } {
		set Ymax $grafdata($gID,$i,2,1)
	    }
	}
	xcDebug -debug "Ymin: $Ymin"
	set grafsize($gID,bar${i}_stipple,1)  gray50
	set grafsize($gID,bar${i}_width,1)    1
    }    

    # arbitrarily we set: Xmin to 0 & Xmax to 1
    set grafdata($gID,Xmin,1)        0.0
    set grafdata($gID,Xmax,1)        1.0
    set grafdata($gID,Xoffset,1)     0.2
    set grafdata($gID,Yoffset,1)     0.05
    set grafdata($gID,N_MXtick,1)    $grafdata($gID,N_point,1)
    set grafdata($gID,N_mXtick,1)    0
    set grafdata($gID,barn,1)        1

    set grafsize($gID,Mtick_size,1)  0
    set grafsize($gID,mtick_size,1)  0

    set yrange [expr $Ymax - $Ymin]
    set grafdata($gID,Ymin,1) \
	    [expr $Ymin - $yrange * $grafdata($gID,Yoffset,1)]
    set grafdata($gID,Ymax,1) \
	    [expr $Ymax + $yrange * $grafdata($gID,Yoffset,1)]
			 
    set grafdata($gID,N_MYtick,1)    6
    set grafdata($gID,N_mYtick,1)    2
    
    set dy [expr ($grafdata($gID,Ymax,1) - $grafdata($gID,Ymin,1)) / \
	    ($grafdata($gID,N_MYtick,1) - 1)]
    for {set i 1} {$i <= $grafdata($gID,N_MYtick,1)} {incr i} {
	set grafdata($gID,Ytick$i,1)        [expr $Ymin + ($i - 1) * $dy]
	set grafdata($gID,Ytick${i}_text,1) \
		[TickFormat $grafdata($gID,Ytick$i,1)]
    }
    set grafsize($gID,canW) 600
    set grafsize($gID,canH) 450	
}

proc grapher_BARGraph {file {args {}}} {
    global grafdata grafsize graf system gID

    # options:
    #  "-Yline" 
    #  "-Yline_text" 
    #  "-Xline" 
    #  "-Xline_text" 
    #  "-Xtitle" 
    #  "-Ytitle" 
    #  "-firstbar" 
    #  "-lastbar"
    if ![file exists $file] {
	tk_dialog .err ERROR "ERROR! \nFile \"$file\" does not exist !!!" \
		error 0 OK
	return
    }
    set gID [NextGrapherID]
    grapher_ReadBARFile $file

    set grafdata($gID,N_Yline,1)        0
    set grafdata($gID,Yline1,1)         {}
    set grafdata($gID,Yline1_text,1)    {}
    set grafsize($gID,Yline1_fill,1)    #ff0000
    set grafsize($gID,Yline1_width,1)   1
    set grafsize($gID,Yline1_stipple,1) @$system(BMPDIR)/dotH.bmp

    set grafdata($gID,N_Xline,1)        0
    set grafdata($gID,Xline1_text,1)    {}
    set grafsize($gID,Xline1_fill,1)    #ff0000
    set grafsize($gID,Xline1_width,1)   1
    set grafsize($gID,Xline1_stipple,1) @$system(BMPDIR)/dotV.bmp

    set grafdata($gID,X_title,1)     "X title"
    set grafdata($gID,Y_title,1)     "Y title"
    set grafdata($gID,N_Xline,1)     [expr $grafdata($gID,N_point,1) - 1]
    set grafdata($gID,firstbar,1)    1
    set grafdata($gID,lastbar,1)     $grafdata($gID,N_point,1)

    xcDebug "grapher_BARGraph:: $args"
    set i 0
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
	    xcDebug "FillEntries Options:: $tag $option"
            switch -- $tag {
		"-Yline" {
		    set grafdata($gID,N_Yline,1)        1
		    set grafdata($gID,Yline1,1)         $option
		}
		"-Yline_text" {
		    set grafdata($gID,Yline1_text,1)    $option
		}
		"-Xline" {
		    set grafdata($gID,N_Xline,1)        1
		    set grafdata($gID,Xline1,1)         $option
		}
		"-Xline_text" {
		    set grafdata($gID,Xline1_text,1)    $option
		}
		"-Xtitle" {
		    set grafdata($gID,X_title,1)    $option
		}
		"-Ytitle" {
		    set grafdata($gID,Y_title,1)    $option
		}
		"-firstbar" {
		    set grafdata($gID,firstbar,1)   $option
		}
		"-lastbar" {
		    set grafdata($gID,lastbar,1)    $option
		}
		default { 
		    tk_dialog .mb_error Error \
			    "ERROR: Bad graph_BARGraph configure option $tag" \
			    error 0 OK 
		    return 0
		}
	    }
	}
    }
    if { $i%2 } {
	tk_dialog .mb_error1 Error \
		"ERROR: You called graph_BARGraph with an odd number of args !" \
		error 0 OK
	return 0    
    }
    
    set fileID [open $file r]

    set grafdata($gID,N_graf)      1
    set grafdata($gID,N_segment,1) 2; # one is min band E, one is max band E
    set grafdata($gID,N_text)      1
    set grafdata($gID,barshadow,1) 0
    set graf($gID,RorigX,1)        0.0
    set graf($gID,RorigY,1)        0.0
    set graf($gID,RsizeX,1)        1.0
    set graf($gID,RsizeY,1)        1.0

    if { $grafdata($gID,N_point,1) < 1 } {
	tk_dialog .dialog WARNING "WARNING: Number of BANDs to plot is 0!!!" \
		error 0 OK
	return
    }
}
