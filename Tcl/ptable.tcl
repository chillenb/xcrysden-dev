#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/ptable.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc ptable {tplw {args {}}} {
    global ptable tcl_platform

    if { $tcl_platform(platform) == "windows" } {
	set bw 2
    } else {
	set bw 1
    }

    # tplw .... according to which toplevel the ptable's toplevel should be
    #           placed by xcToplevel command
    # args:
    #       -command command .... which command to execute when an 
    #                             element button is pressed; should set 
    #                             the ptable(result) variable

    set oldgrab [grab current]
    # if ptable(ptable_script) var exists we are running ptable as utility
    # so make no toplevel
    if { ![info exists ptable(ptable_script)] } {
	set t [xcToplevel [WidgetName] "XCrySDen: Periodic Table of Elements" \
		"Periodic Table" [winfo toplevel $tplw] 10 10]
    } else {
	set t .
	wm title . "XCrySDen: Periodic Table of Elements"
	wm iconname . "Periodic Table"
    }
    catch { grab $t }

    set i 0
    set command DummyProc
    set variable ptable(dummy_var)
    foreach option $args {
	incr i
	# odd cycles are tags, even options
        if { $i%2 } {
            set tag $option
        } else {
            switch -- $tag {
                "-command"  {set command $option}
		"-variable" {set variable $option}
		default { 
		    tk_dialog [WidgetName] Error \
			    "ERROR: Bad \"ptable\" configuration option $tag" \
			    error 0 OK 
		    return 0
		}

	    }
	}
    }
    if { $i%2 } {
	tk_dialog [WidgetName] \
		"ERROR: You called \"ptable\" with an odd number of args !" \
		error 0 OK
	return 0
    }

    set row(0)   {{} 1a  2a  3b  4b  5b  6b  7b  8_3  1b  2b  3a  4a  5a  6a  7a  8a}
    set row(1)   {I    H  x16 He}
    set row(2)   {II  Li  Be x10 B   C   N   O   F   Ne}
    set row(3)   {III Na  Mg x10 Al  Si  P   S   Cl  Ar}
    set row(4)   {IV  K   Ca  Sc  Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr}
    set row(5)   {V   Rb  Sr  Y   Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe}
    set row(6)   {VI  Cs  Ba  La  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn}
    set row(7)   {VII Fr  Ra  Ac  Rf  Ha  Sg  Ns  Hs  Mt x9}
    
    set row(9)   {Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu}
    set row(10)  {Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr}
    
    set col(0) "#ffffbb"
    set col(1) "#ff6688"
    set col(2) $col(1)
    foreach i {3 4 5 6 7 8 9 10 11 12} {
    	 set col($i) "#77bbff"
    }
    foreach i {13 14 15 16 17 18} {
    	 set col($i) "#ffcccc"
    }
    
    foreach el {3,13 4,13 5,13 6,13   4,14 5,14 6,14   5,15 6,15   6,16} {
    	 set col($el) $col(1)
    }
    
    . configure -bg "#ffffff"
    set f  [frame $t.f  -bg "#ffffff"]
    set f1 [frame $f.f1 -bg "#ffffff"]
    set f2 [frame $f.f2 -bg "#ffffff"]
    
    set font [font create]
    font configure $font \
    	     -family Times \
    	     -size   20 \
    	     -weight bold \
    	     -underline 1
    
    label $f.l -text "Periodic Table of Elements" -relief groove -bd 2 -font $font
    pack $f.l -side top -expand 1 -padx 10 -ipadx 5 -ipady 2
    focus $f.l
    #set f [frame .f]
    pack $f $f1 $f2 -ipadx 10 -ipady 10 -side top
    # take care of row(0)
    set ic 0
    foreach item $row(0) {
    	 set sticky {}
    	 set span 1
    	 if [string match *_* $item] {
    	     set c [string last _ $item]	
    	     set span [string range $item [expr $c + 1] end]
    	     set sticky ew
    	     set item  [string range $item 0 [expr $c - 1]]
    	     puts stdout "$c $span $item"
    	 }
    	 if { $ic != 0 } {
    	     puts stdout "$ic $span $item"
    	     set b [button $f1.0_$ic \
    		     -text $item \
    		     -relief ridge -bd 2 \
    		     -width $bw \
    		     -bg $col($ic) -disabledforeground "#000000" \
    		     -state disabled]
    	     grid $b
    	     grid configure $b -column $ic -row 0 \
    		     -columnspan $span -sticky $sticky -in $f1 -pady 5
    	     if { $span > 1 } {
    		 incr ic [expr $span - 1]
    	     }
    	 }
    	 incr ic
    }
    
    for {set ir 1} {$ir <= 7} {incr ir} { 
    	 set ic 0
    	 foreach item $row($ir) { 
    	     if { $ic == 0 } {
    		 set b [button $f1.${ir}_0 \
			 -command [list eval $command $t $item] \
    			 -text $item \
    			 -relief ridge -bd 2 \
    			 -width $bw \
    			 -bg $col($ic) -disabledforeground "#000000" \
    			 -state disabled]
    		 grid $b
    		 grid configure $b -column 0 -row $ir -in $f1 -padx 5
    		 incr ic
    		 continue
    	     }
    	     if { [string index $item 0] == "x" } {	    
    		 set n [string range $item 1 end]
    		 incr ic $n
    		 continue
    	     }
    	     
    	     set color $col($ic)
    	     if [info exists col($ir,$ic)] {
    		 set color $col($ir,$ic)
    	     }
    	     set b [button $f1.${ir}_${ic} \
		     -command [list eval $command $t $item] \
    		     -text $item \
    		     -width $bw \
    		     -bg $color]
    	     grid $b
    	     grid configure $b -column $ic -row $ir -in $f1
    	     incr ic
    	 }
    }
    	     
    set color "#99bbff"
    button $f2.9_0 -text "Lantanides" \
    	     -bg $color -state disabled \
    	     -relief ridge -bd 2 \
    	     -bg $color \
    	     -disabledforeground "#000000"
    
    button $f2.10_0 -text "Actinides" \
    	     -bg $color -state disabled \
    	     -relief ridge -bd 2 \
    	     -bg $color \
    	     -disabledforeground "#000000"
    grid  $f2.9_0 $f2.10_0
    grid configure $f2.9_0  -column 0 -row 9  -columnspan 3 -sticky we -in $f2
    grid configure $f2.10_0 -column 0 -row 10 -columnspan 3 -sticky we -in $f2
    
    set color "#99ddff"
    for {set ir 9} {$ir <= 10} {incr ir} {
	set ic 3
	foreach item $row($ir) { 
	    set b [button $f2.${ir}_$ic \
		    -command [list eval $command $t $item] \
		    -text $item \
		    -width $bw \
		    -bg $color]
	    grid $b
	    grid configure $b -column $ic -row $ir -in $f2
	    incr ic
	}
    }
    
    DefaultButton $f.b -text Close -command [list ptableClose $t]
    pack $f.b -side top -padx 10 -pady 10
    
    tkwait window $t
    if { $oldgrab != {} } {
	catch { grab $oldgrab }
    }

    upvar #0 $variable var
    set var $ptable(result)
    return $ptable(result)
}


proc ptableClose t {
    global ptable

    set ptable(result) {}
    CancelProc $t
}


proc ptableSelectElement {t el} {
    global ptable

    set ptable(result) $el
    CancelProc $t
}


proc ptableSelectAtomicNumber {t el} {
    global ptable

    set ptable(result) [Aname2Nat $el]
    CancelProc $t
}
    

proc scroll_ptableSelectNAT {t globarrayname arrel i} {
    upvar #0 $globarrayname array

    # globarrayname .... name of global array
    # arrel         .... name of array element
    # i             .... i-th element of array with element arrel
    
    set array($arrel,$i) [ptable $t -command ptableSelectAtomicNumber]
}

#lappend auto_path "/home/tone/prog/XCrys/Mesa"
#ptable . -command ptableSelectElement
