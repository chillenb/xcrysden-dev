#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/number.tcl                                       
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc check_var { varlist foclist } {
    global err ok
    # $varlist - list of elements, which are lists themselfs
    #            element = {varname vartype}
    #            ^^^^^^^^^^^^^^^^^^^^^^^^^^^
    #                      if $vartype is not specefied, than default is REAL
    #                      (look the: "proc number { varname {type real} }")
    #
    # $foclist - list of widgets to focus (if an error occur)
    
    set n 0
    foreach elem $varlist {
	set err 0
	set varname [lindex $elem 0]
	set vartype [lindex $elem 1]
	if { $vartype == "text" } {
	    DummyProc
	} else {
	    number $varname $vartype
	    if { $err } { 
		focus [lindex $foclist $n]
		break 
	    }
	    incr n
	}
    }
    if { $err == 0 } {
	set ok 1
	return 1
    } else {
	return 0
    }    
}


proc number { varname {type real} } {
    global err 
    upvar #0 $varname var 

    # $varname - name of variable
    # $var - value of variable
    # $err - boolean pointer for error determining
    # $type - type of number (default is REAL)
    # type of number: int - integer
    #             intlist - list of integers
    #              posint - positive integer
    #               fract - fractional number [-1,1]
    #                real - real number; (all numbers are real 
    #                                       (excluded imaginary))
    #             posreal - positive real number;
    #                 nat - atomic number (0-300; lokk in C95 manual
    #                                      0 is ghost atom)
    #        intrange a b - integer interval [a,b]
    set err 0
    #puts stdout "VAR:: \"$var\""
    if { ! [info exists var] } { return }
    # may be $var is not specified at all
    if { $type == "intlist" } {
	foreach int $var {
	    if { [catch {expr abs($var)}] || \
		    $var != int($var) || [string match *.* $var] } {
		dialog .number2 ERROR "ERROR !\nYou have specified \
			a non-integer number instead of integer number \
			for \"$varname\" variable.\n\
			TRY AGAIN !" error 0 OK
		set err 1 
		return [expr 1 - $err]
	    }
	}
    }

    if { $var == "" } {
	dialog .number1 ERROR "ERROR !\nYou forget to specify \
		\"$varname\" variable. PLEASE DO IT \!" error 0 OK
	set err 1
    } elseif { [catch {expr abs($var)}] } {
	# this CATCH specify if $var is a number;
	# if we get 1 --> not number, else number
	# string is not a number
	dialog .number1 ERROR "ERROR !\nYou have specified a character \
		instead of number for \"$varname\" variable.\
		TRY AGAIN \!" error 0 OK
	set err 1
    } else {
	# string is a number
	# varify if number is a right one
	switch -glob -- $type {
	    int { 
		if { $var != int($var) || [string match *.* $var] } {
		    dialog .number2 ERROR "ERROR !\nYou have specified \
			    a non-integer number instead of integer number \
			    for \"$varname\" variable.\n\
			    TRY AGAIN !" error 0 OK
		    set err 1 
		}
	    }
	    posint { 
		if { $var != int($var) || $var < 0 || \
			[string match *.* $var] } {
		    dialog .number2 ERROR "ERROR !\nYou have specified \
			    a non-positive integer instead of positive \
			    integer for \"$varname\" variable.\n\
			    TRY AGAIN !" error 0 OK
		    set err 1 
		}
	    }
	    fra* { 
		if { $var < -1 || $var > 1 } {
		    dialog .number2 ERROR "ERROR !\nYou should specify \
			    a number between \[-1,1\] for \"$varname\" \
			    variable. TRY AGAIN !" error 0 OK
		    set err 1 
		}
	    }
	    nat* { 
		if { $var != int($var) || $var < 0 || $var > 399 || \
			[string match *.* $var] } {
		    dialog .number2 ERROR "ERROR !\nYou specify wrong Atomic \
			    Number; Atomic Number should be between [0-99],\
			    [100-199] and [200-299], [300-399].\n\
			    TRY AGAIN !" error 0 OK
		    set err 1 
		}
	    }
	    posreal { 
		if { $var < 0.0 } {
		    dialog .number2 ERROR \
			    "ERROR !\nYou specify a negative real number \
			    instead of positive one.\n TRY AGAIN !" error 0 OK
		    set err 1 
		}
	    }
	    intrange* { 
		set a [lindex $type 1]
		set b [lindex $type 2]
		if { $var != int($var) || $var < $a || $var > $b } {
		    dialog .number2 ERROR \
			    "ERROR !\nYou specify a number that is either \
			    non-integer or out of range.\n TRY AGAIN !" \
			    error 0 OK
		    set err 1 
		} 
	    }		
	}
    }
    return [expr 1 - $err]
}


##############################################################################
# similar as previous, but it return 1 if "varname" is the right type
# else it return 0
proc xcNumber { varname {type real} } {
    global err 
    upvar #0 $varname var 

    # $varname - name of variable
    # $var - value of variable
    # $err - boolean pointer for error determining
    # $type - type of number (default is REAL)
    # type of number: int - integer
    #              posint - positive integer
    #               fract - fractional number [-1,1]
    #                real - real number; (all numbers are real 
    #                                       (excluded imaginary))
    #                 nat - atomic number (0-98; 98 is a program limitaition;
    #                                      0 is ghost atom)

    set err 0
    puts stdout "VAR:: \"$var\""
    # may be $var is not specified at all
    if { $var == "" } {
	return 0
    } elseif { [catch {expr abs($var)}] } {
	return 0
    } else {
	# string is a number
	# varify if number is a right one
	switch -glob -- $type {
	    int* { 
		if { $var != int($var) || [string match *.* $var] } { 
		    return 0	    
		}
	    }
	    pos* { 
		if { $var != int($var) || $var < 0 || \
			[string match *.* $var] } { return 0 }
	    }
	    fra* { 
		if { $var < -1 || $var > 1 } { return 0 }
	    }
	    nat* { 
		if { $var != int($var) || $var < 0 || $var > 399 || \
			[string match *.* $var] } { return 0 }		
	    }
	}
    }
    return 1
}




