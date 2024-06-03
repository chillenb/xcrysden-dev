#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/xctr.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc xctr_proc {name args body} {
    global xctr
   
    #set xctr(procName) $name
    #set xctr(argsVars) $args

    set _body [subst -nocommands {	
	global env xctr

	set __args__ ""
	if { [list $args] != [list ""] } {
	    foreach __arg__ [list $args] {
		set __var__ [lindex \$__arg__ 0]
		upvar 0 \$__var__ __value__
		#puts stderr "XCTRvar> \$__var__"
		lappend __args__ \$__value__	    
	    }
	}
	
	if { [info exists env(XCRYSDEN_TRACE)] && (1 >= [info level]) } {
	    if { \$xctr(recording) } {
		puts  $xctr(traceID) [concat $name \$__args__]
		flush $xctr(traceID)
		puts  stderr [concat XCTRexec> $name \$__args__ , level == [info level]]
	    }
	} elseif { [info exists env(XCRYSDEN_FULLTRACE)] } {
	    set level "; level == [info level]"
	    puts  $xctr(fulltraceID) [concat $name \$__args__ \$level]	    	   
	    #puts  $xctr(fulltraceID) "[format %[expr [info level]*3]s { }]$name $args"
	    flush $xctr(fulltraceID)
	} 
	$body
    }]

    tcl_proc $name $args $_body
}

if { [info exists env(XCRYSDEN_TRACE)] || [info exists env(XCRYSDEN_FULLTRACE)] } {
    #set xctr(traceID)     [open [file join $env(XCRYSDEN_SCRATCH) xc_trace.[pid]]      w]
    #set xctr(fulltraceID) [open [file join $env(XCRYSDEN_SCRATCH) xc_fulltrace.[pid]] w]
    set xctr(traceID)     stderr
    set xctr(fulltraceID) stderr
    rename proc      tcl_proc
    rename xctr_proc proc
}

