#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/cygwin.tcl
# ------                                                                    #
# Copyright (c) 2004 by Anton Kokalj                                        #
#############################################################################

#
# this file contains a few very dirty hacks needed to run under CYGWIN
#

if { $xcrys(platform) == "windows" } {
    # testing ...
    rename exec _tcl_exec

    proc exec {args} {
	global env

    	# first try a normal exec
	set result {}
    	if { [catch {set result [uplevel 1 [list eval _tcl_exec $args]]}] } {
    	    # try to execute via "launch.sh" wrapper
    	    set result [uplevel 1 [list eval _tcl_exec sh $env(XCRYSDEN_TOPDIR)/scripts/launch.sh] $args]
    	}
	return $result	
    }    
}

if { $xcrys(platform) == "cygwin" }  {
    rename pwd  _tcl_pwd    
    rename exec _tcl_exec
    rename tk_getSaveFile _tk_getSaveFile
    rename tk_getOpenFile _tk_getOpenFile

    proc exec {args} {
	global env
	
    	# first try a normal exec
	set result {}
    	if { [catch {set result [uplevel 1 [list eval _tcl_exec $args]]}] } {
    	    # try to execute via "launch.sh" wrapper
    	    set result [uplevel 1 [list eval _tcl_exec sh $env(XCRYSDEN_TOPDIR)/scripts/launch.sh] $args]
    	}
	return $result	
    }   
 
    proc pwd {} {
	global env

	# Don't use Tcl pwd command, it will return C:\cygwin\...,
	# and CYGWIN can't handle it. Use a pwd.sh wrapper instead.

	return [exec $env(XCRYSDEN_TOPDIR)/scripts/pwd.sh]
    }    

    proc tk_getSaveFile {args} {
	set sfile [eval _tk_getSaveFile $args]
	return [cygwin_unixpath $sfile]
    }

    proc tk_getOpenFile {args} {
	set sfile [eval _tk_getOpenFile $args]
	return [cygwin_unixpath $sfile]
    }

    proc cygwin_unixpath {path} {
	if { [regexp -- {^[A-Z]:/} $path] } {
	    set drive [string tolower [string index $path 0]]
	    regsub {^[A-Z]:} $path /cygdrive/$drive path
	}
	return $path
    }
}
