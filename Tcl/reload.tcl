#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/reload.tcl
# ------                                                                    #
# Copyright (c) 1996--2019 by Anton Kokalj                                  #
#############################################################################


#------------------------------------------------------------------------
# This proc is called before a reloadable file is opened, so that we
# register what to do when reoad button is pressed
#------------------------------------------------------------------------
proc reloadRegister {filedir type} {
    global reload

    set reload(capable) 1
    set reload(file) $filedir

    set allowed_types {
        XSF XSF2 XYZ GZMAT PDB addOption:reload openExtStruct:reload WIEN2k
    }
    if { [lsearch $allowed_types $type] > -1 } {        
        set reload(type) $type
    } else {
        error "unsupported reload type $type; must be one of [join $allowed_types {, }]"
    }
}
#------------------------------------------------------------------------
# This proc is the deafult command for reloading the structure
#------------------------------------------------------------------------
proc reloadUpdateDefault {reloadCmd filedir {can .mesa}} {

    if { [info exists select(selection_mode)] && $select(selection_mode) } {
	tk_dialog [WidgetName] "WARNING !!!" "WARNING: structure can not be reloaded while in selection mode. First exit from the selection mode, then reload" warning 0 OK
	return
    }

    $reloadCmd $filedir $can 1
}

#------------------------------------------------------------------------
# This proc return 1 if we are reload-capable and 0 otherwise
#------------------------------------------------------------------------
proc reloadCapable {} {
    global reload
    
    if { [info exists reload(capable)] && $reload(capable) } {
        return 1
    }
    return 0
}

#------------------------------------------------------------------------
# This proc register's reload button 
#------------------------------------------------------------------------
proc reloadButton {buttonWidget} {
    global reload

    set reload(buttonWidget) $buttonWidget
    reloadButtonQueryState
}

#------------------------------------------------------------------------
# This proc enables reload button
#------------------------------------------------------------------------
proc reloadButtonEnable  {} { reloadButtonState_ normal }

#------------------------------------------------------------------------
# This proc disables reload button
#------------------------------------------------------------------------
proc reloadButtonDisable {} { reloadButtonState_ disabled }

#------------------------------------------------------------------------
# This proc enables (disables) reload button if reload-capable (reload-incapable)
#------------------------------------------------------------------------
proc reloadButtonQueryState {} {
    expr { [reloadCapable] ? [reloadButtonEnable] : [reloadButtonDisable] }
}

# internal proc ... (do not call directly)
proc reloadButtonState_ {state} {
    global reload    
    if { ! [info exists reload(buttonWidget)] } {
        return
    }
    if { [winfo exists $reload(buttonWidget)] } {
        $reload(buttonWidget) configure -state $state
    }
}

#------------------------------------------------------------------------
# This proc is a reload-command; it is executed when reload button is pressed
#------------------------------------------------------------------------
proc reloadCmd {{togl .mesa}} {
    global reload

    if { ! [info exists reload(type)] } {
        error "reloadRegister or reloadRegisterCmd must be called before reloadCmd"
    }

    switch -- $reload(type) {
        XSF {
            xsfUpdate $reload(file) $togl         
        }
        XSF2 {
            reloadUpdateDefault xsf2Open $reload(file) $togl
        }
        XYZ {
            reloadUpdateDefault xyzOpen $reload(file) $togl
        }
        GZMAT {
            reloadUpdateDefault gzmat $reload(file) $togl
        }
        addOption:reload -
        openExtStruct:reload {
            reloadUpdateDefault $reload(type) $reload(file) $togl
        }
        PDB {
            pdfUpdate $reload(file) $togl
        }
        WIEN2k {
            wnOpenSFileUpdate $reload(file) $togl
        }
        default {
            error "unsupported reload type $type"
        }
    }
}

