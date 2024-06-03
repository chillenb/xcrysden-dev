#############################################################################
# Author:                                                                   #
# ------                                                                    #
#  Anton Kokalj                                  Email: Tone.Kokalj@ijs.si  #
#  Department of Physical and Organic Chemistry  Phone: x 386 1 477 3523    #
#  Jozef Stefan Institute                          Fax: x 386 1 477 3811    #
#  Jamova 39, SI-1000 Ljubljana                                             #
#  SLOVENIA                                                                 #
#                                                                           #
# Source: $XCRYSDEN_TOPDIR/Tcl/mpegParam.tcl
# ------                                                                    #
# Copyright (c) 1996-2003 by Anton Kokalj                                   #
#############################################################################

proc mpegCreateParamFile {output_file input_dir input_files} {
    global system gifAnim myParam

    if { $gifAnim(frame_files_format) == "PPM" } {
	set input_convert "INPUT_CONVERT    cat *"
    } else {
	set input_convert "INPUT_CONVERT    jpegtopnm *"
    }

    if { [info exists myParam(MPEG_ENCODE_PARAM_FILE)] } {
	# encode_param file parameters are specified in the definition file
	
	set encode_param [subst $myParam(MPEG_ENCODE_PARAM_FILE)]

    } else {
	# load the dafult definition

	set encode_param [subst {
# ------------------------------------------------------------------------
# Please edit this ppmtompeg (MPEG_ENCODE) parameter file to suit your needs
# ------------------------------------------------------------------------

PATTERN          IBBPBBPBBPBBPBBP
OUTPUT           $output_file
BASE_FILE_FORMAT PPM

# ------------------------------------------------------------------------
# Put here the appropriate image conversion/nahdling program.
# Native format is PPM, therefore if your frame-files are in 
# PPM format put:
#
# INPUT_CONVERT    cat *
#
# if your files are in JPEG format, then the following would do:
#
# INPUT_CONVERT    jpegtopnm *
# ------------------------------------------------------------------------
$input_convert
GOP_SIZE         16
SLICES_PER_FRAME 1

#INPUT_DIR        $input_dir
INPUT_DIR      /
INPUT
[join $input_files \n]
END_INPUT

PIXEL           FULL
RANGE           10
PSEARCH_ALG     EXHAUSTIVE
BSEARCH_ALG     CROSS2
IQSCALE         7
PQSCALE         10
BQSCALE         15
REFERENCE_FRAME DECODED
BUFFER_SIZE     327680
FRAME_RATE      23.976
BIT_RATE        10000000

# There are many more options, see the users manual for examples....
# ASPECT_RATIO, USER_DATA, GAMMA, IQTABLE, etc.
	}]
    }
	 
    if { $gifAnim(temp_files_dir) == "pwd" } {
	set dir $system(PWD)
    } else {
	set dir $system(SCRDIR)
    }	    

    set file [file join $dir mpeg_encode.param]
    WriteFile $file $encode_param w
    
    if { $gifAnim(edit_param) } {
	xcEditFile $file foreground
    }
    return $file
}	
