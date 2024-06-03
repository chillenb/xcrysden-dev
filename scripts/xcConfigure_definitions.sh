# ------------------------------------------------------------------------
#   FUNCTIONS
# ------------------------------------------------------------------------

type_p() {
    which $1
     #| awk '{print $NF}'
}
CustomDir () {
    $ECHO_n "
   I am about to create the \$HOME/.xcrysden/custom-definitions file,
   which is XCRYSDEN's customization file.

   Press <Enter> to continue ... "

   read ans
}

ExeWarning () {
    if [ ! -x "$1" ]; then
	echo "
   WARNING: file $1 is not executable !!!"
    fi
}

BackupName() {
    # usage: $0
    # return backup name such that the file with this name does not yet exists

    head=$1
    file=$1
    i=1
    while test -f $file
    do
	file=$head.$i
	i=`expr $i + 1`
    done
    echo $file
}
	
QueryProgram() {
    # usage: QueryProgram msg1 mesg2 deafult variable
    $ECHO_n "$1 ([y]es/[n]o): "
    read ans; 
    if [ \( "$ans" != "y" \) -a \( "$ans" != "n" \) ]; then
	echo "
   Please answer y or n !!!
"
	QueryProgram "$1" "$2" "$3" "$4"
	return
    elif [ "$ans" = "y" ]; then
	$ECHO_n "$2 (default: $3): "
	read prog;
	if [ "$prog" = "" ]; then
	    prog=$3
	fi
	ExeWarning "$prog"
	eval $4=\"$prog\"
    else
	unset $4
	#eval $4=""
    fi
}

GetYesNoAnswer() {
    # Usage: $0 text variable
    $ECHO_n "$1 ([y]es/[n]o): "
    read ans_

    if test \( "$ans_" != "y" \) -a \( "$ans_" != "n" \) ; then
	echo "
   Please answer y or n !!!
"
	GetYesNoAnswer "$1" "$2"
	return
    else
	eval $2=\"$ans_\"
    fi
}

CustomDef () {
    #
    # Query the CRYSTAL program
    #
    echo "
   CRYSTAL is an electronic structure program for periodic systems.
   (http://www.crystal.unito.it/)
"
    GetYesNoAnswer "   Do you have a CRYSTAL package" cXX

    echo "" > $tempdir/customFile.$$

    if [ "$cXX" = "y" ]; then	
        cry_def=`type_p crystal`
        if test "x$cry_def" = "x"; then
	    cry_def=/usr/local/bin/crystal
	fi
	$ECHO_n "
   Specify CRYSTAL's crystal module (default: $cry_def):  "
	read cry
	if [ "$cry" = "" ]; then
	    cry=$cry_def
	fi
	ExeWarning "$cry"

	pro_def=`type_p properties`
        if test "x$pro_def" = "x"; then
	    pro_def=/usr/local/bin/properties
	fi
	$ECHO_n "
   Specify CRYSTAL's properties module (default: $pro_def):  "
	read pro
	if [ "$pro" = "" ]; then
	    pro=$pro_def
	fi
	ExeWarning "$pro"
	echo "
# ------------------------------------------------------------------------
#  do we have CRYSTAL package
# ------------------------------------------------------------------------

set system(c95_exist) 1


# ------------------------------------------------------------------------
#  CRYSTAL modules
# ------------------------------------------------------------------------

set system(c95_crystal)    $cry
set system(c95_properties) $pro

" >> $tempdir/customFile.$$
    else
	echo "
# ------------------------------------------------------------------------
#  do we have CRYSTAL package
# ---------------------------

set system(c95_exist) 0


# ------------------------------------------------------------------------
#  CRYSTAL modules
# ------------------------------------------------------------------------

#set system(c95_crystal)    /full/path/to/crystal
#set system(c95_properties) /full/path/to/properties

" >> $tempdir/customFile.$$
    fi

    echo "

# ------------------------------------------------------------------------
# xcrysden can use several encoder programs for creating animated GIF
# (convert, gifsicle, whirlgif) and AVI/MPEG movies
# (mencoder/ppmtompeg)
# ------------------------------------------------------------------------

# what program we use for animated-gif encoding
#set xcMisc(gif_encoder)   convert

# what program we use for AVI/MPEG encoding
#set xcMisc(movie_encoder) mencoder


# ------------------------------------------------------------------------
# NOTICE: Starting from version 1.6, xcrysden tries to automatically
# find various external packages, nevertheless user can still set
# them explicitly, as shown below (if you want to do so, uncomment
# corresponding lines)
# ------------------------------------------------------------------------

#set xcMisc(gifsicle)  /full/path/to/gifsicle
#set xcMisc(ppmtompeg) /full/path/to/ppmtompeg


# ------------------------------------------------------------------------
#  An image conversion program: we need PPM to PNG/JPG/GIF/... conversion.
#  The \"convert\" program of ImageMagick (http://www.imagemagick.org/) is
#  a convinient choice.
# 
#  It is possible to specify the command-line options. For example:
#
#  set xcMisc(ImageMagick.convert) \"/usr/bin/convert \\
#                             -quality 90 -border 3x3 -bordercolor black\"
# 
#  We can also specify convert options separately as:
#
#  set xcMisc(ImageMagick.convertOptions) \"-quality 90 -antialias \\
#                             -blur 1x1 -trim -bordercolor white \\
#                             -border 20x20 -bordercolor black -border 3x3\"
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# How many decimal digits should the \"Measurer\" report for distances
# and angles (including dihedral). The syntax is:
#
# set select(dist_precision) number-of-decimal digits; # for distances
# set select(angl_precision) number-of-decimal digits; # for angles
# 
# Example:
#
# set select(dist_precision) 4
# set select(angl_precision) 3
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Custom setting of the atomic radii. The syntax is:
#
# set atmRad(atomic_number) radius
#
# Example:
#
# set atmRad(1) 0.5; # custom radius for Hydrogen
# set atmRad(8) 1.2; # custom radius for Oxygen
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------
# Custom setting of the atomic colors. The syntax is:
#
# set atmCol(atomic_number) {red gren blue}
#
# The components (red,gren,blue) must be in range [0,1]
#
# Example:
#
# set atmCol(1) {0.5 0.5 0.5}; # custom color for Hydrogen
# set atmCol(8) {0.0 1.0 0.0}; # custom color for Oxygen
# ------------------------------------------------------------------------



# -------------------------------------------------------------------------
# Custom setting for a variety of molecular display parameters, such as
# ball-factors, specefill scale factors, tessellation factors, etc.
#
# Below are the default values. If you would like to change the
# default for a particular parameter, then uncomment the appropriate
# line and set the value according to your needs.
# ------------------------------------------------------------------------

## spacefill scale factor
#set myParam(ATRAD_SCALE)     1.40
#
## tesselation factor
#set myParam(TESSELLATION)    15.0 
#
## RGB color of unibonds (each compoenent must be within [0,1])
#set myParam(UNIBONDCOLOR)    {1.00 1.00 1.00} 
#
## Perspective Fovy, Front and Back parameters. The smaller the Fovy
## the larger the perception of perspective. Front and Back parameters
## determine the front and back clipping planes. The smaller the Back
## parameter the more the structure is clipped from the back side. The
## Front parameter is counter-intuitive, meaning the smaller it is the
## more the structure is clipped from the front side.
#
#set myParam(PERSPECTIVEFOVY)  2.5
#set myParam(PERSPECTIVEFRONT) 0.65
#set myParam(PERSPECTIVEBACK)  3.0
#
## ball-factor
#set myParam(BALLF)           0.4
#
## rod-factor
#set myParam(RODF)            0.6
#
## line-width of wireframe display-mode (in pixels)
#set myParam(WFLINEWIDTH)     1
#
## line-width of pointline display-mode (in pixels)
#set myParam(PLLINEWIDTH)     1
#
## line-width of crystal cell's frames
#set myParam(FRAMELINEWIDTH)  1
#
## Lighting-Off outline width
#set myParam(OUTLINEWIDTH)     1
#
## Lighting-On wire line width
#set myParam(WF3DLINEWIDTH)    1
#
## point-size of pointline display-mode (in pixels)
#set myParam(PLRADIUS)        6
#
## chemical connectivity factor
#set myParam(COV_SCALE)       1.05
#
## RGB color of crystal frame (each compoenent must be within [0,1])
#set myParam(FRAMECOL)        {0.88 1.00 0.67} 
#
## line-width of crystal frame
#set myParam(FRAMELINEWIDTH)  1 
#
## rod-factor of crystal frame
#set myParam(FRAMERODF)       0.1
#
## RGB background of XCRYSDEN display window 
## (each compoenent must be within [0,1])
#set myParam(BACKGROUND)      {0.00 0.00 0.00}
#
## maximum number of cells per direction for CRYSTALs
#set myParam(CRYSTAL_MAXCELL)  10
#
## maximum number of cells per direction for SLABs
#set myParam(SLAB_MAXCELL)     20
#
## maximum number of cells for POLYMERs
#set myParam(POLYMER_MAXCELL)  50
#
## default atomic-label's font (in X11 XLFD format)
#set myParam(ATOMIC_LABEL_FONT) -adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1
#
## default atomic-label's bright and dark color (in clamped-float RGB format)
#set myParam(ATOMIC_LABEL_BRIGHTCOLOR) {1.0 1.0 1.0}
#set myParam(ATOMIC_LABEL_DARKCOLOR)   {0.0 0.0 0.0}
#
#
## this are the parameters for the \"mpeg_encode\" program:
##--BEGIN::
#set myParam(MPEG_ENCODE_PARAM_FILE) {
#PATTERN          IBBPBBPBBPBBPBBP
#OUTPUT           \$output_file
#BASE_FILE_FORMAT PPM
#INPUT_FORMAT     UCB
#\$input_convert
#GOP_SIZE         16
#SLICES_PER_FRAME 1
#INPUT_DIR        \$input_dir
#INPUT
#\$input_files
#END_INPUT
#PIXEL           FULL
#RANGE           10
#PSEARCH_ALG     LOGARITHMIC
#BSEARCH_ALG     CROSS2
#IQSCALE         8
#PQSCALE         10
#BQSCALE         25
#REFERENCE_FRAME ORIGINAL
#BIT_RATE        1000000
#BUFFER_SIZE     327680
#FRAME_RATE      30
#}
##--END
#
# ------------------------------------------------------------------------



#
# ------------------------------------------------------------------------
#  Here go custom user-specified options
# ------------------------------------------------------------------------
#  In order to add an --unknown option to the Xcrysden allowed options,
#  do the following
#  
# Usage: 
#       addOption option converterProgram description
#
# Arguments:
#       option           ... option to add to XCRYSDEN options
#       converterProgram ... program that converts from an "unknown" to XSF format;
#                            this program must be supplied by the user !!!
#       description      ... description of the options that will appear in the
#                            help message (i.e. xcrysden --help).
# Example: 
#       addOption --unknown /home/tone/utils/unknown2xsf {
#               load structure from unknown file format
#       }
#
" >> $tempdir/customFile.$$

    if [ "x$cXX" = "xy" ]; then
	echo "

   You have specified the following definitions:
   ---------------------------------------------
"
	echo "   CRYSTAL's crystal module:      $cry
   CRYSTAL's properties module:   $pro
"
    fi

    GetYesNoAnswer "   Is this correct" cor

    if [ \( "$cor" != "y" \) -a \( "$cor" != "yes" \) ]; then
	$CLEAR
	echo "
   You answered NO. 

   Therefore, please answer again to all the questions !!!
"
	CustomDef
	return
    fi

    $CLEAR

    if test -f $HOME/.xcrysden/custom-definitions ; then
	backup_name=`BackupName $HOME/.xcrysden/custom-definitions`
	echo "
   backing-up existing custom-definitions file to $backup_name"
	cp $HOME/.xcrysden/custom-definitions $backup_name
    fi

    $ECHO_n "   
   creating new $HOME/.xcrysden/custom-definitions file ... "
    cp $tempdir/customFile.$$ $HOME/.xcrysden/custom-definitions
    if [ $? -eq 0 ]; then
	echo OK
    else
	echo "
   ERROR: failed to create $HOME/.xcrysden/custom-definitions file !!!"
    fi

    echo "
   Press <Enter> to continue ..."
    read ans

    $CLEAR

#    $ECHO_n "
#Would you like to edit ~/.xcrysden/custom-definitions file now ([y]es/[n]o):"
#    read edit
#    if [ \( "$edit" = "y" \) -o \( "$edit" = "yes" \) ]; then
#	${EDITOR:-vi} $HOME/.xcrysden/custom-definitions
#    fi


    if test -d $HOME/Desktop; then	
	GetYesNoAnswer "   Would you like to add xcrysden icon to the desktop" desktop
	if test "$desktop" = "y"; then
	    $ECHO_n "   
   creating xcrysden dekstop icon via $HOME/Desktop/xcrysden.desktop file ... "
	    echo "
[Desktop Entry]
Encoding=UTF-8
Name=xcrysden
Exec=$XCRYSDEN_TOPDIR/xcrysden
Type=Application
Icon=$XCRYSDEN_TOPDIR/images/xcrysden.png
" > $HOME/Desktop/xcrysden.desktop
	    if [ $? -eq 0 ]; then
		echo OK
	    else
		echo "
   ERROR: failed to create $HOME/Desktop/xcrysden.desktop file !!!"
	    fi
	fi
    fi


    #
    # print message about the programs xcrysden uses
    #
    echo "

   XCRYSDEN uses a few external utility programs. These are the following:

   1.) the \"Open Babel\" program is needed to display molecular structure from 
       Gaussian Z-matrix file.

   2.) the ImageMagick's suit of programs is used to print-to-file (PNG/JPG/GIF/...), 
       to create animated-GIFs, and to screen dump windows.

   3.) the \"Mplayer's mencoder\" program is used to create AVI/MPEG movies
       from the animation snapshots.

   4.) xcrysden may also use other programs such as gifsicle, xwd, ppmtompeg, ...

   Links
   --
   Open BABEL:  http://openbabel.org/
   IMAGEMAGICK: http://www.imagemagick.org/
   MENCODER:    http://www.mplayerhq.hu/

   Press <Enter> to continue ...
" | $MORE
    read ans
}


# ------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------

if [ \( -d /tmp \) -a \( -w /tmp \) ]; then
    tempdir=/tmp
else
    # set tempdir to $HOME
    tempdir=$HOME
fi

CustomDir

echo "   

   Answer the following question(s) please:
   ----------------------------------------"

if [ ! -d $HOME/.xcrysden ]; then
    mkdir $HOME/.xcrysden
fi

CustomDef

if [ -f $tempdir/customFile.$$ ]; then
    rm -f $tempdir/customFile.$$
fi

$CLEAR
echo "

   Please report BUGS to Tone.Kokalj@ijs.si. 

   TERMS OF USE:
   -------------
   XCRYSDEN is released under the GNU General Public License.
   
   Whenever graphics generated by XCrySDen are used in scientific
   publications, it shall be greatly appreciated to include an
   explicit reference. The preferred form is the following:

   [ref] A. Kokalj, Comp. Mater. Sci., Vol. 28, p. 155, 2003.
         Code available from http://www.xcrysden.org/.
   
      
   XCRYSDEN tips: 
                  - try: xcrysden --help

                  - occasionally clean the xcrysden scratch directory
                    (you can use the \"xc_cleanscratch\" utility)

   Press <Enter> to continue ...
" | $MORE
read ans
$CLEAR

if [ \( "$xcv" != "" \) ]; then
    $CLEAR
    echo "
BEWARE: please edit manually your profile (i.e. $profile) to make the 
        XCRYSDEN section as:
"| cat - $install_file | $MORE
fi
  
echo "

FINAL NOTICE: before running xcrysden please load (source) your profile, i.e.:

   - for bash:           source ~/.bashrc
   - for csh or tcsh:    source ~/.cshrc
   - for sh:             . ~/.profile

Then type: xcrysden
"
