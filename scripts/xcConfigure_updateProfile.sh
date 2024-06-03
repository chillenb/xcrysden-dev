# ------------------------------------------------------------------------
#   FUNCTIONS
# ------------------------------------------------------------------------
ScratchDir () {
    $ECHO_n "   Please specify XCRYSDEN_SCRATCH directory (default: $HOME/xcrys_tmp): "
    read ans
    if [ ! -d ${ans:=$HOME/xcrys_tmp} ]; then
	mkdir ${ans:=$HOME/xcrys_tmp}
    fi
    if [ $? -gt 0 ]; then
	echo "ERROR: can't create directory $ans"
	ScratchDir
    fi
    xc_scratch=${ans:=$HOME/xcrys_tmp}
}

GetProfile () {
    $ECHO_n "   Please specify your profile (default $1):"
    read ans
    if [ "$ans" = "" ]; then
	ans=$1
	if [ ! -w "$ans" ]; then
	    echo "   $ans is not writable; please specify an alternative profile"
	    GetProfile $profile
	fi
    fi
    if [ \( ! -f "$ans" \) -o \( ! -w "$ans" \) ]; then
	echo "   profile $ans does not exists or is not writable; Try again !!!"
	GetProfile $profile
    fi
    profile=$ans
}

PrintTwoProfileBlocks() {
    echo "For sh-based shells:
"
    cat $tempdir/xcInstall-sh.$$
    echo "Foc Csh-based shells:
"
    $tempdir/xcInstall-csh.$$
}


# ------------------------------------------------------------------------
#   MAIN
# ------------------------------------------------------------------------
ver=`cat $XCRYSDEN_TOPDIR/version`

$CLEAR
echo "
 
   ========================================================================

       * * * *      XCRYSDEN's Configuration/Setup Utility      * * * * 

   ------------------------------------------------------------------------

                          XCRYSDEN version: $ver

       PLEASE FOLLOW INSTRUCTIONS DURING THE CONFIGURATION PROCEDURE !!!
   =========================================================================
   
   Press <Enter> to continue ..." | $MORE
read ans
$CLEAR


#
# get the scratch directory
#
xc_scratch=$XCRYSDEN_SCRATCH
if [ "$XCRYSDEN_SCRATCH" = "" ]; then
    ScratchDir
fi
XCRYSDEN_SCRATCH=$xc_scratch


#
# get the user-shell
#
if test -z $SHELL; then
    echo "
WARNING: the SHELL variable does not exists. This is weird !!!
"
fi
user_shell=`echo $SHELL | awk '{n=split($0,shell,"/"); print shell[n]}'`


# for Csh-based shells
echo "
#------------------------------------------------------------------------
# this is for XCRYSDEN $VER; added by XCRYSDEN installation on
# $date
#------------------------------------------------------------------------
setenv XCRYSDEN_TOPDIR  $XCRYSDEN_TOPDIR
setenv XCRYSDEN_SCRATCH $xc_scratch
set path = (\$XCRYSDEN_TOPDIR \$path \$XCRYSDEN_TOPDIR/scripts \$XCRYSDEN_TOPDIR/util)
" > $tempdir/xcInstall-csh.$$

# for sh-based shells
echo "
#------------------------------------------------------------------------
# this is for XCRYSDEN $ver; added by XCRYSDEN installation on
# $date
#------------------------------------------------------------------------
XCRYSDEN_TOPDIR=$XCRYSDEN_TOPDIR
XCRYSDEN_SCRATCH=$xc_scratch
export XCRYSDEN_TOPDIR XCRYSDEN_SCRATCH
PATH=\"\$XCRYSDEN_TOPDIR:\$PATH:\$XCRYSDEN_TOPDIR/scripts:\$XCRYSDEN_TOPDIR/util\"
" > $tempdir/xcInstall-sh.$$


#
# get the profile and corresponding install_file
#
if test \( x"$user_shell" = x"csh" \) -o \( x"$user_shell" = x"tcsh" \) ; then
    # Csh and Tcsh
    profile=$HOME/.cshrc
    install_file=$tempdir/xcInstall-csh.$$
elif test x"$user_shell" = x"bash" ; then
    # Bash
    profile=$HOME/.bashrc
    install_file=$tempdir/xcInstall-sh.$$
elif test -n $user_shell ; then
    # assuming some Sh-based shell
    profile=$HOME/.profile
    install_file=$tempdir/xcInstall-sh.$$
else
    profile=""
    install_file=/dev/null
    echo "         You will have to edit your profile (\$HOME/.bashrc or 
         \$HOME/.profile or \$HOME/.csh) manually and add the 
         following:
" 
    PrintTwoProfileBlocks
fi


#
# write (if possible) to profile
#
grep_result=
if test -n $profile ; then
    if [ -f $profile ]; then
	grep_result=`grep XCRYSDEN_TOPDIR $profile`
    fi

    if [ \( "$xcv" = "" \) -o \( "$grep_result" = "" \) ]; then
        # update profile
	
	if [ -f $profile ]; then
	    if [ ! -w $profile ]; then
		echo "
   ERROR: cannot write to $profile because it is write-protected !!!
          You will have to edit your profile manually and add the 
          following:
"
		cat $install_file
	    fi
	    cp  $profile $profile.orig
	else
	    touch $profile.orig
	fi
	
	$ECHO_n "
   Updating profile file $profile ... "
	cat $profile.orig $install_file > $profile
	if [ $? -eq 0 ]; then
	    echo OK
	else
	    echo FAILED
	fi
	echo "
   The following record was inserted into $profile:
" | cat - $install_file   
    else
	echo "
========================================================================
WARNING ***
WARNING *** XCRYSDEN_TOPDIR  enviromental variable is already defined.
WARNING ***
   
You will have to edit your profile (i.e. \$HOME/.chsrc or
\$HOME/.bashrc or \$HOME/.profile) manually and make the 
XCRYSDEN section as shown below:
" | cat - $install_file | $MORE

	echo "

Press <Enter>  to continue ..."

	read ans
	$CLEAR
    fi
fi 
