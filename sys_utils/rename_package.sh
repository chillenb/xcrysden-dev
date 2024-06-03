#!/bin/sh

# Usage: $0 package.tar.gz [platform]
#
# Purpose: rename xcrysden tar.gz package to proper name for the
# www.xcrysden.org/download

Usage() {
    echo "
Usage: $0 package.tar.gz [platform]
"
    exit 1
}

nargc=$#
orig_package=$1
platform=$2

prefix=xcrysden

if test $nargc -lt 1; then
    Usage
fi

# strip XCrySDen- from the $orig_package
tmp1=${orig_package#xcrysden-}

# check if the package is a binary one
is_bin=`echo $orig_package | grep bin`

if test ! -z $is_bin; then
    # it's binary package
    version=${tmp1%-bin*.tar.gz}

    tmp2=${tmp1#${version}-bin-}
    bin_type=${tmp2%.tar.gz}

    new_package="$prefix-$version-$platform-$bin_type.tar.gz"

    if test $nargc != 2; then
	echo "
ERROR: to rename the binary package, please specify the platform"
	Usage
    fi
else
    # it's source package, simply rename XCrySDen- to $prefix-
    new_package=$prefix-$tmp1
fi
  
echo "
Package renamed to: $new_package
"
mv $orig_package $new_package