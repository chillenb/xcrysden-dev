
for cmd in xcrysden pwo2xsf pwi2xsf ptable unitconv
do
    cat > $cmd.wrapper << EOF 
#!/bin/sh -f

# simple wrapper to run $cmd
XCRYSDEN_TOPDIR="$prefix/share/$xcrysden"
XCRYSDEN_LIB_BINDIR="$prefix/lib/$xcrysden"
export XCRYSDEN_TOPDIR
export XCRYSDEN_LIB_BINDIR

EOF
done

cat >> xcrysden.wrapper <<EOF
. "\$XCRYSDEN_TOPDIR/xcrysden"
EOF

cat >> pwo2xsf.wrapper <<EOF
. "\$XCRYSDEN_TOPDIR/scripts/pwo2xsf.sh"
EOF

cat >> pwi2xsf.wrapper <<EOF
. "\$XCRYSDEN_TOPDIR/scripts/pwi2xsf.sh"
EOF

cat >> ptable.wrapper <<EOF
"\$XCRYSDEN_TOPDIR/util/ptable"
EOF

cat >> unitconv.wrapper <<EOF
"\$XCRYSDEN_TOPDIR/util/unitconv" \$@
EOF
