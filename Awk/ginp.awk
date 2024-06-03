#!/bin/sh

cat "$1" | awk '
BEGIN {endline=0; n=1}
($1 ~ /OPTCOOR/ || $1 ~ /OPTBERNY/ || $1 ~ /OPTGEOM/)  && (n != 1) {if ( endline == 0 ) {print "STOP"}
                           endline=1}
($1 ~ /END/)  && (n != 1) {if ( endline == 0 ) {print "STOP"}
                           endline=1}
($1 ~ /STOP/) && (n != 1) {if ( endline == 0 ) {print "STOP"}
                           endline=1}
/a*/ { n++; if ( endline == 0 ) {print} }
'
