#!/bin/sh

# ------------------------------------------------------------------------
#****** ScriptingExamples/plot_all_MO.sh ***
# NAME
# plot_all_MO.sh --- plots all molecular orbitals (MO) in alike manner
#
# USAGE
# ./plot_all_MO.sh
#
# WARNING
# do not execute this script as: xcrysden -s plot_all_MO.sh 
#                                (this will not work)
#
# PURPOSE
# This shell scripts shows how one can plot several (or all) molecular
# orbitals of a given molecule in an automatic fashion, where the
# display parameters for all MOs are the same. 
#
# To construct the script as this one the following steps have to be
# done:
# 1.) the 1st MO has to be displayed in xcrysden in usual way, and
#     when the display is properly set, then save it as: 
#     File->Save State and Structure.
#
# 2.) encapsulate such saved file in s simple shell-script as this one
#
#
# 3.) replace the STRUCTURE-PART of the saved file, and request
#    therein to load a given file instead (see below)
#
#
# 4.) replace all $ characters with \$ in the saved script except for
#     the $ character for the $file variable
#
# 5.) at the end of the saved file add: scripting::printToFile; exit
#
# 6.) exucute the shell-script
#
#
# AUTHOR
# Anton Kokalj
#
# CREATION DATE
# Fri May 13 11:07:01 CEST 2005
#
# SOURCE


# this is the Bourne-chell script

for file in CO_homo.xsf.gz CO_lumo.xsf.gz 
do
    xcrysden --xsf $file --script MO-state.xcrysden --print ${file%.xsf.fz}.png
done

#****
# ------------------------------------------------------------------------
