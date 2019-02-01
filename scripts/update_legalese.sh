#!/usr/bin/env bash
#
# Run this to update the legalese in the header
# of source code files.

function usage () {
   printf "update_legalese.sh [-nrh]\n"
   printf "\n"
   printf "options:\n"
   printf "   -h      Print this message\n"
   printf "   -n      Do not change files, print what would change\n"
   printf "   -r      Remove existing headers\n"
   exit 1
}

flag_dry=""
flag_rm=""
case "$1" in
    -n) flag_dry="-n";;
    -r) flag_rm='-r' ;;
    -h) usage ;;
    *) ;;
esac

rdir=$(realpath $(dirname $0)/.. )

printf "Modifying legalese in files below $rdir\n"

cmd="$(dirname $0)/annotate_source ${flag_dry} ${flag_rm}"

$cmd $rdir/idaes "*.py" "~__init__.py"
# $cmd $rdir/bin "*"
# $cmd $rdir/examples "*.py"
