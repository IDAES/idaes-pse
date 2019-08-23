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

cmd="$(dirname $0)/annotate_source ${flag_dry} ${flag_rm}"

for subdir in idaes apps
do
    path="${rdir}/${subdir}"
    printf "Modifying files below $path ...\n"
    $cmd "$path" "*.py" "~__init__.py"
    printf "DONE\n"
done
