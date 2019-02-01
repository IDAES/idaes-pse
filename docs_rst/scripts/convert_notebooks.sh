#!/bin/bash
# Convert Jupyter notebooks in the given directory
# into a form displayable in the docs

function usage () {
    printf "usage: convert_notebooks.sh SOURCE TARGET\n"
    printf "   SOURCE is the source directory containing Jupyter notebooks\n"
    printf "   TARGET is the directory to place the converted notebooks\n"
    exit 1
}

function convert_files () {
    logfile=/tmp/dmf-jnc.log
    printf "Converting notebooks in $1\n"
    printf "Logs in %s\n" ${logfile}
    for f in ${1}/*.ipynb; do
        if [ $(basename "$f" .ipynb) != "Untitled" ]; then
            printf "  convert $f\n"
            jupyter nbconvert "$f" >${logfile} 2>&1
            o=$(dirname "$f")/$(basename "$f" .ipynb).html
            # fall back to normal python Pygments lexer
            sed 's/code:: ipython3/code:: python/g' "$o" > _foo_
            /bin/mv _foo_ "$o"
            /bin/mv "$o" "$2"
            printf "    move converted file to $2/$(basename $o)\n"
        else
           printf "Skipping Untitled notebook..\n"
        fi
    done
}

# If num args is not 2, show usage.
# Note that this will include "-h" option
[ $# -ne 2 ] && usage
# Get source and target from args
src=$1
shift
tgt=$1

convert_files $src $tgt

