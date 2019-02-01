#!/bin/bash
describe_program () {
	printf "Build documentation in development mode.\n"
	printf "Runs a clean build of the HTML documentation.\n"
	printf "\n"
	printf "Output is redirected and warnings and errors\n"
	printf "are split into files 'docs.warnings' and\n"
	printf "'docs.errors', respectively.\n"
}

usage () {
	printf "usage: $0 [-hk]\n"
	printf "  -h    print this help message\n"
	printf "  -k    keep docs.stderr and docs.stdout\n"
	exit 0
}
keep=0
while [ "x$1" != x ]; do
	a=$1
	shift
	case $a in
		-h)	# help
			describe_program
			printf "\n"
			usage
			;;
		-k)	# keep
			keep=1
			;;
		*)
			printf "unrecognized argument '$a'\n";
		   	usage
			;;
	esac
done

make clean >/dev/null
printf "Rebuilding API docs\n"
make apidoc > docs.stdout 2> docs.stderr
printf "Extracting schemas\n"
./extract_schemas.py
printf "Rebuilding HTML pages\n"
make html >> docs.stdout 2>> docs.stderr

grep WARNING docs.stderr > docs.warnings
numwarn=$( wc -l docs.warnings | cut -f1 -d' ' )
grep ERROR docs.stderr > docs.errors
numerr=$( wc -l docs.errors | cut -f1 -d' ' )
[ $keep -eq 1 ] || rm -f docs.stderr docs.stdout
if [ $numerr -eq 0 ]; then
	printf "No errors\n"
	rm -f docs.errors
else
	printf "%d errors in docs.errors\n" $numerr
	cat docs.errors
	printf "----\n"
fi
printf "%s warnings in docs.warnings\n" $numwarn
