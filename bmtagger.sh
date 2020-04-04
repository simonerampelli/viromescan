#!/bin/bash

set -o pipefail

TEMP=`getopt -o hV1:2:o:T:b:d:q:x:A:C:X --long help,debug,ref:,extract,old-srprism -n 'bmtagger' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

eval set -- "$TEMP"

accession=''
quality=0
reference=''
blastdb=''
srindex=''
bmfiles=''
input1=''
input2=''
inp1='' # to extract unmasked fasta
inp2='' # to extract unmasked fasta
output=/dev/stdout
extract=''
done=0
debug=0

#--sra-read-clipinfo"
#bmoptions="-TP -lN5 -z -L 1000000 -c 32 -s 4"
bmoptions="-TP -lN5 -z -L 1000000 -c 64 -s 16"
bmoptions_oldsrprism="-TP -lN5 -z -L 1000000 -c 32 -s 4"
bmoptions_debug="-R"
bmoptions_debug_sra="--dump-as-fasta"
blastnopts="-outfmt 6 -word_size 16 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -use_index true"
srprismopts_oldsrprism="-b 100000000 -n 2 -R 0 -r 1 -M 7168"
srprismopts="-b 100000000 -n 5 -R 0 -r 1 -M 7168 -p false"
	#--trace-level info"
	#-M 15360 

: ${TMPDIR:=/tmp}
: ${BMFILTER:=bmfilter}
: ${BLASTN:=blastn}
: ${SRPRISM:=srprism}
: ${EXTRACT_FA:=extract_fullseq}

function check_exec () {
	if [ "$1" == "" ] ; then return 0; fi
	if ! which "$1" >&2 ; then echo "FATAL: Failed to find $1" >&2 ; exit 20 ; fi
}

function show_help () {
	echo "usage: bmtagger [-hV] [-q 0|1] [-C config] -1 input.fa [-2 matepairs.fa] -b genome.wbm -d genome-seqdb -x srindex [-o blacklist] [-T tmpdir] [-X]"
	echo "usage: bmtagger [-hV] [-q 0|1] [-C config] -1 input.fa [-2 matepairs.fa] --ref=reference [-o blacklist] [-T tmpdir] [-X]"
	echo "usage: bmtagger [-hV] [-q 0|1] [-C config] -A accession [--ref=reference] [-b genome.wbm] [-d genome-seqdb] [-x srindex] [-T tmpdir]"
	echo "use --ref=name to point to .wbm, seqdb and srprism index if they have the same path and basename"
	echo "use --extract or -X to generate fasta or fastq files which will NOT contain tagged sequences (-o required)"
	echo "use --debug to leave temporary data on exit"
	echo "use --old-srprism to use options for older version of srprism (interferes with config file)"
	done=1
}

function show_version () {
	echo "version 1.1.0"
	done=1
}

function finalize () {
	if [ $debug == 0 ] ; then rm -fr "$TMPDIR"/bmtagger.$tmpstr.* ; fi
	if [ "$1" == "" ] ; then exit 100 ; else exit "$1" ; fi
}

function parse_config () {
	local required="$1"
	local file="$2"
	if [ -e "$file" ] ; then
		echo "Config: Reading $file" >&2
		source $file
		rc="$?"
		if [ $rc != 0 ] ; then
			echo "Error: Failed to read file $file" >&2
			exit 22
		fi
	elif [ $required != 0 ] ; then
		echo "Error: failed to read $file" >&2
		exit 23
	else 
		echo "Info: no $file found" >&2
	fi
}

parse_config 0 ./bmtagger.conf

while true ; do
	case "$1" in
	-h|--help) show_help ; shift ;;
	-V) show_version ; shift ;;
	--debug) debug=1 ; shift ;;
	--ref) reference="$2" ; shift 2 ;;
	-1) input1="-1$2" ; inp1="$2" ; shift 2 ;;
	-2) input2="-2$2" ; inp2="$2" ; shift 2 ;;
	-o) output="$2" ; shift 2 ;;
	-A) if [ "$accession" == "" ] ; then accession="$2" ; else echo "ERROR: Can't use -A multiple times" >&2 ; exit 1 ; fi ; shift 2 ;;
	-C) parse_config 1 "$2" ; shift 2 ;;
	-b) bmfiles="$bmfiles -b$2" ; shift 2 ;;
	-d) blastdb="$2" ; shift 2 ;;
	-x) srindex="$2" ; shift 2 ;;
	-q) quality="$2" ;  shift 2 ;;
	-T) TMPDIR="$2" ; shift 2 ;;
	-X) extract=1 ; shift ;;
	--extract) extract=1; shift ;;
	--old-srprism) bmoptions=$bmoptions_oldsrprism; srprismopts=$srprismopts_oldsrprism; shift ;;
	--) break ;;
	*) echo "Unknown option ``$1 $2''" >&2 ; exit 1 ;;
	esac
done

echo "Using following programs:" >&2
check_exec "$BMFILTER"
check_exec "$SRPRISM"
check_exec "$BLASTN"
check_exec "$EXTRACT_FA"

if [ $done == 1 ] ; then exit 0 ; fi
if [ $debug != 0 ] ; then 
	bmoptions="$bmoptions $bmoptions_debug" 
	if [ "$accession" != "" ] ; then bmoptions="$bmoptions $bmoptions_debug_sra" ; fi
fi 

if [ ! -d "$TMPDIR" ] ; then
	echo "FATAL: $TMPDIR is not directory" >&2 
	exit 21 
fi

spotId_only=""

if [ "$accession" != "" ] ; then
	if [ -n "$extract" ] ; then
		echo "ERROR: options -A and -X (or --extract) can't be used together" >&2
		exit 1
	fi
	if [ "$input1:$input2" != ":" ] ; then
		echo "ERROR: One should not use -A with -1, -2" >&2 
		exit 1 
	fi
	if [ "$output" == "/dev/stdout" ] ; then output=. ; fi
	test -d "$output" || mkdir -p "$output"
	test -d "$output" || { echo "ERROR: failed to create directory [$output]" >&2 ;  exit 100 ; }
	output="$output/"$(basename "$accession" .sra)".blacklist"
	accession="-A$accession"
	spotId_only="-I"
fi

case "$output" in
	-) output="/dev/stdout" ; tmpout="$output" ;;
	/dev/*) tmpout="$output" ;;
	*) tmpout="$output~" ;;
esac

if test -n "$extract" ; then 
	case "$output" in 
		/dev/*) echo "ERROR: -o is required and can't point to /dev/* if -X (or --extract) is used" >&2 ; finalize 1 ;;
	esac
fi

test -z "$srindex" && srindex="$reference"
test -z "$blastdb" && blastdb="$reference"
test -z "$bmfile"  && bmfile="$reference".wbm

echo "MAIN SCRIPT IS $0 (PID=$$)" >&2
trap finalize INT TERM USR1 USR2 HUP

tmpstr=`date '+%s'`.`hostname -s`.$$

if test -z $bmfiles ; then
	echo "FAILED: bmfilter needs bitmask files as argument" >&2
	finalize 1
fi

echo "RUNNING bmfilter" >&2
#echo "$BMFILTER $accession $input1 $input2 -q $quality $bmfiles $spotId_only -o "$TMPDIR"/bmtagger.$tmpstr $bmoptions"
time $BMFILTER $accession $input1 $input2 -q $quality $bmfiles $spotId_only -o "$TMPDIR"/bmtagger.$tmpstr $bmoptions
rc=$?

if [ $rc != 0 ] ; then
	echo "FAILED: bmfilter with rc=$rc" >&2 ; 
	finalize 2
fi

function align_long () {
	test -s "$1".fa || return 0
	echo "RUNNING align_long for '$1'" >&2
	time $BLASTN \
		-task megablast \
		-db "$blastdb" \
		-query "$1".fa \
		-out   "$1".bn \
		-index_name "$blastdb" \
		$blastnopts
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: blastn for $1" >&2 ; finalize 3 ; fi
	## blastn filter criteria: $4 = hitLength, $3 = %id, $1 = readID
	time awk '($4 >= 90 || ($4 >= 50 && $3 >= 90)) { print $1 }' "$1".bn > "$1".lst
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: awk for blastn results for $1" >&2 ; finalize 4 ; fi
	return $rc
}

function align_short () {
	test -s "$1".fa || return 0
	echo "RUNNING align_short for '$1'" >&2
	time $SRPRISM search \
		-I "$srindex" \
		-i "$1".fa \
		-o "$1".srprism \
		-T "$TMPDIR" \
		-O tabular \
		$srprismopts
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: srprism for $1" >&2 ; finalize 5 ; fi
	## srprism filter criteria: everything, $2 = readID
	time awk '{ print $2 }' "$1".srprism > "$1".lst
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: awk for srprism results for $1" >&2 ; finalize 6 ; fi
	return $rc
}

function append () {
	test -s "$2" || return 0
	time uniq "$2" >> "$1"
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: cat $2 >> $1 with rc=$?" >&2 ; finalize 7 ; fi
	return $rc
}

function extract_fa () {
	test -z "$EXTRACT_FA" && return 0
	test -s "$1".lst || return 0
	test -s "$1"2.fa || return 0
	time "$EXTRACT_FA" "$1".lst -remove -fasta -single "$1"2.fa > "$1"2x.fa
	rc=$?
	if [ $rc != 0 ] ; then echo "FAILED: $EXTRACT_FA with rc=$?" >&2 ; finalize 8 ; fi
	return $rc
}

align_short "$TMPDIR"/bmtagger.$tmpstr.short
extract_fa  "$TMPDIR"/bmtagger.$tmpstr.short
align_short "$TMPDIR"/bmtagger.$tmpstr.short2x

align_long "$TMPDIR"/bmtagger.$tmpstr.long
extract_fa "$TMPDIR"/bmtagger.$tmpstr.long
align_long "$TMPDIR"/bmtagger.$tmpstr.long2x

awk '($2 == "H") { print $1 }' "$TMPDIR"/bmtagger.$tmpstr.tag >> "$tmpout"
rc=$?
if [ $rc != 0 ] ; then echo "FAILED: awk for tagfile" >&2 ; finalize 8; fi

append "$tmpout" "$TMPDIR"/bmtagger.$tmpstr.short.lst
append "$tmpout" "$TMPDIR"/bmtagger.$tmpstr.short2x.lst
append "$tmpout" "$TMPDIR"/bmtagger.$tmpstr.long.lst
append "$tmpout" "$TMPDIR"/bmtagger.$tmpstr.long2x.lst

function extract_result() {
	case "$1" in
	*.gz) time gzip -dc "$1" | "$EXTRACT_FA" "$tmpout" -remove "$3" "$4" /dev/stdin > "$2"~ ;;
	*.bz2) time bzip2 -dc "$1" | "$EXTRACT_FA" "$tmpout" -remove "$3" "$4" /dev/stdin > "$2"~ ;;
	*) time "$EXTRACT_FA" "$tmpout" -remove "$3" "$4" "$1" > "$2"~ ;;
	esac
	rc=$?
	if [ $rc != 0 ] ; then 
		echo "FAILED to extract result to $2" >&2
		finalize 9
	fi
	mv "$2"~ "$2"
	if [ $rc != 0 ] ; then 
		echo "FAILED to move [$2~] to [$2]" >&2
		finalize 9
	fi
}

if test -n "$extract" ; then
	if [ "$quality" == 0 ] ; then mode="-fasta" ; suffix=".fa" ; else mode="-fastq" ; suffix=".fastq" ; fi
	if test -n "$inp2" ; then
		id1=$(head -1 "$inp1" | awk '{print $1}')
		id2=$(head -1 "$inp2" | awk '{print $1}')
		if test "$id1" == "$id2" ; then
			mate1=-single
			mate2=-single
		else
			mate1=-mate1
			mate2=-mate2
		fi
		extract_result "$inp1" "${output}_1$suffix" "$mode" "$mate1"
		extract_result "$inp2" "${output}_2$suffix" "$mode" "$mate2"
	else
		extract_result "$inp1" "$output$suffix" "$mode" -single
	fi
	rm "$tmpout"
else
	if test "$tmpout" != "$output" ; then
		mv "$tmpout" "$output"
		rc=$?
		if [ $rc != 0 ] ; then 
			echo "FAILED to move [$tmpout] to [$output]" >&2 
			finalize 9
		fi
	fi
fi

echo "DONE $0 (PID=$$)" >&2
finalize 0

