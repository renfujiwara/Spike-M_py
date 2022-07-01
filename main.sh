#!/bin/sh
mode=$1
datdir=$2
outdir=$3
seqdir=$4

if [ "$#" -ne 4 ]; then
echo " Usage: cmd [mode] [datdir] [outdir] [seqdir]"
exit 1
fi

IFS='-'; set -- $mode 

seqfnORG=$datdir"/sequence.dat"

yamlfn=$outdir"/params.yaml"

python3 main.py -s $seqfnORG -o $outdir -p $yamlfn #> $outdirC"log.txt"
