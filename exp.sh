#!/bin/sh
DATADIR="./_data/"
OUTDIR="./_out/"

while read data || [ -n "${data}" ]; do

outdir=$OUTDIR$data
datdir=$DATADIR$data
seqdir=$outdir"/seq"
mode="1-1-0"

# if [ ! -d $seqdir ]; then
#     mkdir $seqdir
# fi

nohup sh main.sh $mode $datdir $outdir $seqdir > $outdir"/log.txt" &

done < ./list_exp.txt