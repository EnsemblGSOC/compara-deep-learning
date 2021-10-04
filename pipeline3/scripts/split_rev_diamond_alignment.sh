#!/bin/sh

input=$1
output=$2
threads=$3
out_dir=$4
SPECIES=$5

mkdir -p $out_dir
echo $out_dir
less $input | cut -d $'\t' -f 1 | cut -d . -f 1 | uniq > $output
xargs -n 1 -P $threads -a $output -I string sh -c "echo string; zgrep string $input | gzip > ${out_dir}string_${SPECIES}_rev.outs.tsv.gz"