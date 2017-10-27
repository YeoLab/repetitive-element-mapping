#!/usr/bin/env bash

mkdir genomeDir
tar -zxvf $4 --strip-components 1 -C ./genomeDir;

filename=$(ls genomeDir/*.fa)
extension="${filename##*.}"
filename="${filename%.*}"

perl \
$1 \
$2 \
$3 \
$filename \
$5 \
$6