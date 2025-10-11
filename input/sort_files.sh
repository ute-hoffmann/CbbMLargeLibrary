#!/bin/bash
OUTDIR="./fastq/" # files will be created to this directory
mkdir $OUTDIR # out directory is created
dir="./fastq/illumina" # script assumes directories with sequencing files in this directory
for file in ${dir}/*/
  do
    file="${file%%/}"
    #echo "${file##.*/}"
    for data in $file/*
       do
#          data="${data}"
          echo "${data}"
          mv ${data} $OUTDIR
       done
  done

