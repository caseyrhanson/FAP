#!/bin/sh

set -e
set -u
set -o pipefail
set -x

module load dx-toolkit
module load dnanexus_ua/1.5.20

for  i in Macrogen/* 
do
  if [[ -d $i ]]
  then  
    project="project-F7vFb2j0ZgxJ3BPB18fkxKQG" #	FAP-DNAseq - Oct 31, 2017)
    folder=$(basename $i)
    fqs=$(ls ${i}/*.fastq.gz)
    fqs=($fqs)
    ${fqs[@]:0:1}
    ua -p ${project} -f "/${folder}" ${fqs[@]} --do-not-compress
  fi
done
