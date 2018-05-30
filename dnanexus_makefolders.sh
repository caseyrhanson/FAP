#!/bin/sh

'''
This script creates folders in DNAnexus based on the folder names in my DNAseq folder
'''

# Order of script is important

set -e
set -u
set -o pipefail


# source ~/.bashrc
# dxtk
cd /home/ahorning/DNAseq
dx select --level ADMINISTER project-F7vFb2j0ZgxJ3BPB18fkxKQG
for sample in $(ls Macrogen | grep -v nohup.out | grep -v 1508KHX-0009_EP-107-NL.pdf); do
	echo dx mkdir ${sample} | bash
	#echo dx mkdir $(find . -name '/home/ahorning/DNAseq/root/${sample}/01_READS/${sample}_R{1,2}.fastq.gz')
done







