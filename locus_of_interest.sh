#!/bin/bash


# This script will interatively go through the large bam files and select portions specific to genetic regions of interest.
# With these smaller regions, I can put them through IGV to create images showing read plots
#input bam files. output is smaller bam files focused on specific genetic loci

cd /home/ahorning/DNAseq/DNAseq_fq_VCF
module load samtools/1.5
#here is a list of the full paths "*markdup.realigned.recalibrated.bam" files in order
path=($(find "$(pwd)" -name "*markdup.realigned.recalibrated.bam" | sort)) ; echo ${#path[@]}
path_nobam=($(find "$(pwd)" -name "*markdup.realigned.recalibrated.bam" | sort | cut -d "." -f 1-4)) ; echo ${#path_nobam[@]}
#sample_gene=($(find "$(pwd)" -name "*markdup.realigned.recalibrated.bam" | sort | cut -d "/" -f 6)) ; echo ${#path_nobam[@]} # sample name

#For APC
#gene="APC"
#locus="chr5:112737888-112846239"

#For TP53
gene="TP53"
locus="chr17:7668402-7687550"

#For KRAS
#gene="KRAS"
#locus="chr12:25209431-25250803"

for index in ${!path[*]}; do 
  # echo "${path[$index]}" 
  # echo "${path_nobam[$index]}"
  # echo "${path_nobam[$index]}_${gene}.bam"
  # echo "${path_nobam[$index]}_${gene}_depth.txt"
	samtools view -b "${path[$index]}" "${locus}" > "${path_nobam[$index]}_${gene}.bam"
	samtools index "${path_nobam[$index]}_${gene}.bam" #this creates a .bam.bai file from the .bam file i just made
	samtools depth "${path_nobam[$index]}_${gene}.bam" > "${path_nobam[$index]}_${gene}_depth.txt"
	awk '{ total += $3; count++ } END { print total/count }' "${path_nobam[$index]}_${gene}_depth.txt" #calculates depth of reads
done
