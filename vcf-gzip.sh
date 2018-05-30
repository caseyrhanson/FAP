#!/bin/bash

#This code will re zip all of the vcf files of interest.

# array=($(find . -iname "*.vcf.gz" | grep -v "g.vcf.gz" | sort)) ; echo ${#array[@]}
# pathtofile=($(find "$(pwd)" -iname "*.vcf.gz" | grep -v "g.vcf.gz" | sort)); echo ${#pathtofile[@]}
# folder=($(find "$(pwd)" -name "*.tnscope.vcf.gz" | sort | cut -d "/" -f 6)); echo ${#folder[@]}
# file=($(find "$(pwd)" -name "*.tnscope.vcf.gz" | sort | cut -d "/" -f 7)); echo ${#file[@]}
# base=($(find "$(pwd)" -name "*.tnscope.vcf.gz" | sort | cut -d "/" -f 7 | cut -d "." -f 1)); echo ${#base[@]}

cd /home/ahorning/DNAseq/DNAseq_fq_VCF

# dir="/home/ahorning/DNAseq/DNAseq_fq_VCF/VCFComparison"
files=($(find "$(pwd)" -iname "*.vcf.gz" | grep -v "g.vcf.gz" | sort))

module load bcftools/1.3.1

for i in ${!files[*]}; do 
	bcftools view "${files[$i]}" -Oz -o "${files[$i]}" # -Oz outputs file format as zipped .vcf. -o is output name
done


# -c, --collapse snps|indels|both|all|some|none
# see Common Options
# -C, --complement
# output positions present only in the first file but missing in the others



# for testing the zipped-ness of a file
# bgzip file.vcf       # or:   bcftools view file.vcf -Oz -o file.vcf.gz
# tabix file.vcf.gz    # or:   bcftools index file.vcf.gz
# Another useful command is htsfile file.vcf.gz which prints the actual file type.

# for copying and rezipping files

# mv file.vcf.gz plain.vcf
# bcftools view -Oz -o compressed.vcf.gz plain.vcf
# htsfile compressed.vcf.gz
# bcftools index compressed.vcf.gz

