#!/bin/bash

cd /home/ahorning/DNAseq/DNAseq_fq_VCF/ #VCFComparison

#create arrays with important file and folder name information. each should have the same number of files
# This will only work if the file paths have no spaces in them
array=($(find "$(pwd)" -name "*.vcf.gz" | grep -v ".g.vcf.gz" | sort)) ; echo ${#array[@]}
pathtofile=($(find "$(pwd)" -name "*.vcf.gz" | grep -v ".g.vcf.gz" | sort)); echo ${#pathtofile[@]}
folder=($(find "$(pwd)" -name "*.vcf.gz" | grep -v ".g.vcf.gz" | sort | cut -d "/" -f 6)); echo ${#folder[@]}
file=($(find "$(pwd)" -name "*.vcf.gz" | grep -v ".g.vcf.gz" | sort | cut -d "/" -f 7)); echo ${#file[@]}
base=($(find "$(pwd)" -name "*.vcf.gz" | grep -v ".g.vcf.gz" | sort | cut -d "/" -f 7 | cut -d "." -f 1)); echo ${#base[@]}

module load bcftools/1.3.1

for index in ${!array[*]}; do 
  mkdir "VCFComparison/${folder[$index]}"
  bcftools stats "${pathtofile[$index]}" > "VCFComparison/${folder[$index]}/${folder[$index]}_vcfstats.txt"
  plot-vcfstats "VCFComparison/${folder[$index]}/${folder[$index]}_vcfstats.txt" -p "VCFComparison/${folder[$index]}/"
done


# array=(  Vietnam  Germany  Argentina)
# array2=(  Asia  Europe  America)

# for index in ${!array[*]}; do 
#   echo "${array[$index]} is in ${array2[$index]}"
# done

# Vietnam is in Asia
# Germany is in Europe
# Argentina is in America


# # array=( $(find blah -mindepth 3 -maxdepth 3 -type d -regex ".*/V[0-9]+/wsdls+") )

# # # loop over it
# # for i in ${array[@]}
# # do
# #     echo $i
# # done

