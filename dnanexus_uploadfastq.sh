#!/bin/sh

'''
This script uploads fastq files from the SCG4 folders into specific folders in DNAnexus based on the folder names in my DNAseq folder

I added if/then statments to test for existence of the files before uploading them too
'''

# Order of script is important

set -e
set -u
set -o pipefail
set -x

cd ~
for sample in JP69 ; do # $(ls DNAseq/Macrogen/ | grep -v nohup.out | grep -v 1508KHX-0009_EP-107-NL.pdf | grep -v EP-107-NL | grep -v EP-11 | grep -v EP-26 ); do
	dx cd ..
	dx tree > uploaded_samples.txt
	dx cd /${sample}
	FASTA1="${sample}_R1.fastq.gz"
	FASTA2="${sample}_R2.fastq.gz"
		echo dx upload --wait -r /home/ahorning/DNAseq/root/${sample}/01_READS/$FASTA1  |  bash
		echo dx upload --wait -r /home/ahorning/DNAseq/root/${sample}/01_READS/$FASTA2  |  bash
done

# if [ grep $FASTA1 uploaded_samples.txt ]; then
# 		echo dx upload --wait -r /home/ahorning/DNAseq/root/${sample}/01_READS/$FASTA1 #| bash
# 	fi

# 	if [ grep $FASTA2 uploaded_samples.txt ]; then
# 		echo dx upload --wait -r /home/ahorning/DNAseq/root/${sample}/01_READS/$FASTA2 #| bash
# 	fi

# dx find data -name EP-107-NL_R1.fastq.gz


# test [ -f /EP-107-NL/EP-107-NL_R1.fastq.gz ] && echo "Found" || echo "Not Found"

# grep   -l, --files-with-matches  print only names of FILEs containing matches

# grep EP-107-NL_R1.fastq.gz uploaded_samples.txt && echo "found"
# grep NBM-R-10_R1.fastq.gz uploaded_samples.txt || echo "not found"

# #tests if file exists
# file="/etc/hosts"
# if [ -e "$file" ]; then
# 	echo "$file found."
# else
# 	echo "$file not found."
# fi

# while [! -f $FASTA1]; do
# 	#statements
# done

# if [ "$PASSWORD" == "$VALID_PASSWORD" ]; then
# 	echo "You have access!"
# else
# 	echo "ACCESS DENIED!"
# fi





