#!/bin/bash


# Filters large database-annotated VCF/TXT files by a mutations Pathogenicity (CLINSIG), and predicted deleterious effect based on fathmm-MKL and DANN algorithms

for i in EP-107-NL EP-11 EP-26 EP-31 EP-57 EP-6 EP-74B EP-84 EP-88B EP-AdenoCa EP_Dec_NL ES24 ES-3 ES_Dec_NL JC-AdenoCa JC-ASC JC-SIG JP31 JP34B JP38 JP3 JP61B JP63 JP69 JP6B JP9 JPAdenoCa JP_Dec_NL-1 NBM-R-10 NBM-R-1; do
	awk '$12 ~ /Pathogenic/ || $60 ~ /D/ || $56 > 0.95' /home/ahorning/DNAseq/DNAseq_fq_VCF/${i}/${i}_myanno.hg38_multianno.txt | cut -f1-17,56,60,95,96 > /home/ahorning/DNAseq/DNAseq_fq_VCF/VCFannotation/${i}_anno_filtered.txt
done


#cut -f1-17,56,60 EP-107-NL_myanno.hg38_multianno.txt

# 60 is fathmm-MKL_coding_pred == D
# 12 is CLINSIG == Pathogenic
# 38 is FATHMM_pred==D
# 60 is fathmm-MKL_coding_pred == D
# 56 is DANN_score < 0.950


#FATHMM information
#https://brb.nci.nih.gov/seqtools/colexpanno.html

#
#http://fathmm.biocompute.org.uk/fathmmMKL.htm



