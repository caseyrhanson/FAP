#!/bin/bash

#this code is for redownloading .vcf.gz files from dnanexus and re-tabixing them

module load bcftools/1.3.1
module load htslib/1.3.2
module load vcftools/0.1.13

# input the list of sample/folder names to redownload them from dnanexus and re-tabix them after download. 
#this will change directories too
for i in NBM-R-10 NBM-R-1; do
	cd /home/ahorning/DNAseq/DNAseq_fq_VCF/"$i"
	dx cd /DNAseq_fq_VCF/${i}
	dx download -f ${i}.vcf.gz
	tabix ${i}.vcf.gz
done



# JP34B_USD16080377_H7VLWALXX_L1.vcf.gz
# JP61B_USD16080378_H7VLWALXX_L1.vcf.gz
# JP6B_USD16080376_H7VLWALXX_L1.vcf.gz
# JPAdenoCa_USD16080379_H7VLWALXX_L1.vcf.gz


# cd /home/ahorning/DNAseq/DNAseq_fq_VCF/JP34B
# dx cd /DNAseq_fq_VCF/JP34B
# dx download -f JP34B_USD16080377_H7VLWALXX_L1.vcf.gz
# tabix JP34B_USD16080377_H7VLWALXX_L1.vcf.gz

# cd /home/ahorning/DNAseq/DNAseq_fq_VCF/JP61B
# dx cd /DNAseq_fq_VCF/JP61B
# dx download -f JP61B_USD16080378_H7VLWALXX_L1.vcf.gz
# tabix JP61B_USD16080378_H7VLWALXX_L1.vcf.gz

# cd /home/ahorning/DNAseq/DNAseq_fq_VCF/JP6B
# dx cd /DNAseq_fq_VCF/JP6B
# dx download -f JP6B_USD16080376_H7VLWALXX_L1.vcf.gz
# tabix JP6B_USD16080376_H7VLWALXX_L1.vcf.gz


# cd /home/ahorning/DNAseq/DNAseq_fq_VCF/JPAdenoCa/
# dx cd /DNAseq_fq_VCF/JPAdenoCa
# dx download -f JPAdenoCa_USD16080379_H7VLWALXX_L1.vcf.gz
# tabix JPAdenoCa_USD16080379_H7VLWALXX_L1.vcf.gz


