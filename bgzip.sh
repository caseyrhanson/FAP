
module load htslib
z 
#BGZF-compressed variant calling data

for sample in EP-AdenoCa EP_Dec_NL ES24 ES-3 ES_Dec_NL JC-AdenoCa JC-ASC JC-SIG JP3 JP31 JP34B JP38 JP61B JP63 JP69 JP6B JP9 JPAdenoCa JP_Dec_NL-1 NBM-R-1 NBM-R-10; do
	echo bgzip -d /home/ahorning/DNAseq/DNAseq_fq_VCF/"${sample}"/"${sample}".g.vcf.gz | bash
done




# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-107-NL/EP-107-NL.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-11/EP-11.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-26/EP-26.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-31/EP-31.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-57/EP-57.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-6/EP-6.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-74B/EP-74B.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-84/EP-84.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-88B/EP-88B.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-AdenoCa/EP-AdenoCa.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/EP_Dec_NL/EP_Dec_NL.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/ES24/ES24.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/ES-3/ES-3.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/ES_Dec_NL/ES_Dec_NL.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-AdenoCa/JC-AdenoCa.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-ASC/JC-ASC.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-SIG/JC-SIG.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP31/JP31.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP34B/JP34B_USD16080377_H7VLWALXX_L1.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP38/JP38.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP3/JP3.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP61B/JP61B_USD16080378_H7VLWALXX_L1.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP63/JP63.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP69/JP69.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP6B/JP6B_USD16080376_H7VLWALXX_L1.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP9/JP9.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JPAdenoCa/JPAdenoCa_USD16080379_H7VLWALXX_L1.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/JP_Dec_NL-1/JP_Dec_NL-1.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/NBM-R-10/NBM-R-10.g.vcf.gz
# /home/ahorning/DNAseq/DNAseq_fq_VCF/NBM-R-1/NBM-R-1.g.vcf.gz
