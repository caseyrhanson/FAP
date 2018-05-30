#!/usr/bin/bash

## uploads the .mapped.g.vcf files to my google bucket in preparation for their addition to Big Query.

module load google-cloud-sdk
module load htslib


echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-107-NL/EP-107-NL.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-11/EP-11.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-26/EP-26.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-31/EP-31.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-57/EP-57.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-6/EP-6.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-74B/EP-74B.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-84/EP-84.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-88B/EP-88B.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP-AdenoCa/EP-AdenoCa.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/EP_Dec_NL/EP_Dec_NL.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/ES24/ES24.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/ES-3/ES-3.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/ES_Dec_NL/ES_Dec_NL.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-AdenoCa/JC-AdenoCa.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-ASC/JC-ASC.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JC-SIG/JC-SIG.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP31/JP31.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP34B/JP34B.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP38/JP38.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP3/JP3.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP61B/JP61B.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP63/JP63.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP69/JP69.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP6B/JP6B.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP9/JP9.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JPAdenoCa/JPAdenoCa.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/JP_Dec_NL-1/JP_Dec_NL-1.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/NBM-R-10/NBM-R-10.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash
echo gsutil cp /home/ahorning/DNAseq/DNAseq_fq_VCF/NBM-R-1/NBM-R-1.mapped.g.vcf gs://gbsc-gcp-lab-snyder-user-ahorning | bash






