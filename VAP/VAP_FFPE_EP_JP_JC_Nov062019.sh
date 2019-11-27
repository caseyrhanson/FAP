VAP Pipeline for A001 and A002 with Blood
WGS


####### Step 1 for preparing for mapping 11/20/19: 

SAMPLE=$1 #EPA7  EPA8  JCA1  JCA3  JPA24  JPA26
ROOT=$2 # /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019
SOMATICINFO=$3 #/home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/SomaticInfoTable_EP-JP-JC_WGS.txt
READPOOL=$4 # /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/hwftp.novogene.com/H202SC19101824/Rawdata
FASTQ1=$5
FASTQ2=$6
CONFIG="/labs/ccurtis2/ruping/tools/seqare/config_hg38_wgs.tsv" #change all the paths
READLEN=151
THREADS=4

Deep WGS #redo with 
  for sample in EPA7 EPA8 JCA1 JCA3 JPA24 JPA26; \
   do echo sbatch aaron_step1_prepare_mapping_slurm.sh $sample \
    /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019 \
    /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/SomaticInfoTable_EP-JP-JC_WGS.txt \
    /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/hwftp.novogene.com/H202SC19101824/Rawdata \
    ${sample}/${sample}_1.fq.gz \
    ${sample}/${sample}_2.fq.gz
    done #| bash


for sample in EPA7 EPA8 JCA1 JCA3 JPA24 JPA26; \
   do 
    bwa mem -r 1.2 -t 4 -R '@RG     ID:1    SM:JCA1 PL:ILLUMINA' /labs/ccurtis2/ruping/annotation/hg38/BWAIndex/genome.fa \
 '<zcat /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/hwftp.novogene.com/H202SC19101824/Rawdata/JCA1/JCA1_1.fq.gz' \
 '<zcat /home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/hwftp.novogene.com/H202SC19101824/Rawdata/JCA1/JCA1_2.fq.gz' \
 | samtools view -bS - >/home/ahorning/DNAseq/FFPE_WGS/Novogene_DeepWGS_October2019/JCA1/02_MAPPING/JCA1.bam
    done #| bash
 


# Shallow WGS
# for sample in A001C004 A001C005 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213 B001A001 B001A006 B001A101 B001A201 B001A301 B001A401 B001A406 B001A501; \
#    do echo sbatch aaron_step1_prepare_mapping_slurm.sh $sample \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002 \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002/SomaticInfoTable_A001-A002_WGS.txt \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/WGS_30X \
#     ${sample}_1.fq.gz \
#     ${sample}_2.fq.gz
#     done #| bash

# #redo on 7/21/19
# Shallow WGS
# for sample in B001A201 B001A301; \
#    do echo sbatch aaron_step1_prepare_mapping_slurm.sh $sample \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002 \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002/SomaticInfoTable_A001-A002_WGS.txt \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/WGS_30X \
#     ${sample}_1.fq.gz \
#     ${sample}_2.fq.gz
#     done #| bash

# Shallow WGS Blood
# for sample in A001_blood A002_blood; \
#    do echo sbatch aaron_step1_prepare_mapping_slurm.sh $sample \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002 \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002/SomaticInfoTable_A001-A002_WGS.txt \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Blood_WGS_FAP-patients \
#     ${sample}_1.fq.gz \
#     ${sample}_2.fq.gz
#     done #| bash


# for sample in EP_blood JP_blood; \
#    do echo sbatch aaron_step1_prepare_mapping_slurm.sh $sample \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002 \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002/SomaticInfoTable_A001-A002_WGS.txt \
#     /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Blood_WGS_FAP-patients \
#     ${sample}_1.fq.gz \
#     ${sample}_2.fq.gz
#     done #| bash

# SAMPLE=$1
# ROOT=$2
# SOMATICINFO=$3
# READPOOL=$4
# FASTQ1=$5
# FASTQ2=$6
# CONFIG="/srv/gsfs0/projects/curtis/ruping/tools/seqare/config_hg38_wgs.tsv"
# READLEN=151
# THREADS=4


make a log directory with the specific ${sample}.log.mutectScan files within it.

####### Step 2 for mutect scanning on 7/25/19: 
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213; \
 do echo sbatch aaron_step2_mutectScan.sh $sample \
  ; done #| bash


##### Step 3 has 2 parts
#Run step3 per chromosome first: 7/28/19
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213;\
  do for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    echo sbatch aaron_step3_samtoolsScan_perchr.sh $sample $i; done;
  done #| bash

#rerun step 3 on 7/29/19:
for sample in A002_blood A001C222;\
  do for i in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
    echo sbatch aaron_step3_samtoolsScan_perchr.sh $sample $i; done;
  done #| bash

#Run step3 merge chromosones after runing samtools scan separately first: 8/1/19
#look closely at the results. This step leaves the per chromosome files so you may need to delete those.
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213;\
  do echo sbatch aaron_step3_samtoolsScan_mergechr.sh $sample ; done #| bash

#Run step2p5 because strelka is needed: 7/30/19
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213; \
 do echo sbatch aaron_step2p5_Strelka.sh $sample \
  ; done #| bash


######## Ran step4 on 8/1/18:
sbatch aaron_step4_mergeRawVariant_allsample.sh


######## Run step5 on 8/1/19: recheck variants.
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213;\
 do for i in mutect samtools; do
  echo sbatch aaron_step5_recheck_slurm.sh $sample $i ; done
 done  #| bash

######## Ran on 8/2/19: variant classificaiton. Includes changes made to the algorithm to account for repeat regions and indels
### Reran on 8/5/19 because of an issue with perl version. replaced module commands with Rupings working copy.
### Reran on 8/6/19
### need to make sure that when it is rerun, to delete certain mutect files which werent completed AND the root/titan/ directory
sbatch aaron_step6_varClassification_allsample_slurm.sh

####### Step 7 (QC wig) ran on 8/2/19
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213;\
 do echo sbatch aaron_step7_QCwig_slurm.sh $sample ; done #| bash

####### Step 8 for Titan ran on 8/6/19: (if it doesn't work make sure the .bed target files dont have "chr" in chromosome line) 
for sample in A001_blood A001C004 A001C005 A001C007 A001C021 A001C102 A001C107 A001C122 A001C210 A001C214 A001C219 A001C222 A002_blood A002C019 A002C020 A002C022 A002C023 A002C102 A002C103 A002C107 A002C122 A002C202 A002C207 A002C209 A002C213;\
 do echo sbatch aaron_step8_CNcall_slurm.sh $sample ; done #| bash


####### Step 2 for mutect scanning on 9/10/19:
####### rerun for A001C222 because the .mod.vcf wasn't created for some reason: 
for sample in A001C222; \
 do echo sbatch aaron_step2_mutectScan.sh $sample \
  ; done #| bash


#### step 2p6 changes the symbolic links creatd by strelka (only perform if necessary)
#rerun step4 on 1/16/19: i changed the symbolic links so needed to rerun the merging. If the symbolic links already worked for merging, skip this step
sbatch aaron_step4_mergeRawVariant_allsample_slurm_WXS.sh 



