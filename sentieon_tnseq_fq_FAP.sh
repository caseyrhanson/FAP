#!/bin/bash
set -e #exit script if any command exits with a nonzero exit status
set -u # exits script if any variable is unset
set -o pipefail
set -x #prints every line run

# Usage: <sentieon_tnseq_fq_FAP.sh>


#This code will iteratively run sentieon TNseq fq on FAP patient samples.
#I will need to change the code for every patient. Add a list of 
#"Normal" and "Tumor" samples for each patient above the "for/loop".


#select FAP project
dx select project-F7vFb2j0ZgxJ3BPB18fkxKQG #takes it back to the root directory

#Samples to analyze. Change this section every time
# EP
normal_sample="EP-107-NL" #only need to change this once per patient
# declare -a tumor_samples=("EP_Dec_NL" "EP-11" "EP-26" "EP-31" "EP-57" "EP-6" "EP-74B" "EP-84" "EP-88B" "EP-AdenoCa") #creates array of tumor/polyp samples
declare -a tumor_samples=("EP-AdenoCa") #creates array of tumor/polyp samples


# ES
# normal_sample="ES_Dec_NL" #only need to change this once per patient
# declare -a tumor_samples=("ES-3" "ES24") #creates array of tumor/polyp samples

# JP   first set
# normal_sample="JP_Dec_NL-1" #only need to change this once per patient
# declare -a tumor_samples=("JP3" "JP31"  "JP38"  "JP63" "JP69" "JP9") #creates array of tumor/polyp samples

#JP    second set. had to adjust some things how the fastq.gz s were input because there are multiple per sample
# normal_sample="JP_Dec_NL-1" #only need to change this once per patient
# declare -a tumor_samples=("JP34B" "JP61B" "JP6B"  "JPAdenoCa") #creates array of tumor/polyp samples

#TNseq algorithm inputs choices
genome_index="Reference Genome Files:/H. Sapiens - GRCh38/GRCh38.no_alt_analysis_set.bwa-index.tar.gz"
genome_fasta="Reference Genome Files:/H. Sapiens - GRCh38/GRCh38.no_alt_analysis_set.fa.gz"
gatk_resource="Apps Data:/sentieon_resources/GRCh38.gatk.resource.bundle.tar.gz"
tn_algo="TNscope"

dx mkdir "${tn_algo}" -p #produces no error message if this folder already exists

  for i in "${tumor_samples[@]}"; do #echo "$i"
    ## Make output folder on DNAnexus project
    dx mkdir "${tn_algo}/${i}_v_${normal_sample}/" -p

    tumor_read1=/"$i"/"$i"_R1.fastq.gz #file location in dnanexus tumor_read1=/EP_Dec_NL/EP_Dec_NL_R1.fastq.gz
    tumor_read2=/"$i"/"$i"_R2.fastq.gz

    normal_read1=/"$normal_sample"/"$normal_sample"_R1.fastq.gz
    normal_read2=/"$normal_sample"/"$normal_sample"_R2.fastq.gz

    job_name="${i}_v_${normal_sample}_TNseqFQtoVCF"

    dx run sentieon_tnseq_fq  -itumor_reads_fastqgzs=$tumor_read1 \
      -itumor_reads2_fastqgzs=$tumor_read2 \
      -inormal_reads_fastqgzs=$normal_read1 \
      -inormal_reads2_fastqgzs=$normal_read2 \
      -igenomeindex_targz="$genome_index" \
      -igenome_fastagz="$genome_fasta" \
      -igatk_resource_bundle="$gatk_resource" \
      -itn_algo="$tn_algo" \
      -ioutput_refined_bam_to_output="recaled.bam" \
      -imark_or_remove_duplicate="mark_duplicate" \
      -imark_as_secondary=true \
      -iignore_decoy=true \
      -itumor_sample="$i" \
      -inormal_sample="$normal_sample" \
      --folder "${tn_algo}/${i}_v_${normal_sample}/" \
      --brief -y --name "$job_name"
  done

# To add FAP project number: project-F7vFb2j0ZgxJ3BPB18fkxKQG
# for help with sentieon_tnseq_fq algorithm
# usage: dx run sentieon_tnseq_fq [-iINPUT_NAME=VALUE ...]

# App: Sentieon TNseq FASTQ to VCF

# Runs Sentieon's TNseq FASTQ to VCF pipeline

# See the app page for more information:
#   https://platform.dnanexus.com/app/sentieon_tnseq_fq

