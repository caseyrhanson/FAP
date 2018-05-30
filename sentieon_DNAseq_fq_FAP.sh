#!/bin/bash
set -e #exit script if any command exits with a nonzero exit status
set -u # exits script if any variable is unset
set -o pipefail
set -x #prints every line run

# App: Sentieon DNAseq FASTQ to VCF
# Runs Sentieon's DNAseq FASTQ to VCF pipeline  #'
# See the app page for more information:
#   https://platform.dnanexus.com/app/sentieon_dnaseq_fq

#This code will iteratively run sentieon DNAseq_fq on FAP patient samples.
# This algorithm will compare the fastq reads to a reference to determine mutations (as opposed to the MuTect2 algorithm)
# Runs Sentieon's Haplotyper as a recreation of GATK HaplotypeCaller (in gVCF mode) and Sentieon's GVCFtyper as GenotypeGVCFs.

# This script will run DNAseq_fq algorithm on every FAP sample.
# It will run HaplotypeCaller (compare each sample .bam to a reference)
#  and create .vcf files for every sample.
# This script will require modification to run a different set of patient samples


# FAP project number: project-F7vFb2j0ZgxJ3BPB18fkxKQG
# For help with sentieon_dnaseq_fq algorithm
# dx run sentieon_dnaseq_fq -h

# usage: dx run sentieon_dnaseq_fq [-iINPUT_NAME=VALUE ...]


# 1) select FAP project
dx select project-F7vFb2j0ZgxJ3BPB18fkxKQG #takes it back to the root directory of the project

# 2) select specific set of patient samples to run. Change this section every time

#EP
##declare -a samples=("EP-107-NL" "EP_Dec_NL" "EP-11" "EP-26" "EP-31" "EP-57" "EP-6" "EP-74B" "EP-84" "EP-88B" "EP-AdenoCa") #creates array of tumor/polyp samples
declare -a samples=("EP_Dec_NL")

#ES
##declare -a samples=("ES_Dec_NL" "ES-3" "ES24") #creates array of tumor/polyp samples

#JP   #first set with single fastq files.
##declare -a samples=("JP_Dec_NL-1" "JP3" "JP31"  "JP38"  "JP63" "JP69" "JP9") #creates array of tumor/polyp samples

#NBM-R   #first set with single fastq files.
##declare -a samples=("NBM-R-1" "NBM-R-10") #creates array of tumor/polyp samples


#JC   #first set with single fastq files.
#declare -a samples=("JC-AdenoCa" "JC-ASC" "JC-SIG") #creates array of tumor/polyp samples



#JP   #second set with multiple fastqs files. May just need to run this one manually.
#declare -a samples=("JP34B" "JP61B" "JP6B"  "JPAdenoCa") #creates array of tumor/polyp samples

#DNAseq algorithm inputs choices

genome_index="Reference Genome Files:/H. Sapiens - GRCh38/GRCh38.no_alt_analysis_set.bwa-index.tar.gz"
genome_fasta="Reference Genome Files:/H. Sapiens - GRCh38/GRCh38.no_alt_analysis_set.fa.gz"
gatk_resource="Apps Data:/sentieon_resources/GRCh38.gatk.resource.bundle.tar.gz"
algo="DNAseq_fq_VCF" 

  #make output folder in the root directory
  dx mkdir "${algo}" -p # -p produces no error message if this folder already exists

for i in "${samples[@]}"; do
  
 ## Make output folder inside the 'algo' folder 
    dx mkdir "${algo}/${i}/" -p

    sample_read1=/"$i"/"$i"_R1.fastq.gz #file location in dnanexus tumor_read1=/EP_Dec_NL/EP_Dec_NL_R1.fastq.gz
    sample_read2=/"$i"/"$i"_R2.fastq.gz

    job_name="${i}_DNAseqFQtoVCF"


  dx run sentieon_dnaseq_fq -ireads_fastqgzs=$sample_read1 \
  -ireads2_fastqgzs=$sample_read2 \
  -igenomeindex_targz="$genome_index" \
  -igenome_fastagz="$genome_fasta" \
  -igatk_resource_bundle="$gatk_resource" \
  -ioutput_refined_bam_to_output="recaled.bam" \
  -imark_or_remove_duplicate="mark_duplicate" \
  -imark_as_secondary=true \
  -iignore_decoy=true \
  --folder "${algo}/${i}/" \
  --brief -y --name "$job_name"
done

#add these arguments?
# Extra options available for Haplotyper:
# --var_type [SNP/INDEL/BOTH]  --min_base_qual QUALITY --prune_factor FACTOR --emit_conf CONFIDENCE --call_conf CONFIDENCE --pcr_indel_model [HOSTILE/AGGRESSIVE/CONSERVATIVE/NONE] 





# dx run sentieon_dnaseq_fq -h

# usage: dx run sentieon_dnaseq_fq [-iINPUT_NAME=VALUE ...]

# App: Sentieon DNAseq FASTQ to VCF

# Runs Sentieon's DNAseq FASTQ to VCF pipeline  #'

# See the app page for more information:
#   https://platform.dnanexus.com/app/sentieon_dnaseq_fq

# Inputs:
#   Reads: -ireads_fastqgzs=(file) [-ireads_fastqgzs=... [...]]
#         An array of files, in gzipped FASTQ format, with the first read mates to be mapped.

#   Reads (right mates): [-ireads2_fastqgzs=(file) [-ireads2_fastqgzs=... [...]]]
#         (Optional) An array of files, in gzipped FASTQ format, with the second read mates to be
#         mapped.

#   Read group information: [-irg_info_csv=(file)]
#         (Optional) A file, in CSV format, with each row describing RG tag information for each
#         (pair) of read file in the following order: read_filename, RG ID, RG LB, RG PU. This file
#         needs to cover and match all FASTQ filenames of "Reads" (reads_fastqgzs).

#   BWA reference genome index: -igenomeindex_targz=(file)
#         A file, in gzipped tar archive format, with the reference genome sequence already indexed
#         with BWA.

#         Suggestions:
#           project-BQpp3Y804Y0xbyG4GJPQ01xv://file-* (DNAnexus Reference Genomes: AWS US-east)
#           project-Bxv84pQ27z999kj0jbQz0FK9://file-* (DNAnexus Reference Genomes: AWS China-north)
#           project-F3zxk7Q4F30Xp8fG69K1Vppj://file-* (DNAnexus Reference Genomes: AWS Germany)
#           project-F4gXb605fKQyBq5vJBG31KGG://file-* (DNAnexus Reference Genomes: AWS Sydney)
#           project-F0yyz6j9Jz8YpxQV8B8Kk7Zy://file-* (DNAnexus Reference Genomes: Azure West US)

#   Reference genome FASTA file: -igenome_fastagz=(file)
#         A gzipped reference genome file that the reads will be mapped against.

#         Suggestions:
#           project-BQpp3Y804Y0xbyG4GJPQ01xv://file-* (DNAnexus Reference Genomes: AWS US-east)
#           project-Bxv84pQ27z999kj0jbQz0FK9://file-* (DNAnexus Reference Genomes: AWS China-north)
#           project-F3zxk7Q4F30Xp8fG69K1Vppj://file-* (DNAnexus Reference Genomes: AWS Germany)
#           project-F4gXb605fKQyBq5vJBG31KGG://file-* (DNAnexus Reference Genomes: AWS Sydney)
#           project-F0yyz6j9Jz8YpxQV8B8Kk7Zy://file-* (DNAnexus Reference Genomes: Azure West US)

#   GATK resource bundle: [-igatk_resource_bundle=(file)]
#         (Optional) A bundle containing known SNPs, mutations and indels.

#         Suggestions:
#           project-B6JG85Z2J35vb6Z7pQ9Q02j8:/sentieon_resources/file-* (Sentieon GATK Bundles: AWS US-east)
#           project-F29KXx020fXF13JyFQ4qQvB4:/sentieon_resources/file-* (Sentieon GATK Bundles: AWS China-north)
#           project-F3zqGV04fXX5j7566869fjFq:/sentieon_resources/file-* (Sentieon GATK Bundles: AWS Germany)
#           project-F4gYG1850p1JXzjp95PBqzY5:/sentieon_resources/file-* (Sentieon GATK Bundles: AWS Sydney)
#           project-F29g0xQ90fvQf5z1BX6b5106:/sentieon_resources/file-* (Sentieon GATK Bundles: Azure West US)

#   Target coordinates: [-itargets_bed=(file) [-itargets_bed=... [...]]]
#         (Optional) A BED file with target coordinates. If given, processing will be restricted to
#         only inside those coordinates. NOTE: Providing this cooridnates has priority over the option
#         of "Ignore the decoy regions". The variant calling will be within the regions provided.

#  Sample Information
#   Read group platform: [-iread_group_platform=(string, default="ILLUMINA")]
#         Choices: ILLUMINA, SOLID, LS454, HELICOS, PACBIO

#         (Optional) This is a string (without spaces) describing the platform/technology used to
#         produce the reads. The read group platform will appear in the read group information in the
#         BAM file (as RG PL tag).

#   Sample ID: [-isample=(string)]
#         (Optional) A string (without spaces) describing the sample. This Sample ID will appear in
#         the read group information in the BAM file (as RG SM tag) and in the sample information in
#         the VCF file. The output files will be prefixed by this string. If not specified, the Sample
#         ID will be set as the prefix of the FASTQ filename.

#  Output Setting
#   What set of BAM file to output?: [-ioutput_refined_bam_to_output=(string, default="realigned.bam+recal.table")]
#         Choices: realigned.bam+recal.table, recaled.bam

#         What set of BAM file to output? Choosing "realigned.bam+recal.table" will not run ReadWriter
#         algorithm and directly run Haplotyper using the realigned BAM file and the recal table as
#         input, thus reducing the total runtime. Alternatively, choosing "recaled.bam" will run
#         ReadWriter and produce recaled.bam before calling variants.

#  Algorithm Options
#   Mark or remove duplicated reads?: [-imark_or_remove_duplicate=(string, default="mark_duplicate")]
#         Choices: mark_duplicate, remove_duplicate

#         Select "mark_duplicate" will mark duplicated reads and keep all reads in the output BAM
#         file, while select "remove_duplicate" will remove the duplicated reads. Default is mark   #"
#         duplicated reads.

#   Mark shorter split hits as secondary?: [-imark_as_secondary=(boolean, default=true)]
#         (Optional) If long reads are split among multiple locations in the genome (because different
#         parts of the same read align to different locations), select this to mark the shorter ones
#         as secondary alignments. This will ensure better compatibility with Picard and GATK, as some
#         of these tools are not designed for multiple primary split alignments. This will supply the
#         '-M' option to 'bwa mem'.

#   Ignore the decoy regions (i.e. ignore the sponge regions for variant calling)?: [-iignore_decoy=(boolean, default=tru
#         (Optional) Selecting this option will ignore the decoy regions (i.e. chromosomes w/ names as
#         "hs37d5", "chrEBV", or "hs38d1") for variant calling, and will reduce the time cost
#         performing the pipeline. NOTE: if the target coordinates BED file is provided, the applet
#         will then only call variants within the target regions.

#   Compression level to use for writing BAM files. Default is "6".: [-ibam_compression_level=(int, default=6)]
#         Choices: 1, 2, 3, 4, 5, 6, 7, 8, 9

#         (Optional) If set as "1", it will perform the least level of compression (produce BAM file
#         with larger size) and will provide the least time cost generating the BAM file.

#   Extra BWA Options: [-iextra_bwa_options=(string, default="-K 10000000")]
#         (Optional) Extra BWA options.

#   Extra Halotyper Options: [-ihaplotyper_algo_options=(string, default="")]
#         (Optional) Extra options that are algo-level paramters for algo haplotyper.

#   Extra GVCFtyper Options (If wish to match GATK3.7s default call confidence, please add "--call_conf 10 --emit_conf 10" here.): [-igvcftyper_algo_options=(string, default="")]
#         (Optional) Extra options that are algo-level paramters for algo gvcftyper. Note Sentieons
#         default values for call_conf and    emit_conf is 30, while GATK3.7 has its default as 10.

# Outputs:
#   Variants: variants_vcf (file)
#         The genotyped variants in VCF format.

#   Variants index: variants_vcftbi (file)
#         The associated TBI file.

#   Variants: variants_gvcf (file)
#         The called variants in block-gzipped GVCF format.

#   Variants index: variants_gvcftbi (file)
#         The associated TBI file.

#   Deduped mappings: mappings_deduped_bam (file)
#         A coordinate-sorted BAM file with dedupped mappings.

#   Deduped mappings index: mappings_deduped_bai (file)
#         The associated BAM index file with dedupped mappings.

#   Reads stats metrics: metrics (array:file)
#         An array of files containing statistical summaries of the data quality.

#   Realinged mappings: [mappings_realigned_bam (file)]
#         A coordinate-sorted BAM file with realinged mappings.

#   Realinged sorted mappings index: [mappings_realigned_bai (file)]
#         The associated BAM index file with realinged mappings.

#   Recalibrated sorted mappings: [mappings_recaled_bam (file)]
#         A coordinate-sorted BAM file with recalibrated mappings.

#   Recalibrated sorted mappings index: [mappings_recaled_bai (file)]
#         The associated BAM index file with recalibrated mappings.

#   Recalibration table data: [recal_table (file)]
#         The table of the several covariate values, number of observations, number of mismatches, and
#         empirical quality score.

