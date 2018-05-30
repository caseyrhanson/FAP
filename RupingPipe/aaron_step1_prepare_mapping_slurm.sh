#!/bin/bash
#SBATCH --job-name=aaron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --workdir=/home/ahorning/DNAseq/scripts/RupingPipe
#SBATCH --account=mpsnyder


export TOOLP=/srv/gsfs0/projects/curtis/ruping/tools/
export PATH=$PATH:$TOOLP/trick/:$TOOLP/cmake/current/bin/:$TOOLP/gmap/current/bin/:$TOOLP/bedtools/current/bin/:$TOOLP/vcftools/current/bin/:$TOOLP/samtools/current/bin/:$TOOLP/bcftools/current/bin/:$TOOLP/htslib/current/bin/:$TOOLP/bwa/current/:$TOOLP/MACS/bin/:/usr/java/latest/bin/:$TOOLP/bowtie2/current/
export PYTHONPATH=/srv/gsfs0/software/python/python-2.7/lib/python2.7/site-packages/:/home/ruping/.local/lib/python2.7/site-packages/:/srv/gsfs0/projects/curtis/ruping/tools/MACS/lib/python2.7/site-packages/:/srv/gsfs0/projects/curtis/ruping/tools/Cython-0.23.2/
export PERL5LIB=/home/ruping/cpan/export/cpan/cpan/build/Text-NSP-1.31-cU1PW2/blib/lib/:/ifs/home/c2b2/ac_lab/rs3412/tools/vcftools/current/perl/:$PERL5LIB
alias python='/srv/gsfs0/software/python/python-2.7/bin/python'
alias R='/srv/gsfs0/projects/curtis/ruping/tools/R/bin/R'
module load perl-scg
#comment the following line if mapping
#module load java/6u35
module load vcftools/0.1.13

SAMPLE=$1
ROOT=$2
SOMATICINFO=$3
READPOOL=$4
FASTQ1=$5
FASTQ2=$6
CONFIG="/srv/gsfs0/projects/curtis/ruping/tools/seqare/config_hg38_wgs.tsv"
READLEN=151
THREADS=4



date
echo "This job runs step 1 of Rupings VAP pipeline on WGS FASTQ files"
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a


#this step does preparation of directories and mapping
perl /srv/gsfs0/projects/curtis/ruping/tools/seqare/DTrace.pl --maxMem 16g --qcOFF --configure $CONFIG --runlevel 1-2 --sampleName $SAMPLE --skipTask recalMD --readpool $READPOOL --FASTQ1 $FASTQ1 --FASTQ2 $FASTQ2 --seqType paired-end,WGS,ignore --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/$SAMPLE.run.log

date


