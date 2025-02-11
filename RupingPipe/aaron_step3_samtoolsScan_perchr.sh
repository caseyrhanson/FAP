#!/bin/bash
#SBATCH --job-name=aaron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=16G
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
module load java/6u35
module load vcftools/0.1.13

#usage example: sbatch aaron_step3_samtoolsScan.sh NBM-R-1 chrY
# runs sample per chromosome. Create a nested for loop to get each sample and each chromosome.

SAMPLE=$1 #file name of samples: EP-107-NL etc.
ROOT="/srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/root"
CONFIG="/srv/gsfs0/projects/curtis/ruping/tools/seqare/config_hg38_wgs.tsv"
SOMATICINFO="/srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/SamplesForSomaticComparison.txt"
READLEN=151
THREADS=1
CHR=$2

######## to retrieve info before run:
date
echo "This job ran $SAMPLE"
echo "I ran on host: $(hostname -s)"

echo "SGE Environment is:"
env | grep "SGE" | sort

echo "My limits are:"
ulimit -a

#this step is for samtools scanning
perl /srv/gsfs0/projects/curtis/ruping/tools/seqare/DTrace.pl --configure $CONFIG --runlevel 4 --sampleName $SAMPLE --chrProcess $CHR --germline samtoolsOnly --seqType paired-end,WGS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/logs/$SAMPLE.run.log.samtoolScan.$CHR

date
