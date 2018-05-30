#!/bin/sh
#
# set the name of the job
#$ -N recheck
#
# set the maximum memory usage (per slot)
#$ -l h_vmem=4G
#
# set the maximum run time
#$ -l h_rt=6:00:00
#
# send mail when job ends or aborts
#$ -m ea
#
# specify an email address
#$ -M kpogrebn@stanford.edu
#
# check for errors in the job submission options
#$ -w e
#
# set to local directory
#$ -cwd
#
# use multiple cores
##$ -pe shm 4

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

SAMPLE="TRIO_52_N_merged"
ROOT="/srv/gsfs0/projects/curtis/kpogrebn/TRIO/fastq/root"
SOMATICINFO="/srv/gsfs0/projects/curtis/kpogrebn/TRIO/fastq/somaticInfo"
RECHECK="/srv/gsfs0/projects/curtis/kpogrebn/TRIO/fastq/root/mutect.snv.table.annotated"
CONFIG="/srv/gsfs0/projects/curtis/ruping/tools/seqare/config_nexterarapidcapture.tsv"
READLEN=76
THREADS=1

#this step is for recheck mutation files
perl /srv/gsfs0/projects/curtis/ruping/tools/seqare/DTrace.pl --configure $CONFIG --runTask recheck --sampleName $SAMPLE --recheck $RECHECK --seqType paired-end,WGS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/logs/$SAMPLE.run.log.recheck
