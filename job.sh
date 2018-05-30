
#!/bin/sh
#
SAMPLE=$1
ROOT=$2
SOMATICINFO=$3
READPOOL=$4
FASTQ1=$5
FASTQ2=$6
CONFIG="/srv/gsfs0/projects/curtis/ruping/tools/seqare/config_hg38_wgs.tsv"
READLEN=151
THREADS=4

#this step does preparation of directories and mapping
echo "perl /srv/gsfs0/projects/curtis/ruping/tools/seqare/DTrace.pl --qcOFF --configure $CONFIG --runlevel 1-2 --sampleName $SAMPLE --skipTask recalMD --readpool $READPOOL --FASTQ1 $FASTQ1 --FASTQ2 $FASTQ2 --seqType paired-end,WGS,ignore --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/$SAMPLE.run.log"
