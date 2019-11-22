#!/bin/bash
#SBATCH --job-name=DNAscope_Aaron
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --chdir=/home/ahorning/DNAseq/scripts
#SBATCH --mem=1000G
#SBATCH --account=mpsnyder
#####SBATCH --partition=nih_s10
####  consider adding this if it keeps messing up --partition=nih_s10


sample=$1
echo $1
input="/home/ahorning/DNAseq/Bulk_WGS-WES_May2019/root_WGS-A001-A002/"$sample"/02_MAPPING/${sample}.sorted.ir.br.rmDup.md.bam"
out='/home/ahorning/DNAseq/Bulk_WGS-WES_May2019/DNAscope_B001'
reference='/labs/mpsnyder/aaron/FAP/DNAseq/Reference/Homo_sapiens_assembly38.fasta'
dbsnp='/labs/mpsnyder/aaron/FAP/DNAseq/Reference/dbsnp_146.hg38.vcf.gz'

given='/home/ahorning/DNAseq/Bulk_WGS-WES_May2019/DNAscope_B001/B001A_WGS_Oct2019_hg38.vcf.qual100dp10.gz'

date
echo "This job runs the sentieon's DNAscope algorithm on all of the VAP .bam files. This run is for ${sample}"
echo "I ran on host: \$(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

cd $out
module load sentieon

#without the given argument
#sentieon driver -i $input \
#-r $reference \
#-t ${SLURM_CPUS_PER_TASK} \
#--algo DNAscope \
#-d $dbsnp \
#--call_conf 10 \
#${out}/${sample}.DNAscope.dbSNP.vcf.gz 2>> ${out}/log.DNAscope.${sample}.txt

#with the given argument
sentieon driver -i $input \
-r $reference \
-t ${SLURM_CPUS_PER_TASK} \
--algo DNAscope \
-d $dbsnp \
--call_conf 10 \
--given $given \
${out}/mergegiven/${sample}.mergegiven.DNAscope.dbSNP.vcf.gz 2>> ${out}/mergegiven/log.mergegiven.DNAscope.${sample}.txt


date





