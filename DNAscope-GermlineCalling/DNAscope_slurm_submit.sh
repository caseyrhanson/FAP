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
#SBATCH --partition=nih_s10
####  consider adding this if it keeps messing up --partition=nih_s10


# usage: stil need to run these
# for sample in EP025 EP088B EP095; do
# 	echo sbatch DNAscope_slurm_submit.sh $sample
# done

# usage: stil need to run these
# for sample in EP088; do
# 	echo sbatch DNAscope_slurm_submit.sh $sample
# done



sample=$1
echo $1

date
echo "This job runs the sentieon's DNAscope algorithm on all of the VAP .bam files. This run is for ${sample}"
echo "I ran on host: \$(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

cd /home/ahorning/DNAseq/WES_batch2_novogene/DNAscope
module load sentieon

#with the given argument
# sentieon driver -i ~/DNAseq/WES_batch2_novogene/root/${sample}/02_MAPPING/${sample}.sorted.ir.br.rmDup.md.bam \
# -r ~/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
# -t ${SLURM_CPUS_PER_TASK} \
# --algo DNAscope \
# -d /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Reference/dbsnp_146.hg38.vcf.gz \
# --call_conf 10 \
# --given /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/WES_batch2_novogene/DNAscope/EP_JP_WES_Nov2018_hg38.vcf.gz \
# ~/DNAseq/WES_batch2_novogene/DNAscope/mergegiven/${sample}.mergegiven.DNAscope.dbSNP.vcf.gz 2>> ~/DNAseq/WES_batch2_novogene/DNAscope/mergegiven/log.mergegiven.DNAscope.${sample}.txt

#without the given argument
sentieon driver -i ~/DNAseq/WES_batch2_novogene/root/${sample}/02_MAPPING/${sample}.sorted.ir.br.rmDup.md.bam \
-r ~/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
-t ${SLURM_CPUS_PER_TASK} \
--algo DNAscope \
-d /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/Reference/dbsnp_146.hg38.vcf.gz \
--call_conf 10 \
~/DNAseq/WES_batch2_novogene/DNAscope/mergegiven/${sample}.mergegiven.DNAscope.dbSNP.vcf.gz 2>> ~/DNAseq/WES_batch2_novogene/DNAscope/mergegiven/log.mergegiven.DNAscope.${sample}.txt

date


