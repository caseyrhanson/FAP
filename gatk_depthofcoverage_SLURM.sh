#!/bin/bash
#SBATCH --job-name=DepthofCoverage_GATK
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=300G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --workdir=/home/ahorning/DNAseq/scripts
#SBATCH --account=mpsnyder


#The below link can tell us more about the DepthOfCoverage analysis
#https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php

module load gatk/3.7

date
echo "This job ran GATKs DepthOfCoverage on all of the root directory .bam files"
echo "I ran on host: $(hostname -s)"

echo "SGE Environment is:"
env | grep "SGE" | sort

echo "My limits are:"
ulimit -a

cd /home/ahorning/DNAseq/root/depthofcoverage_test
echo $(pwd)
java -Xmx16384M -jar /srv/gsfs0/software/gatk/gatk-3.7/GenomeAnalysisTK.jar \
	-T DepthOfCoverage \
	-R /srv/gsfs0/shared_data/RefGenomes/GATK_Resource_Bundle/hg38/Homo_sapiens_assembly38.fasta \
	-o depthtest_wobases \
	--omitDepthOutputAtEachBase \
	-I /home/ahorning/DNAseq/root/input_bams.list >> /home/ahorning/DNAseq/root/depthofcoverage_test/log.txt

date


# T name of thing
# R reference fasta
# o file_name_base
# I input_bams_list


# ERROR MESSAGE: An error occurred because you did not provide enough memory to run this program. 
# You can use the -Xmx argument (before the -jar argument) to adjust the maximum heap size provided to Java. Note that this is a JVM argument, not a GATK argument.

# depth-cover is several orders of magnitude faster than other tools. Execution time and memory 
# consumption depend on the size of the BAM file. In a modern desktop computer, it can process 
# 15-30 million reads per minute.

# Recommended memory allocation is 15% - 50% of the size of the BAM file - i.e.: if your BAM is 10Gb,
# execute like: 

     # java -Xmx2g -XX:+UseParallelGC -jar depth-cover.jar ARGUMENTS

# The JVM option -XX:+UseParallelGC is not mandatory, but usually it is a good idea (particularly, if 
# memory is a scarce resource).

# If the BAM file is indexed and bigger than 4 Gb, and there is enough memory, depth-cover will read
# the chromosomes in parallel, performing up to 50% faster. 

# The parallel reader makes intensive use of CPU and RAM. 
# If you use a shared computer, you might want to disable it with the --ignore-index option.

