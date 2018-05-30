#!/bin/bash
#
# set the name of the job
#$ -N aaron 
#
# Any following resource requirements are "hard" rather than "soft/suggested"  #### maybe keep
#####$ -hard                                                                       #### maybe keep
#
# set the maximum memory usage (per slot)
#$ -l h_vmem=8G                                 #maybe increase to 8G
#
# set the maximum run time
#$ -l s_rt=108:00:00 #changed to s_rt from h_rt
#
# send mail when job ends or aborts
#$ -m ea
#
# specify an email address
#$ -M ahorning@stanford.edu
#
# check for errors in the job submission options
#$ -w e
#
# set to local directory
#$ -cwd
#
# use multiple cores
#$ -pe shm 8                                    # maybe 12?
#
# Reserve space
#$ -R y

#The below link can tell us more about the DepthOfCoverage analysis
#https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_coverage_DepthOfCoverage.php

module load gatk/3.7

date
echo "This job ran GATKs DepthOfCoverage on all of the TNscope "*markdup.realigned.recalibrated.bam" files"
echo "I ran on host: $(hostname -s)"

echo "SGE Environment is:"
env | grep "SGE" | sort

echo "My limits are:"
ulimit -a

java -Xmx16384M -jar /srv/gsfs0/software/gatk/gatk-3.7/GenomeAnalysisTK.jar \
	-T DepthOfCoverage \
	-R /home/ahorning/DNAseq/Reference/Homo_sapiens_assembly38.fasta \
	-o /srv/gsfs0/projects/snyder/aaron/FAP/DNAseq/TNscope/GATK_DepthofCoverage \
	-I /home/ahorning/DNAseq/TNscope/input_bams.list

date
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

