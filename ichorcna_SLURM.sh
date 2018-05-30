#!/bin/bash
#SBATCH --job-name=aaron_ichorcna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --workdir=/home/ahorning/DNAseq/scripts/RupingPipe
#SBATCH --account=mpsnyder


# Running IchorCNA to take over for Chris's cfdnapipeline. For some reason it didnt' finish. 
# Using mostly Chris's parameters.
# Added "normal" 0.1 because these are tumor samples basically.

#change working directories to get the sample input ncbi wig files.
cd /home/ahorning/curtis_aaron/data_miseq_DNAseq_NBM15-55-40/analysis/cfdnapipeline-output_hg19/
sampleIDs=( $(ls mapped/ | cut -f 1 -d . | sort | uniq) ) #NBM_15-1_S11_L001_R1_001 NBM_15-2_S84_L001_R1_001 NBM_15-3_S24_L001_R1_001 NBM_40-2_S58_L001_R1_001 NBM_40-3_S71_L001_R1_001 NBM_40-6_S45_L001_R1_001 NBM_55-2_S19_L001_R1_001 NBM_55-3_S32_L001_R1_001 NBM_55-4_S6_L001_R1_001

#run the ichorCNA pipeline which was installed according to instructions online

module load R/3.4.1

for sample in ${sampleIDs[*]}; do
	#500Kb or 500kb
	Rscript /home/ahorning/software/ichorCNA/scripts/runIchorCNA.R --id "${sample}" \
	  --WIG "/home/ahorning/curtis_aaron/data_miseq_DNAseq_NBM15-55-40/analysis/cfdnapipeline-output_hg19/mapped/${sample}.dedup.cov.500Kb.ncbi.wig" --ploidy "c(2,3)" --normal "c(0.1,0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
	  --gcWig "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/gc_hg19_500kb.wig" \
	  --mapWig "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/map_hg19_500kb.wig" \
	  --centromere /home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
	  --normalPanel "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds" \
	  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
	  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
	  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir "/home/ahorning/curtis_aaron/data_miseq_DNAseq_NBM15-55-40/analysis/cfdnapipeline-output_hg19/ichor-cna/${sample}.dedup.filt.bam" > "${sample}_500Kb_ichor_stdout.txt" 2> "${sample}_500Kb_ichor_stderror.txt"

	#1Mb or 1000kb
	Rscript /home/ahorning/software/ichorCNA/scripts/runIchorCNA.R --id "${sample}" \
	  --WIG "/home/ahorning/curtis_aaron/data_miseq_DNAseq_NBM15-55-40/analysis/cfdnapipeline-output_hg19/mapped/${sample}.dedup.cov.1Mb.ncbi.wig" --ploidy "c(2,3)" --normal "c(0.1,0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
	  --gcWig "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/gc_hg19_1000kb.wig" \
	  --mapWig "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/map_hg19_1000kb.wig" \
	  --centromere /home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
	  --normalPanel "/home/ahorning/software/R/x86_64-pc-linux-gnu-library/3.4/ichorCNA/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds" \
	  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
	  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
	  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir "/home/ahorning/curtis_aaron/data_miseq_DNAseq_NBM15-55-40/analysis/cfdnapipeline-output_hg19/ichor-cna/${sample}.dedup.filt.bam" > "${sample}_1Mb_ichor_stdout.txt" 2> "${sample}_1Mb_ichor_stderror.txt"	
done
