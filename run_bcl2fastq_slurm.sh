#!/bin/bash
#SBATCH --job-name=aaron_HiSeq_bcl2fastq
#SBATCH --nodes=1 #good
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 #adjust to suit my needs
#SBATCH --time=4-00:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=ALL ##BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --workdir=/home/ahorning/ATACseq/data_miseq_ATAC_NBM15-55-40/
#SBATCH --account=mpsnyder

#module load bcl2fastq2/2.20.0 #this is the newer version though
module load bcl2fastq2/2.17.1.14

set -e

#usage: ./home/ahorning/ATACseq/scripts/run_bcl2fastq_slurm.sh /home/ahorning/curtis_aaron/data_miseq_omniATACseq_NBM15-55-40/ SampleSheet_052518.csv

data_dir=$1 #The first argument should be the path to the data_miseq_ATAC_NBM15-55-40 file or one like it
runfolders="${data_dir}Files"
sample_sheet=$2

echo "The .csv sample sheet should be found here: ${runfolders}"

#runfolders="/home/ahorning/ATACseq/data_miseq_ATAC_NBM15-55-40/Files" #

for runfolder in $runfolders; do

    bcl2fastq \
        --ignore-missing-bcls \
        --ignore-missing-filter \
        --ignore-missing-positions \
        --ignore-missing-controls \
        --barcode-mismatches 0 \
        --loading-threads 1 \
        --demultiplexing-threads 2 \
        --processing-threads 2 \
        --writing-threads 1 \
        -R "${runfolder}" \
        --sample-sheet "${runfolder}/${sample_sheet}" \
        -o "${runfolder}/Analysis/Fastqs" \
        -i "${runfolder}/Data/Intensities/BaseCalls" #2>&1 ${runfolder}_bcl2fastq.log.txt  #-i [ --input-dir ] arg (=<runfolder-dir>/Data/Intensities/BaseCalls/)

#     fastqdir=${runfolder}/Analysis/Fastqs
# 
#     gsutil -m cp -r ${fastqdir} gs://gbsc-gcp-lab-curtis-atac-projects/fastqs/${runfolder}
#     gsutil cp ${runfolder}_bcl2fastq.log.txt \
#         gs://gbsc-gcp-lab-curtis-atac-projects/fastqs/bcl2fastqlogs/${runfolder}/${runfolder}_bcl2fastq.log.txt
# 
#     rm -rf ${runfolder}
#     rm ${runfolder}_bcl2fastq.log.txt

done;



