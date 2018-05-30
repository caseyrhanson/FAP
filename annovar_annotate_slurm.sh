#!/bin/bash
#
# set the name of the job
#SBATCH --job-name=aaron_annovar3
#
# set the maximum memory usage for job
#SBATCH --mem=6G                                 #maybe increase to 8G
#
# set the maximum run time
#SBATCH --time=108:00:00 
#
# send mail when job ends or aborts
#SBATCH --mail-type=BEGIN,FAIL,END
#
# specify an email address
#SBATCH --mail-user=ahorning@stanford.edu
#
# use multiple cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

date
echo "This job runs annovar on 4 databases (refGene,clinvar_20170905,cosmic70,dbnsfp33a) on all of the DNAseq_fq_VCF .vcf files"
echo "I ran on host: $(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

module load annovar/20170717 

cd /home/ahorning/DNAseq/DNAseq_fq_VCF



#file_path=($(find "$(pwd)" -name "*.vcf*" | grep -v ".g.vcf.gz" | grep -v ".tbi" | sort)); echo ${#file_path[@]} # vcf.gz file path and names
files=( $(find "$(pwd)" -name "*.vcf" | sort) )


#file_path_base=($(find "$(pwd)" -name "*.vcf*" | grep -v ".g.vcf.gz" | grep -v ".tbi" | sort | cut -d. -f1)); echo ${#file_path_base[@]} # file path without .vcf
#dir_path=($(find "$(pwd)" -name "*.vcf*" | grep -v ".g.vcf.gz" | grep -v ".tbi" | sort | cut -d "/" -f1-6)); echo ${#dir_path[@]} # path to directory which contains files
#sample=($(find "$(pwd)" -name "*.vcf*" | grep -v ".g.vcf.gz" | grep -v ".tbi" | sort | cut -d "/" -f6)) ; echo ${#sample[@]} # sample name (title of directory)

# consider grep -v "chrUn" for the .vcf file

for file in ${files[*]} ; do
  directory=$(dirname $file)
  sample=${directory##*/} 
  pushd $directory > /dev/null 2>&1 
 
  table_annovar.pl $file /srv/gsfs0/software/annovar/annovar-20170717/humandb/ \
                   -buildver hg38 \
                   -out "${sample[$index]}_myanno" \
                   -remove \
                   -protocol refGene,clinvar_20170905,cosmic70,dbnsfp33a \
                   -operation gx,f,f,f \
                   -nastring .  \
                   -polish  \
                   -xref /srv/gsfs0/software/annovar/annovar-20170717/example/gene_xref.txt \
                   -vcfinput &
done

# Wait until all those sruns we forked above have finished.
wait

#consider including this database in the future dbnsfp33a (filter based)


date

# $ANNOVAR/example/
# $ANNOVAR/humandb/

# Filter the rows
# Func.refGene - filter for exonic. Not ncRNA_exonic
# ExonicFunc.refGene - nonsynonymous and stopgain

# fathmm-MKL_coding_pred - D (likely damaging)
#   fathmm-MKL_coding_pred  FATHMM-MKL  predicting the effects of both coding and non-coding variants using nucleotide-based HMMs Classifier based on multiple kernel learning
#     D: Deleterious; Score >= 0.5: D 
#     T: Tolerated ; Score < 0.5: T


# Columns of Annotations from dbNSFP Database Pediction Algorithm/Conservation Score  Description Method  Categorical Prediction  Author(s)
# SIFT_pred 
# SIFT_score  SIFT  Sort intolerated from tolerated P(An amino acid at a position is tolerated | The most frequentest amino acid being tolerated) D: Deleterious (sift<=0.05);
# T: tolerated (sift>0.05)  Pauline Ng, Fred Hutchinson 
# Cancer Research Center, Seattle, Washington
# Polyphen2_HDIV_pred 
# Polyphen2_HDIV_score  Polyphen v2 Polymorphism phenotyping v2 D: Probably damaging (>=0.957), 
# P: possibly damaging (0.453<=pp2_hdiv<=0.956), 
# B: benign (pp2_hdiv<=0.452) Probablistic Classifier Training sets: HumDiv Havard Medical School/td>

# Polyphen2_HVAR_pred
# Polyphen2_HVAR_score  Polyphen v2 Polymorphism phenotyping v2 Machine learning Training sets: HumVar  
# D: Probably damaging (>=0.957), 
# P: possibly damaging (0.453<=pp2_hdiv<=0.956); 
# B: benign (pp2_hdiv<=0.452) 
# Shamil Sunyaev
# Havard Medical School

# LRT_pred 
# LRT_score LRT Likelihood ratio test LRT of H0: each codon evolves neutrally vs H1: the codon evovles under negative selection D: Deleterious; 
# N: Neutral;
# U: Unknown
# Lower scores are more deleterious Sung Chung, Justin Fay Washington University
# MutationTaster_pred 
# MutationTaster_score  MutationTaster  Bayes Classifier  A: (""disease_causing_automatic""); 
# D: (""disease_causing""); 
# N: (""polymorphism [probably harmless]""); 
# P: (""polymorphism_automatic[known to be harmless]"
# higher values are more deleterious"   Markus Schuelke
# the Charité - Universitätsmedizin Berlin
# MutationAssessor_pred 
# MutationAssessor_score  MutationAssessor  Entropy of multiple sequence alighnment H: high; 
# M: medium; 
# L: low; 
# N: neutral. 
# H/M means functional and L/N means non-functional higher values are more deleterious



# # code for downloading clinvar_20170905 to my local computer.
# ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20170905 humandb/


# Func.refGene, Gene.refGene, GeneDetail.refGene, ExonicFunc.refGene, AAChange.refGene






# ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/

# #./annotate_variation.pl -buildver hg38 -downdb cosmic70 humandb/ downloaded the cancer_gene_census manually
# sftp ahorning@stanford.edu@sftp-cancer.sanger.ac.uk
# pwd is bettyboop

#  cancer_gene_census.csv
#  clinvar_20170905
