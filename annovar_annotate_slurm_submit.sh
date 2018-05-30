#!/bin/bash
# written Initially by Aaron, then rewritten by John Hanks
# Run for a bunch of vcf files.

cd /home/ahorning/DNAseq/DNAseq_fq_VCF

files=( $(find "$(pwd)" -name "*.vcf" | sort) )

for file in ${files[*]} ; do
  directory=$(dirname $file)
  sample=${directory##*/} 
  pushd $directory > /dev/null 2>&1 
 
  cat << EOF | sbatch # --hold
#!/bin/bash
#SBATCH --job-name=${sample}_myanno
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=ahorning@stanford.edu
#SBATCH --chdir=${directory}
#SBATCH --mem=32G

date
echo "This job runs annovar on 4 databases (refGene,clinvar_20170905,cosmic70,dbnsfp33a) on all of the DNAseq_fq_VCF .vcf files"
echo "I ran on host: \$(hostname -s)"

echo "SLURM Environment is:"
env | grep "SLURM" | sort

echo "My limits are:"
ulimit -a

module load annovar/20170717 

table_annovar.pl $file /srv/gsfs0/software/annovar/annovar-20170717/humandb/ \
                   -buildver hg38 \
                   -out "${sample}_myanno" \
                   -remove \
                   -protocol refGene,clinvar_20170905,cosmic70,dbnsfp33a \
                   -operation gx,f,f,f \
                   -nastring .  \
                   -polish  \
                   -xref /srv/gsfs0/software/annovar/annovar-20170717/example/gene_xref.txt \
                   -vcfinput

date
EOF

done


# Notes from John about the script rewrite.
# That script basically build s a job script for each file and sends that to sbatch. Note the `--hold` which is there so I could check the jobs before they ran. You probably don't want that in there, but it is good for debugging.

# Also these
# ```  directory=$(dirname $file)
#   sample=${directory##*/} ```
# replace some of your finds to build the arrays.

# The `dirname` command is a standard command, there is a corresponding command `basename` which will pull the filename from teh end of a path.

# For the sample name, that is using bash parameter expansion, see https://www.tldp.org/LDP/abs/html/parameter-substitution.html

# In this case I'm doing a replace to chop out the last directory in the path.

# You could also do $(basename $(dirname $directory) ) to get that.

# oh, I left the pushd in there but it's not needed. But in the original approach where cd was used to cahnge to the sample directory, I (purely personal taste) prefer to do that with pushd to the directory then a popd afterwards to create a block which has all the code executed in that directory.

# It's just syntax that helps me when debugging later.




