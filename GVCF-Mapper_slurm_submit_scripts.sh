#!/bin/bash
# Paul Billings Ross asked me to remap the VCF files I have with gvcf-mapper-cl.py found https://github.com/StanfordBioinformatics/googva/blob/master/gvcf-mapper-cl.py
# John Hanks provided the initial code which i used for the Annovar vcf annotation.
# Run for a bunch of vcf files.

cd /home/ahorning/DNAseq/DNAseq_fq_VCF

#files=( $(find "$(pwd)" -name "*.g.vcf" | sort) )

files=( $(find "$(pwd)" -name "*.g.vcf" | sort | grep -v "mapped") )


for file in ${files[*]} ; do
  directory=$(dirname $file) #/home/ahorning/DNAseq/DNAseq_fq_VCF/NBM-R-1
  sample=${directory##*/} #NBM-R-1
  pushd $directory > /dev/null 2>&1 
  echo "Creating/submitting job script for ${sample}"
 
  # Note: The heredoc lines start with a TAB, not spaces. <<- ignores leading 
  # tab characters so they can be used for indenting but not get passed in the actual output
  # This helps with readability.
  cat > ${sample}.sbatch <<- EOF 
	#!/bin/bash
	#SBATCH --job-name=${sample}_gvcfmapper
	#SBATCH --nodes=1
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task=1
	#SBATCH --time=2-00:00:00
	#SBATCH --mail-type=BEGIN,FAIL,END
	#SBATCH --mail-user=ahorning@stanford.edu
	#SBATCH --chdir=${directory}
	#SBATCH --mem=256G
	#SBATCH --account=mpsnyder

	date
	echo "This job runs the gvcf mapper python script on all of the DNAseq_fq_VCF g.vcf files. It merges the reference matching blocks in a gVCF and output a new gVCF."
	echo "I ran on host: \$(hostname -s)"

	echo "SLURM Environment is:"
	env | grep "SLURM" | sort

	echo "My limits are:"
	ulimit -a

	../../scripts/gvcf-mapper-cl.py -g $file -o "${sample}.mapped.g.vcf"

	# table_annovar.pl $file /srv/gsfs0/software/annovar/annovar-20170717/humandb/ \
	#                    -buildver hg38 \
	#                    -out "${sample}_myanno" \
	#                    -remove \
	#                    -protocol refGene,clinvar_20170905,cosmic70,dbnsfp33a \
	#                    -operation gx,f,f,f \
	#                    -nastring .  \
	#                    -polish  \
	#                    -xref /srv/gsfs0/software/annovar/annovar-20170717/example/gene_xref.txt \
	#                    -vcfinput

	date
	EOF
  sbatch ${sample}.sbatch
  popd
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




