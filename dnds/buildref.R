This vignette shows how to perform driver discovery and selection analyses with dNdScv
in cancer sequencing data. By default, dNdScv assumes human data from the GRCh37/hg19 assembly, 
but the latest version of the package provides a function (buildref) to generate the necessary 
reference file to run dNdScv on others species or assembly (see the relevant tutorial). 
Although designed for cancer genomic studies, dNdScv can be also used to quantify selection 
in other resequencing studies, such as SNP analyses, mutation accumulation studies in bacteria 
or for the discovery of mutations causing developmental disorders using data from human trios.



library(dndscv)
#build new reference
#first download the hg38 data from ensembl: http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html

# cds_table = read.delim("~/aarons_FAP_github_repository/reference/cds_biomart/cds_file_biomart.txt",header = T,sep = "\t")
# chr.index = which(cds_table$Chromosome.scaffold.name %in% c(1:22,"X","Y"))
# cds_table_chr = cds_table[chr.index,]
# table(cds_table_chr$Chromosome.scaffold.name)
# write.table(x = cds_table_chr,file = "~/aarons_FAP_github_repository/reference/cds_biomart/cds_file_biomart_chr1-22-X-Y-only.txt",
#             append = F,quote = F,sep = "\t",row.names = F)

path_cds_table = "~/aarons_FAP_github_repository/reference/cds_biomart/cds_file_biomart_chr1-22-X-Y-only.txt"
path_genome_fasta = "~/aarons_FAP_github_repository/reference/fasta/Homo_sapiens_assembly38.fasta"
buildref(cdsfile=path_cds_table, genomefile=path_genome_fasta, outfile = "aarons_FAP_github_repository/dnds/hg38_refcds.rda", excludechrs="MT")


