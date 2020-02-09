#https://bioconductor.org/packages/3.7/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf

rm(list=ls())

library(dplyr);library(tidyr);library(stringr)


# setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")
setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")

load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")

source("~/aarons_FAP_github_repository/vapReduce_cover,ccf,somatic.R")
pt = "JP"
maf = vapReduce_cover_ccf_somatic(data = G,
                                  patient = pt,mafOrccf = "ccf")[[1]]
maf = maf%>%select(-"id")
maf$Change = paste0(maf$ref,">",maf$alt)
x = maf
maf = x

#Add "chr" to chromosome names
maf$chr = as.character(maf$chr)
index.chr = which(!grepl(pattern = "chr",maf$chr))
maf$chr[index.chr] = paste0("chr",maf$chr[index.chr])


#### Create the phred input table
vars = which(grepl(pattern = "d$",x = colnames(maf))) #choose columns of sample names
phred = maf%>%
  select(Chromosome = "chr", Position = "pos", Change, Gene = "geneName",vars)%>%
  select(-contains("blood"))

#remove the "d" in sample names
colnames(phred) = str_remove_all(string = colnames(phred),pattern = "d$")
head(phred)


### Create the mutant read data input table
vars = which(grepl(pattern = "_ccf_adjusted_altc$",x = colnames(maf))) #choose columns of sample names
mutreads = maf%>%
  select(Chromosome = "chr", Position = "pos", Change, Gene = "geneName",vars)%>%
  select(-contains("blood"))

#remove the "d" in sample names
colnames(mutreads) = str_remove_all(string = colnames(mutreads),pattern = "_ccf_adjusted_altc$")
head(mutreads)


# Make sure phred and mutreads have the same columns
# phred should only use columns found in mutreads because although every sample has Depth information,
# Not every sample has CCF information (which can be converted to CCF-adjusted Altc Values)
phred = phred[,colnames(phred) %in% colnames(mutreads)]

#Treeomics - https://github.com/johannesreiter/treeomics
#phredcoverage.txt  
# Chromosome	Position	Change	Gene	Pam01N3	LiM 1	LiM 2	NoM 1	NoM 2	PT 1	PT 2
# # Chromosome	Position	Change	Gene	Pam01N3	Pam01PT4C 	Pam01PT1C	Pam01PT3C 	Pam01PT2C 	Pam01PT7 	Pam01PT8 
# chr12p13.31	9220824	A>-	A2M	44	242	126	340	642	268	117
# chr16q22.1	70299546	G>A	AARS	383	34	28	353	431	336	139

#mutatant_reads.txt
# Chromosome	Position	Change	Gene	Pam03N3	LiM 1	LiM 2	LiM 3	LiM 4	LiM 5	LuM 1	LuM 2	LuM 3	PT 10	PT 11
# chr22q11.22	23243367	T>C	abParts	4	7	7	0	0	1	14	0	2	3	13
# chr5p15.32	5232549	G>A	ADAMTS16	0	183	221	17	56	13	49	6	2	66	242
# chr5p15.31	7706895	G>A	ADCY2	0	257	285	14	169	11	8	4	4	118	252

#Save the files

if (pt %in% c("A001","A002")) {
  path = "/Users/ahorning/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002"
  write.table(mutreads,file = paste0(path,"/treeomics/",pt,"_WGS_mutant_reads_CCF_Adjusted.txt"),
              sep = "\t", append = FALSE, row.names = FALSE)
  write.table(phred,file = paste0(path/"/treeomics/",pt,"_WGS_phred_CCF_Adjusted.txt"),
              sep = "\t", append = FALSE, row.names = FALSE)
} else if (pt %in% c("EP","JP")) {
  path = "/Users/ahorning/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP"
  write.table(mutreads,file = paste0(path,"/treeomics_WGS_output/",pt,"_WGS_mutant_reads_CCF_Adjusted.txt"),
              sep = "\t", append = FALSE, row.names = FALSE)
  write.table(phred,file = paste0(path,"/treeomics_WGS_output/",pt,"_WGS_phred_CCF_Adjusted.txt"),
              sep = "\t", append = FALSE, row.names = FALSE)
} else { NULL }

