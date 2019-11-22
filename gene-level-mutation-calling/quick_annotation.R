library(dplyr);library(tidyr)
setwd("~/aarons_FAP_github_repository/gene-level-mutation-calling/")


clinvar = read.table(file = "clinvar_result_APC.txt",header = T,sep = "\t")

tumor = read.table(file = "A001_blood_APC.sorted.ir.br.rmDup.md.bam.vcf",sep = "\t")
colnames(tumor) = c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",	"A001C007")
tumor$POS = as.numeric(as.character(tumor$POS))
str(tumor)

clinvar = separate(clinvar,GRCh38Location,into = c("POS","stop"),sep = " - ")
clinvar$POS = as.numeric(as.character(clinvar$POS))

str(clinvar)


left_join(tumor,clinvar,"POS")
