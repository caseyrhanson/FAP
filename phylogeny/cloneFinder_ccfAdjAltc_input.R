# Example input
# SNVID	Wild	Mut	R2:ref	R2:alt	R3:ref	R3:alt	R4:ref	R4:alt	
# S1	T	A	70	30	70	30	50	50	
# S2	C	A	70	30	70	30	50	50
# S3	A	T	70	30	70	30	50	50
# S4	A	T	70	30	70	30	50	50
# S5	C	T	70	30	70	30	40  60CPCOLS <- c("#1f78b4", "#33a02c", "#e31a1c", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")

rm(list=ls())

library(stringr);library(GenomicRanges);library(IRanges);library(dplyr);library(tidyr)

setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
source("~/aarons_FAP_github_repository/vapReduce_cover,ccf,somatic.R")
############ Prep for CloneFinder

# FAP CCF Datasets
# load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")

pt = "A001"

sampAB = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
data = vapReduce_cover_ccf_somatic(data = sampAB,patient = pt,mafOrccf = "ccf")
snvs = data[[1]]
samples = data[[4]]
normal = data[[3]]

#Gather the d columns, not the "id" column though and not the Normal sample either
d = snvs[,grepl(pattern = "d$",x = colnames(snvs)) & !grepl(pattern = "id",colnames(snvs)) & !grepl(normal,colnames(snvs))]

# Make Altc Columns
ccf.altc = snvs[,grep(pattern = "_ccf_adjusted_altc",x = colnames(snvs))]
colnames(ccf.altc) = str_replace_all(string = colnames(ccf.altc),
                                     pattern = "_ccf_adjusted_altc",
                                     replacement = ":alt")

#Make Ref Columns
ref = d-ccf.altc #Total depth subtracted by the CCF-adjusted AltC values
colnames(ref) = str_replace_all(string = colnames(ref),pattern = "d",replacement = "\\:ref")

colnames(sampAB)
snvid = sampAB%>%unite(col= "SNVID" ,c("chr","pos","id","ref","alt"),sep=":")%>%select(SNVID)
wild.mut = sampAB%>%select(Wild = "ref",Mut = "alt")

input = data.frame(cbind(snvid,wild.mut,ref,ccf.altc))
colnames(input) = str_replace_all(string = colnames(input),pattern = "\\.",replacement = "\\:")

head(input)

#Write input without quotes
write.table(x = input,file = paste0("clonefinder/",pt,"_CCF_adjusted_Altc.txt"),
            append = F,sep = "\t",row.names = F,quote = F)

