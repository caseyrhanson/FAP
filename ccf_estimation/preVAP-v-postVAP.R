rm(list=ls())
library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2);library(GGally);library(vioplot)


#read in A001 CCF data table
pt = "A001";setwd("~/Bulk_A001_A002/")
path.ccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
data.ccf.A001 = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
#make a new united column: mut.ID
data.ccf.A001 = unite(data.ccf.A001,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
any(duplicated(data.ccf.A001$mut.ID)) #FALSE only if all the mut.ID are distinct


#read in A002 CCF data table
pt = "A002";setwd("~/Bulk_A001_A002/")
path.ccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
data.ccf.A002 = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
#make a new united column: mut.ID
data.ccf.A002 = unite(data.ccf.A002,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
any(duplicated(data.ccf.A002$mut.ID)) #FALSE only if all the mut.ID are distinct


#read in preCCF data table
path.preccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt")
data.preccf = read.table(file = path.preccf,header = T,sep = "\t",stringsAsFactors = F)
#make a new united column: mut.ID
data.preccf = unite(data.preccf,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
any(duplicated(data.preccf$mut.ID)) #FALSE only if all the mut.ID are distinct


A001.muts = data.ccf.A001$mut.ID
names(A001.muts) = rep("A001",length(A001.muts))
A002.muts = data.ccf.A002$mut.ID
names(A002.muts) = rep("A002",length(A002.muts))
preccf.muts = data.preccf$mut.ID
names(preccf.muts) = rep("PreCCF",length(preccf.muts))


# longest = length(preccf.muts)
# A001.muts = c(A001.muts,rep("nothing",longest-length(A001.muts)))
# A002.muts = c(A002.muts,rep("nothing",longest-length(A002.muts)))
# 
# x = data.frame("PreCCF" = preccf.muts,"A001" = A001.muts,"A002" = A002.muts,stringsAsFactors = F)
# write.table(x,"~/Desktop/muts.csv",append = F,sep = ",",col.names = T,row.names = F)

overlap = length(A001.muts[A001.muts %in% A002.muts])

length(A002.muts[A002.muts %in% A001.muts])

length(preccf.muts[preccf.muts %in% A001.muts])

length(preccf.muts[preccf.muts %in% A002.muts])

length(preccf.muts) - length(A001.muts) - length(A002.muts) + overlap
