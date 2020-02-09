rm(list=ls())

library(dplyr);library(tidyr);library(stringr)
# detach(name = "package:VariantAnnotation",unload = T)

setwd("~/Bulk_A001_A002/")
source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")
source("~/DNAseq_WGS/scripts/find.sample.and.tumor.names.R")
source("~/aarons_FAP_github_repository/phylogeny/MaxParsimony_inspiredbyWillCrossTutorial_input.R")

#read in data
data = read.delim(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt",header = T,sep = "\t")
colnames(data)
sample.names = find.sample.and.tumor.names(vap = data,patient = "A001",maforccf = "ccf")$sample.names

#read in cancer drivers list
drivers = read.delim(file = "~/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv",
           header = T,sep = ",",skip = 3)
drivers = unique(drivers$Gene)

#filter data by cancer drivers for working with a smaller set of data
data.drivers = data[data$geneName %in% "KRAS",]


data.drivers$somatic
data.drivers
x = maxpars_input(pt = "A001",vap = data.drivers)
x

sample.names = find.sample.and.tumor.names(vap = data.drivers,
                                           patient = "A001",
                                           maforccf = "ccf")$sample.names


altc
refc
d
ccf


head(data.drivers)
for (i in sample.names) {
  data.drivers[,c("chr","pos","id","ref","alt","geneName","somatic",paste0(i,))]
}
colnames(data.drivers)









    