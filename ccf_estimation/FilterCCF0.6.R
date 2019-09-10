#This code is for filtering out non-clonal mutations. 
#We will only keep clonal mutations around which have CCF >= 0.6
rm(list=ls())

library("GenVisR");library(dplyr);library(tidyr);library(stringr);library(reshape2);library(ggplot2)
setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")


pts=c("EP","JP")
i=1
for (i in 1:length(pts)) {
  file<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",
                                 pts[i],".txt"),
                   header = T,sep = "\t")
  
  
  
  ccf.samps.bool = grepl(pattern = "ccf$",x = colnames(file)) #which columns end in ccf 
  
  
  #find max CCF across all columns which have CCF at end of name
  max = do.call(pmax, file[ccf.samps.bool])
  
  #Determine if the max value of the row across all CCF columns is >= 0.6 
  max.bool = max >= 0.6
  
  #select only columns which have a clonal mutation (one of the mutations is CCF >= 0.6)
  small.file = file[max.bool,];length(max.bool);nrow(file)
  
  write.table(small.file,
              file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.clonal0.6.",
                            pts[i],".txt")
              ,append = F,
              sep = "\t")
}




