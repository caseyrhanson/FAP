# rm(list=ls())
library(stringr);library(RColorBrewer)
setwd("/Users/ahorning/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002")

#Make a list of the samples which are in the VAP output file
load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
ccf.data = list(ptA001,ptA002,ptF,ptG)
ccf.sets = length(ccf.data)

# gather samples names for A001
samples = c()
for (set in 1:ccf.sets) {
  ccf = ccf.data[[set]]
  sample.index = which(grepl("ccf$",colnames(ccf)))
  samples = c(samples,str_remove_all(colnames(ccf)[sample.index],pattern = "ccf$"))
}

samples

#Make Ploidy and Purity Table
pp.matrix = matrix(data = 0,nrow = length(samples),ncol = 2)
rownames(pp.matrix) = samples
colnames(pp.matrix) = c("Purity", "Ploidy")


for (sample in samples) {
  if (grepl(pattern = "A001|A002",x = sample)) {
    seg.file = list.files(path = "~/Bulk_A001_A002/titan/segments/",pattern = sample)
    segments = read.delim(paste0("~/Bulk_A001_A002/titan/segments/",seg.file))
    
  } else if (grepl("EP|JP",sample)) {
    seg.file = list.files(path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/titan/",
                          pattern = paste0(sample,"_"))
    seg.file = grep(pattern = "segments.txt",seg.file,value = T)
    segments = read.delim(paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/titan/",
                                 seg.file))
  }
  purity = 1- mean(segments$normalproportion)
  ploidy = mean(segments$ploidy)
  
  pp.matrix[sample,"Purity"] = purity
  pp.matrix[sample,"Ploidy"] = ploidy
}
pp.matrix
save(pp.matrix,file = "~/aarons_FAP_github_repository/ccf_estimation/titan_PloidyAndPurity_data.Rdata")

