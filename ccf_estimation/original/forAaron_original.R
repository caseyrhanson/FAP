library(GenomicRanges)
# source("~/DNAseq/scripts/RupingPipe/forAaron_source.R")
source("~/aarons_FAP_github_repository/ccf_estimation/original/forAaron_source_original.R")


setwd("~/DNAseq/root/")
#  loading data and prepare object "d"  #############################################################

d = read.delim("mutect.snv.res.filtered.classified.founds.nopara.somatic.table", header=T)
d$link = NULL
d$mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad = NULL
dg = d[which(!is.na(d$germline)),]
d = d[which(is.na(d$germline)),]
d = d[which(d$chr != "X" & d$chr != "Y"),]
d = data.frame(d, dron=0)
colnames(d) = gsub("\\.", "-", colnames(d))


#  generate mutation Table for EP #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^EP", samplesall)]
snormal = "EP-107-NL"
samples = setdiff(samples, snormal)
EP = getSampMutMulti(samples = samples,normal = snormal, d)

save.image(file = "getSampMutMulti_output_EP.RData")
load("getSampMutMulti_output_EP.RData")

#  adjust for MAF and get CCF  ########################################################################

EP = adjust.ccf.titan.multi(EP, samples, 0.05, titanPath="./titan/", correctColname=TRUE)
save.image(file = "EP_adjust.ccf.titan.multi.RData")
write.table(x = EP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.adjustedMAF.EP",
          append = F,sep="\t", col.names = T, row.names = F)



#  generate mutation Table for JP #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^JP", samplesall)]
snormal = "JP_Dec_NL-1"
samples = setdiff(samples, snormal)
JP = getSampMutMulti(samples, snormal, d)

save.image(file = "getSampMutMulti_output_JP.RData")
load("getSampMutMulti_output_JP.RData")

#  adjust for MAF and get CCF  ########################################################################

JP = adjust.ccf.titan.multi(JP, samples, 0.05, titanPath="./titan/", correctColname=TRUE)
save.image(file = "JP_adjust.ccf.titan.multi.RData")
write.table(x = JP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.adjustedMAF.JP",
            append = F,sep="\t", col.names = T, row.names = F)
