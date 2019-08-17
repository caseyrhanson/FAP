rm(list=ls())
library(GenomicRanges)

setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")

source("forAaron_source.R")

#source("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/forAaron_source.R")
#setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP")

setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/WES_Nov2018/")
#  loading data and prepare object "d"  #############################################################

d = read.delim("mutect.snv.res.filtered.classified.founds.nopara.somatic.table", header=T)
d$link = NULL
d$mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad = NULL
dg = d[which(!is.na(d$germline)),]
d = d[which(is.na(d$germline)),]
d = d[which(d$chr != "X" & d$chr != "Y"),]
d = data.frame(d, dron=0)
colnames(d) = gsub("\\.", "-", colnames(d))

#"EP_2normals" 

#  generate mutation Table for EP #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^EP", samplesall)]
snormal = "EP_2normals"
samples = setdiff(samples, snormal)
EP = getSampMutMulti(samples, snormal, d)

save.image(file = "getSampMutMulti_output_EP.RData")
load("getSampMutMulti_output_EP.RData")

#  adjust for MAF and get CCF  ########################################################################

# calculte CCF for nondiploid cases: Titan normal proportion versus Ploidy plot
samples.nondiploid = c("EP025",  "EP027",  "EP088",  "EP088B", "JP002B", "JP032B", "JP045",  "JP065") #samples from WES which need CCF adjustment 
samples.nondiploid = samples.nondiploid[str_detect(samples.nondiploid,"EP")] #only adjusts samples which need it
EP = adjust.ccf.titan.multi(EP, samples.nondiploid, 0.05, titanPath="./titan/", correctColname=TRUE)

# calculate CCF for Diploid cases: Titan normal proportion versus Ploidy plots
samples=samples[!(samples %in% samples.nondiploid)]
EP = adjust.ccf.titan.multi.fordiploidsamples(sampAB = EP,samples = samples)
save.image(file = "EP_adjust.ccf.titan.multi.RData")
write.table(x = EP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.EP.txt",
          append = F,sep="\t", col.names = T, row.names = F)

# write.table(x = EP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.preCCF.EP.txt",
#             append = F,sep="\t", col.names = T, row.names = F)


#"JP_2normals"
#  generate mutation Table for JP #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^JP", samplesall)]
snormal = "JP_2normals"
samples = setdiff(samples, snormal)
JP = getSampMutMulti(samples, snormal, d)

save.image(file = "getSampMutMulti_output_JP.RData")
load("getSampMutMulti_output_JP.RData")

#  adjust for MAF and get CCF  ########################################################################
# calculte CCF for nondiploid cases: Titan normal proportion versus Ploidy plot
samples.nondiploid = c("EP025",  "EP027",  "EP088",  "EP088B", "JP002B", "JP032B", "JP045",  "JP065") #samples from WES which need CCF adjustment 
samples.nondiploid = samples.nondiploid[str_detect(samples.nondiploid,"JP")] #only adjusts samples which need it
JP = adjust.ccf.titan.multi(JP, samples.nondiploid, 0.05, titanPath="./titan/", correctColname=TRUE)

# calculate CCF for Diploid cases: Titan normal proportion versus Ploidy plots
samples=samples[!(samples %in% samples.nondiploid)]
JP = adjust.ccf.titan.multi.fordiploidsamples(sampAB = JP,samples = samples)
save.image(file = "JP_adjust.ccf.titan.multi.RData")
write.table(x = JP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.JP.txt",
            append = F,sep="\t", col.names = T, row.names = F)

####### Info about output of table#######
# `samplename`+`ccf` is the ccf-and-cna-adjusted Variant Allele Frequency. 
#       ccf cutoff of 0.1 may be good threshold.
# `ccfsd` is the standard deviation of ccf estimate (confidence in the ccf)
### parameters estimated in Titan:
# `pu` is purity 
# `pa` is cellular prevalence
# `nt` is total copy number
# `nb` is minor copy number
# `seg` is segment id
# `AGP` is abnormal genotype percentage
# `mafa` is VAF
# `ccf` is CCF. Should plot mostly with this.(CCF is VAF multiplied by two until 1)
