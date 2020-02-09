rm(list=ls())
library(GenomicRanges);library(stringr)

# setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")

source("~/aarons_FAP_github_repository/ccf_estimation/forAaron_source.R")

#setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP")

setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
#  loading data and prepare object "d"  #############################################################

#read in file labelled "mutect.snv.res.filtered.classified.founds.nopara.somatic.table". Not the "simplified" file.
d = read.delim("mutect.snv.res.filtered.classified.founds.nopara.somatic.table", header=T,sep = "\t")
d$link = NULL
d$mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad = NULL
dg = d[which(!is.na(d$germline)),]
d = d[which(is.na(d$germline)),]
d = d[which(d$chr != "X" & d$chr != "Y"),]
d = data.frame(d, dron=0)
colnames(d) = gsub("\\.", "-", colnames(d))

#A001_redo if needed
#  generate mutation Table for A001 #############################################################
    samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
    samples = samplesall[grepl("^A001", samplesall)]
    snormal = "A001_blood"
    samples = setdiff(samples, snormal)
    # A001 = getSampMutMulti(samples, snormal, d)
    # 
    # save.image(file = "getSampMutMulti_output_A001.RData")
load("getSampMutMulti_output_A001.RData")

#  adjust for MAF and get CCF  ########################################################################

# calculte CCF for nondiploid cases: Titan normal proportion versus Ploidy plot
# samples.nondiploid = c("EP025",  "EP027",  "EP088",  "EP088B", "JP002B", "JP032B", "JP045",  "JP065") #samples from WES which need CCF adjustment 
# samples.nondiploid = samples.nondiploid[str_detect(samples.nondiploid,"EP")] #only adjusts samples which need it


#onlly include samples which have columns before the "geneName" column. 
#These samples likely have cancer cells for the following analysis
index.genename.col = which(grepl("geneName",colnames(d)))
left.cols = colnames(d)[-c(index.genename.col:length(colnames(d)))]
left.cols.samples = left.cols[which(grepl("maf",left.cols))]
left.cols.samples.nomaf = str_replace_all(left.cols.samples,pattern = "maf",replacement = "")

samples = samples[samples %in% left.cols.samples.nomaf]

A001 = adjust.ccf.titan.multi(sampAB = A001,samples = samples,
                              t =  0.05, titanPath = "./titan/segments/",
                              correctColname = TRUE)

save.image(file = "adjust.ccf.titan.multi.A001.RData")
write.table(x = A001,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt",append = F,sep = "\t")

# # calculate CCF for Diploid cases: Titan normal proportion versus Ploidy plots
# samples=samples[!(samples %in% samples.nondiploid)]
# EP = adjust.ccf.titan.multi.fordiploidsamples(sampAB = EP,samples = samples)
# save.image(file = "EP_adjust.ccf.titan.multi.RData")
# write.table(x = EP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.EP.txt",
#           append = F,sep="\t", col.names = T, row.names = F)
# 
# # write.table(x = EP,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.preCCF.EP.txt",
# #             append = F,sep="\t", col.names = T, row.names = F)


#"A002"
#  generate mutation Table for A002 #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^A002", samplesall)]
snormal = "A002_blood"
samples = setdiff(samples, snormal)
A002 = getSampMutMulti(samples, snormal, d)

save.image(file = "getSampMutMulti_output_A002.RData")
load("getSampMutMulti_output_A002.RData")

A002 = adjust.ccf.titan.multi(sampAB = A002,samples = samples,
                              t =  0.05, titanPath = "./titan/segments/",
                              correctColname = TRUE)


# #  adjust for MAF and get CCF  ########################################################################
# # calculte CCF for nondiploid cases: Titan normal proportion versus Ploidy plot
# samples.nondiploid = c("EP025",  "EP027",  "EP088",  "EP088B", "A002002B", "A002032B", "A002045",  "A002065") #samples from WES which need CCF adjustment 
# samples.nondiploid = samples.nondiploid[str_detect(samples.nondiploid,"A002")] #only adjusts samples which need it
# A002 = adjust.ccf.titan.multi(A002, samples.nondiploid, 0.05, titanPath="./titan/", correctColname=TRUE)
# 
# # calculate CCF for Diploid cases: Titan normal proportion versus Ploidy plots
# samples=samples[!(samples %in% samples.nondiploid)]
# A002 = adjust.ccf.titan.multi.fordiploidsamples(sampAB = A002,samples = samples)
save.image(file = "adjust.ccf.titan.multi.A002.RData")
write.table(x = A002,file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A002.txt",
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
# `AGP` is anueploid genotype percentage
# `mafa` is VAF
# `ccf` is CCF. Should plot mostly with this.(CCF is VAF multiplied by two until 1)


# Testing filtration steps in getSampMutMulti() function with EP
#read in file labelled "mutect.snv.res.filtered.classified.founds.nopara.somatic.table". Not the "simplified" file.
path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table"
d = read.delim(path, header=T,sep = "\t");rm(path)
d$link = NULL
d$mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad = NULL
dg = d[which(!is.na(d$germline)),]
d = d[which(is.na(d$germline)),]
d = d[which(d$chr != "X" & d$chr != "Y"),]
d = data.frame(d, dron=0)
colnames(d) = gsub("\\.", "-", colnames(d))

#A001_redo if needed
#  generate mutation Table for A001 #############################################################
samplesall = as.character(gsub("maf", "", colnames(d)[which(grepl("maf$", colnames(d)))]))
samples = samplesall[grepl("^EP", samplesall)]
snormal = "EP-107-NL"
samples = setdiff(samples, snormal)
EP = getSampMutMulti(samples = samples,normal = snormal, d)

