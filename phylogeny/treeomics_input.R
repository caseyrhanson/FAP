#https://bioconductor.org/packages/3.7/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf

rm(list=ls())

library(dplyr);library(tidyr);library(stringr)


setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")


############ Prep for Treeomics
pts=c("JP","EP")
pt="EP"
for (pt in pts) {
  sampAB<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.preCCF.q100.",pt,".txt"), header = T, sep = "\t")
  
  #select ref and alt columns and melt the table. Modify sample names too.
  snvs<-select(sampAB,
               Chromosome = chr,
               Position = pos,
               ref, alt,
               geneName, geneLoc, functionalClass,
              ends_with("d"), ends_with("altc"),-id,-CADD_phred,-Polyphen2_HVAR_pred,-ends_with("ccfSD"))%>%
    unite(col = "Change",c("ref","alt"),sep = ">")
  #Make new geneName column to specify exonic versus non-exonic mutations
  for (i in 1:nrow(snvs)) {
    if (snvs[i,"geneLoc"]=="exonic") {
      snvs$Gene[i] = paste0(snvs[i,"geneName"])
    } else {
      snvs$Gene[i] = paste0(snvs[i,"geneName"],"[",snvs[i,"geneLoc"],"]")
    }
  }
  
  #Make mutation reads table  
  mutreads = select(snvs,Chromosome,	Position,	Change,	Gene, ends_with("altc"))
  colnames(mutreads)[5:length(colnames(mutreads))] = str_sub(string = colnames(mutreads)[5:length(colnames(mutreads))],
                                                             start = 1,
                                                             end = -5)
  colnames(mutreads)[5:length(colnames(mutreads))] = str_replace_all(string = colnames(mutreads)[5:length(colnames(mutreads))],
                                                                 pattern = "\\.",replacement = "_")    
  colnames(mutreads)  
  
    #Make phred table
  phred = select(snvs,Chromosome,	Position,	Change,	Gene, ends_with("d"))
  colnames(phred)[5:length(colnames(phred))] = str_sub(string = colnames(phred)[5:length(colnames(phred))],
                                               start = 1,
                                               end = -2)
  colnames(phred)[5:length(colnames(phred))] = str_replace_all(string = colnames(phred)[5:length(colnames(phred))],
                                                           pattern = "\\.",
                                                           replacement = "_")
  phred = phred[colnames(mutreads)]
  
  write.table(phred,file = paste0(pt,"_WES_phredcoverage.txt"),
              sep = "\t", append = F, row.names = F)
  
  write.table(mutreads,file = paste0(pt,"_WES_mutant_reads.txt"),
              sep = "\t", append = F, row.names = F)
}

  
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

##Check VAFs - histogram or violin plot
pts=c("JP","EP")
pt="EP"

sampAB<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt"), header = T, sep = "\t")
  
  #select ref and alt columns and melt the table. Modify sample names too.
  snvs<-select(sampAB,
               Chromosome = chr,
               Position = pos,
               ref, alt,
               geneName, geneLoc, functionalClass,
               ends_with("d"), ends_with("altc"), ends_with("ccf"),-id,-CADD_phred,-Polyphen2_HVAR_pred,-ends_with("ccfSD"),-mergeCCF)
colnames(snvs)

# x = snvs[intersect(which(snvs$JP6Baltc>=7),which(snvs$JP31altc>=7)),]
# snvs$JP6Baltc>=7
# snvs$JP31altc>=7
# 
# b = x[intersect(which(x$JP6Bd>=25),which(x$JP31d>=25)),]
# 
# write.csv(x = b,file = "treeomics_WGS_output/allVariants/JP6BvJP31_othermutslikeKRAS.csv",append = F,sep = "\t")



z = ggplot(snvs,snvs[colnames(snvs)[which(grepl("ccf$",colnames(snvs)))]])
geom_violin(z)

vioplot(snvs[,which(grepl("ccf$",colnames(snvs)))])
colnames(snvs)
do.call(vioplot,c(snvs[,which(grepl("ccf$",colnames(snvs)))],list(names=colnames(which(grepl("ccf$",colnames(snvs)))))))


library(reshape2)
mdat <- melt(t(as.matrix(snvs[,grepl("ccf$",colnames(snvs))])))
ggplot(mdat,aes(x=factor(Var1),y=value))+geom_violin()+
  labs(x="variable")+ylim(0,0.25)
hist(mdat$value)
