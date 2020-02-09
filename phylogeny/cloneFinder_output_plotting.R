#Take output from 

rm(list=ls())

library(dplyr);library(tidyr);library(stringr);library(GenomicRanges); library(RColorBrewer)

setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")

############ Prep for Treeomics

pts=c("A001","A002")
pt="A002"
file.type = "ccf"


sampAB = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
sampAB$Wild = sampAB$ref #make new column just renamed
sampAB$Mut =  sampAB$alt #make new column just renamed
#select refc and altd columns and unite the SNVID column
snvs<-select(sampAB,
             Chromosome = chr,
             Position = pos,id,
             Wild, Mut,
             ref, alt,
             ends_with("refc"), ends_with("altc"),-CADD_phred,-Polyphen2_HVAR_pred,-ends_with("ccfSD"),-starts_with("bp"))%>%
  unite(col = "SNVID",c("Chromosome","Position", "id", "ref", "alt" ),sep = ":")


# clonal.pops = read.table("clonefinder/A001_HETsegmentswoNonHets_w_mutIDasSNVIDsnv_CloneFinder.txt",
#                          header = T,sep = "\t",row.names = 1)
clonal.pops = read.table("clonefinder/A002_HETsegmentsonly_20dCut_3altCut_HETsegmentswoNonHetssnv_CloneFinder.txt",
                         header = T,sep = "\t",row.names = 1)
clonal.pops
#Fill in the table for each sample with the abundance of types of copy number changes observed in each sample
m = t(round(as.matrix(clonal.pops),3))


meds <- scale(m, FALSE, rep(1,dim(m)[2])) * 100
meds = rbind(meds,"Undetermined" = 100-colSums(meds))
RColorBrewer::display.brewer.all()
colors = brewer.pal(n = dim(meds)[1],name = "Set3")

# svg(paste0("clonefinder/PercentClonePerSample_",pt,".svg"),width = 10,height = 10)
pdf(paste0("clonefinder/PercentClonePerSample_",pt,".pdf"),width = 10,height = 10)
x = barplot(height = meds,
            col = colors,
            cex.axis = 1,
            ylab = "Percentage of Subclone within Sample",
            xaxt = "n",
            main = paste0("Percentage of Clones per Sample from pt ",pt),
            axes = T,ylim = c(0,100))
legend("topright",legend = rownames(m),fill = colors)
x
text(cex=1, x=x, y=-7,labels = colnames(m), xpd = T,srt = 30)
dev.off()

for (samp in 1:length(colnames(meds))) {
  # svg(filename = paste0(pt,"_subclones_",colnames(meds)[samp],".svg"))
  pdf(paste0("clonefinder/",pt,"_subclones_",colnames(meds)[samp],".pdf"),width = 5,height = 5)
  pie(meds[,samp],labels = rownames(meds),col = colors,main = colnames(meds)[samp])
  dev.off()
}


