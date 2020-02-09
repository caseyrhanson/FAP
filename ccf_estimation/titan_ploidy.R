rm(list=ls())
library(stringr);library(RColorBrewer)
setwd("/Users/ahorning/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002")

#Make a list of the samples which are in the VAP output file
load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
ccf.data = list(ptA001,ptA002,ptF,ptG)
ccf.sets = length(ccf.data)

# gather samples names for A001
pt="A001"

gather.sample.names = function(pt){
  snps = read.table(paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"),header = T,sep = "\t",stringsAsFactors = F)[1:10,]
  
  pt.cols = colnames(snps)[which(grepl(pattern = pt,colnames(snps)))]
  
  pt.samples.d = pt.cols[which(str_ends(string = pt.cols,pattern = "d"))]
  pt.samples = str_sub(pt.samples.d,end = -2)
  return(pt.samples)
}

pts = c("A001","A002")
#Make Ploidy and Copy Number Changes plots
for (pt in pts) {
  sample.names = gather.sample.names(pt)
  # i=8
  #Remove blood sample
  sample.names.no.blood = sample.names[!grepl(pattern = "blood",x = sample.names)]

  #Create empty vector for each of the sample names. Name each spot
  ploidys = c(rep(x = "NULL",length(sample.names.no.blood)))
  names(ploidys) = sample.names.no.blood
  ploidys
  
  #Fill in empty vector with ploidy values from each sample
  for (i in 1:length(sample.names.no.blood)) {
    segment.file = list.files("titan/segments/")[grepl(pattern = sample.names.no.blood[i],x = list.files("titan/segments/"))]
    x = read.table(paste0("titan/segments/",segment.file),header = T,sep = "\t")
    ploidys[i] = as.numeric(unique(x$ploidy))
  }
  
  #Create Plot for Ploidy values
  par(mar = c(5,4,4,2))
  layout(1)
  layout.show()
  svg(filename = paste0("Ploidy_",pt,".svg"))
  plot(ploidys,main = paste0("Ploidy Values for ",pt," Samples"),xaxt = "n",ylab = "Ploidys",ylim = c(0,4))
  text(x = c(1:length(ploidys)),y = as.numeric(ploidys),labels = names(ploidys),pos = 3,cex = 0.6,srt = 45)
  text(x = c(1:length(ploidys)),y = as.numeric(ploidys),labels = ploidys,pos = 1,cex = 0.6,srt = 45)
  dev.off()
  
  
  ##### Create barplots for the Types of Copy Number Changes in each sample
  
  # Determine How many of each type of CNV there are in each sample.
  # Descriptions of states: homozygous deletion (HOMD), hemizygous deletion LOH (DLOH), copy neutral LOH (NLOH),
  # diploid heterozygous (HET), amplified LOH (ALOH), gain/duplication of 1 allele (GAIN), 
  # allele-specific copy number ampli- fication (ASCNA), balanced copy number amplification (BCNA), 
  # unbalanced copy number amplification (UBCNA). State -1 represents the outlier state (OUT).
  
  
  #Create an empty matrix with columns as sample names and rows and columns as types of copy number changes
  m = matrix(0,nrow = 10,ncol = length(sample.names.no.blood))
  colnames(m) = sample.names.no.blood
  rownames(m) = c("HOMD", "DLOH","NLOH","HET", "ALOH", "GAIN", "ASCNA", "BCNA","UBCNA","OUT")
  m = data.frame(m,stringsAsFactors = F)
  
  #Fill in the table for each sample with the abundance of types of copy number changes observed in each sample
  for (i in 1:length(sample.names.no.blood)) {
    segment.file = list.files("titan/segments/")[grepl(pattern = sample.names.no.blood[i],x = list.files("titan/segments/"))]
    segments = read.table(paste0("titan/segments/",segment.file),header = T,sep = "\t")
    cnvs = table(segments$LOHcall)
    for (l in 1:length(cnvs)) {
      r = names(cnvs)[l]
      c = names(ploidys[i])
      m[r,c] = as.numeric(cnvs[l])
    }
  }
  
  m = as.matrix(m)
  svg(filename = paste0("CopyNumberChange-Types_",pt,".svg"),height = 8, width = 10)
  layout.matrix <- matrix(c(1,2), nrow = 2, ncol = 1)
  layout.matrix
  layout(mat = layout.matrix,
         heights = c(3, 1), # Heights of the two rows
         widths = c(2, 2)) # Widths of the two columns
  layout.show(2)
  par(mar = c(5,4,4,2))
  colors = brewer.pal(n = length(rownames(m)),name = "Set3")
  x = barplot(height = m,
              col = colors,
              cex.axis = 0.6,
              xaxt = "n",
              main = paste0("Copy Number Changes in Samples from ",pt),
              axes = T)
  legend("topright",legend = rownames(m),fill = colors)
  x
  text(cex=0.6, x=x-.1, y=-80,labels = colnames(m), xpd = T,srt = 70)
  
  par(mar = c(1,1,1,1))
  plot(1,1,type="n",axes = F,ann = F)
  text(1,1,cex = 0.9,paste0( "Descriptions of states: homozygous deletion (HOMD), hemizygous deletion LOH (DLOH), copy neutral LOH (NLOH),
                             diploid heterozygous (HET), amplified LOH (ALOH), gain/duplication of 1 allele (GAIN),
                             allele-specific copy number ampli- fication (ASCNA), balanced copy number amplification (BCNA),
                             unbalanced copy number amplification (UBCNA). State -1 represents the outlier state (OUT)."))
  dev.off()
}

