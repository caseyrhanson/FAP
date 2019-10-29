# Read in libraries and set seed

rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR") #already installed on Aaron's macbookpro
# Load the GenVisR package
library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2);library(GGally);library(vioplot)
set.seed(426)


#################### Make CCF v Depth scatter plots for every samples

#read in the sample table with the sample names and its associated tissue types
setwd("~/aarons_FAP_github_repository/")
source("recent_Annotation.R")
recent = recent_Annotation("no")
recent = data.frame(recent,stringsAsFactors = F)


##### For A001 and A002 graphs
# pts = c("A001","A002"); setwd("~/Bulk_A001_A002/");pt = "A002"
# file<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")

#### For EP and JP
pts = c("EP","JP"); setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/");pt = "JP"
file<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt"), header = T, sep = "\t")

#Create a new folder for the Depth v CCF figures
ifelse(file.exists("Depth_v_CCF"),NA,dir.create("Depth_v_CCF"))

#de-identify the EP and JP values
if (pt=="EP") {
  colnames(file) = str_replace(colnames(file),paste0(pt,"."),"F_")
} else if (pt=="JP") {
  colnames(file) = str_replace(colnames(file),pt,"G_")
}


#Rename samples to include Stage and Size (mm)
#find sample names
samples = str_remove(string = colnames(file)[str_detect(string = colnames(file),pattern = "ccf$")],pattern = "ccf")
#remove "-" from names in recent Sample Names column
recent$SampleName = str_remove_all(string = recent$SampleName,pattern = "-")

#use the VAP names and replace the "." with "_"
EP.JP.index = which(grepl(pattern = "EP|JP",x = recent$SampleName))
recent$SampleName[EP.JP.index] = recent$VAP_Names[EP.JP.index]
recent$SampleName = str_replace_all(string = recent$SampleName,pattern = "\\.",replacement = "_")

#Change EP and JP to F and G in the recent table
recent$SampleName = str_replace(string = recent$SampleName,pattern = "EP",replacement = "F")
recent$SampleName = str_replace(string = recent$SampleName,pattern = "JP",replacement = "G_")
recent$SampleName = str_replace_all(string = recent$SampleName,pattern = "__",replacement = "_")

for (sample in samples) {
  sample.index = recent$SampleName %in% sample
  stage = recent$Stage..Polyp..Normal..AdCa.[recent$SampleName == sample]
  stage = stage[!is.na(stage)]
  size = recent$Size..mm.[recent$SampleName == sample]
  size = size[!is.na(size)]
  recent$SampleName[recent$SampleName == sample] = paste0(sample,"_",stage,"_",size)
}




polyp1 = "G_3"
for (polyp1 in samples) {
  file.sample = cbind(file[,1:5],file[,str_detect(string = colnames(file),pattern = polyp1)])
  ccf.index = which(str_detect(string = colnames(file.sample),pattern = paste0(polyp1,"ccf")))[1]
  file.sample = file.sample[file.sample[,ccf.index]>0,]
  
  
  polyp1ccf = file.sample[,colnames(file.sample) %in% paste0(polyp1,"ccf")]
  polyp1d = file.sample[,colnames(file.sample) %in% paste0(polyp1,"d")]
  
  svg(filename = paste0("Depth_v_CCF/",polyp1,"_Depth_v_CCF.svg"),)
  
  layout(mat = matrix(data = c(3,1,4,2),nrow = 2),widths = c(2,1),heights = c(1,2))
  layout.show(4)
  
  col_vaf = rep("color",nrow(file.sample))
  altc.index = file.sample[,colnames(file.sample) %in% paste0(polyp1,"altc")]
  col_vaf[altc.index<5] = "pink"
  col_vaf[altc.index<=2] = "grey"
  col_vaf[altc.index>=5] = "green"
  
  par(mar = c(5,4,4,2))
  plot(polyp1ccf, polyp1d, xlab = paste0(polyp1," CCF"), ylab = paste0(polyp1," Depth"),ylim = c(0,100),
       col = "darkgrey",pch = 21, bg = col_vaf)
  
  par(mar = c(5,4,4,2))
  vioplot(polyp1d[polyp1d<=100])#,ylim = c(0,100),ylab = paste0(polyp1," Depth"), xlab = polyp1)
  myjitter = jitter(x = rep(1.5,length(polyp1d)),amount = 0.2)
  # points(myjitter,polyp1d, col = "black", bg = rgb(red = .109,green = .109,blue = .109,alpha = 0.4),pch = 21)
  points(myjitter,polyp1d, col = "black", bg = col_vaf,pch = 21, cex = 0.8)
  rug(x = polyp1d,col = "purple",side = 2)
  
  par(mar = c(1,4,1,2))
  vioplot(polyp1ccf,ylab = polyp1,xlab = paste0(polyp1," CCF"), horizontal = T)
  myjitter = jitter(x = rep(1.5,length(polyp1ccf)),amount = 0.1)
  points(polyp1ccf,myjitter,col = "black", bg = col_vaf,pch = 21, cex = 0.8,xpd = NA)
  
  par(mar = c(1,4,1,2))
  plot(1,1,type = "n",ann = F, axes = F,frame.plot = F)
  sample.name = recent$SampleName[grep(pattern = polyp1,x = recent$SampleName)]
  text(1,1.2,paste0(sample.name[1], "\n Depth v CCF"),font = 2, cex = 1.5,xpd = NA)
  legend("bottomleft",legend = c("<= 2 atlc","< 5 altc", ">= 5 altc"),pch = 21,
         border = "grey",pt.bg = c("grey","pink","green"),)
  
  dev.off()
}



# https://www.reed.edu/data-at-reed/resources/R/loops_with_ggplot2.html
# for making an iterative for loop for plotting








