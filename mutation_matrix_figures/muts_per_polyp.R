rm(list=ls())
library(plot.matrix)
#Gather data from individual samples

pt = "A001"
source("~/aarons_FAP_github_repository/phylogeny/MaxParsimony_inspiredbyWillCrossTutorial_input.R")
data.binary = maxpars_input(pt = pt)

data = read.delim(paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"),
                  header = T)

#total muts
colSums(data.binary)
barplot(as.matrix(data.binary),xaxt = "n")
title(paste0(pt,": Total number of SNVs per Tissue"))
text(cex=1, x=x-.25, y=-1000, colnames(data.binary), xpd=TRUE, srt=45)
text(cex=1, x=x-.25, y=-2000, colSums(data.binary), xpd=TRUE, srt=45)

#exonic muts
data.binary.exonic = data.binary[data$geneLoc=="exonic",]
barplot(as.matrix(data.binary.exonic),xaxt = "n")
title(paste0(pt,": Total number of Exonic SNVs per Tissue"))
text(cex=1, x=x-.25, y=-10, colnames(data.binary), xpd=TRUE, srt=45)
text(cex=1, x=x-.25, y=-20, colSums(data.binary.exonic), xpd=TRUE, srt=45)


#C->T Mutations per polyp
data.binary.CtoT = data.binary[data$ref=="C" & data$alt=="T",]
barplot(as.matrix(data.binary.CtoT),xaxt = "n")
title(paste0(pt,": Total number of C->T SNVs per Tissue"))
text(cex=1, x=x-.25, y=-200, colnames(data.binary), xpd=TRUE, srt=45)
text(cex=1, x=x-.25, y=-300, colSums(data.binary.CtoT), xpd=TRUE, srt=45)


#Shared Mutations between samples
samples = head(colnames(data.binary),-1)
x = matrix(data = NA,nrow = length(samples),ncol = length(samples))
rownames(x) = samples;colnames(x) = samples

row = 1
col = 1
for (row in 1:length(samples)) {
  for (col in 1:length(samples)) {
    samples[row]
    
    x[row,col] = 
  }
}



data.binary = data.frame(data.binary,
                         functionalClass = data$functionalClass,
                         geneLoc = data$geneLoc,
                         geneName = data$geneName)
head(data.binary)
sample.index = grep(pattern = pt,x = colnames(data.binary))



