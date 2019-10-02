rm(list=ls())
library(dplyr);library(tidyr);library(yarrr);library(stringr)
layout(1)

source("~/aarons_FAP_github_repository/recent_Annotation.R")
recent = recent_Annotation(go_online = "no")
recent = recent[recent$DNA_WGS!="",]
recent$SampleName = str_remove_all(recent$SampleName,"-")
recent$Size..mm.[recent$Stage..Polyp..Normal..AdCa.=="Normal"] = "" #Normal gets no size

# recent$Stage..Polyp..Normal..AdCa.
# recent$Location

drivers = read.csv("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv",
                   sep = ",",skip = 3,header = T)
driver.colors = data.frame("geneName" = sort(unique(drivers$Gene)),
           "Colors" = piratepal(palette = "basel",length.out = length(unique(drivers$Gene))))
# pt = "A001"
# patients = c("A001","A002");setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
# pt = "EP";setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")
patients = c("EP","JP");setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")

for (pt in patients) {
  if (pt == "EP" | pt == "JP") {
    ccf.table = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt"),
                           header = T,sep = "\t",stringsAsFactors = F)
    recent$SampleName = recent$VAP_Names
  }
  if (pt == "A001" | pt== "A002") {
    ccf.table = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"),
                           header = T,sep = "\t",stringsAsFactors = F)
  }
  ccf.table.cols = select(ccf.table,1:5,geneName,geneLoc,functionalClass,CADD_phred,ends_with("ccf"),-contains("bpccf"),-mergeCCF)
  ccf.cols = which(grepl("ccf",colnames(ccf.table.cols)))
  ccf.cols.name = grep("ccf",colnames(ccf.table.cols),value = T)
  ccf.cols.name.woccf = str_remove(string = ccf.cols.name,pattern = "ccf")

  #prepare to make at most 12 plots
  
  #CCF Table of just Driver genes
  ccf.table.cols.drivers = ccf.table.cols[(ccf.table.cols$geneName %in% drivers$Gene),]
  ccf.table.cols.drivers = filter(ccf.table.cols.drivers,geneLoc=="exonic")%>%
    filter(functionalClass!="synonymous SNV")
  
  #Save Image
  svg(filename = paste0(pt,"_ccf_Histogram_annotated_drivers.svg"),width = 12)
  
  par(mfrow=c(3,4),mar = c(4, 4, 3, 1))
  layout.show(n = 12)
  # i=8;layout(1);layout.show(n=1)
  for (i in 1:length(ccf.cols)) {
    #Select Sample CCF
    plot.me = ccf.table.cols[,ccf.cols[i]]
    
    names.me = ccf.table.cols.drivers[,c(6,9,ccf.cols[i])]
    nonzero.index = names.me[,ccf.cols.name[i]]>0
    names.me.nonzero = names.me[nonzero.index,]
    names.me.nonzero = names.me.nonzero[names.me.nonzero$CADD_phred>=20,]
    print(length(names.me.nonzero$geneName))
    
    ## Get clinical info about samples
    sample.index = which(recent$SampleName==ccf.cols.name.woccf[i])
    stage = recent$Stage..Polyp..Normal..AdCa.[sample.index]
    size = recent$Size..mm.[sample.index]
    location = recent$Location[sample.index]
    
    #remove EP and JP from the Identifiers and replace with F and G
    ccf.cols.name.woccf[i] = str_replace_all(string = ccf.cols.name.woccf[i],pattern = "EP",replacement = "F")
    ccf.cols.name.woccf[i] = str_replace_all(string = ccf.cols.name.woccf[i],pattern = "JP",replacement = "G")
    
    #Make Histogram
    h = hist(plot.me[plot.me>=0.1],
             xlim = c(0,1),breaks = 30,
             ylim = c(0,2000), 
             main = paste0(ccf.cols.name.woccf[i],"_",stage,"_",size,"_",location),
             xlab = "Cancer Cell Fraction",
             ylab = "Frequency",bty = "l")
    text.line = sort(h$counts,decreasing = T)[2] #height of labels is 2nd highest bar
    names.me.nonzero$y = text.line + 20*names.me.nonzero$CADD_phred #increase height of labels
    
    #add color column for each gene
    genes.colors = driver.colors[driver.colors$geneName %in% names.me.nonzero$geneName,]
    names.me.nonzero = full_join(names.me.nonzero,genes.colors,"geneName")
    
    #Is each gene a COAD driver gene?
    # gene = 1
    if (nrow(names.me.nonzero)!=0) {
      #create new column
      names.me.nonzero$COAD = NA
      for (gene in 1:nrow(names.me.nonzero)) {
        gene.index = which(grepl(pattern = names.me.nonzero$geneName[gene],x = drivers$Gene))
        if (sum(grepl(pattern = "COAD",x = drivers$Cancer[gene.index]))>0) {
          names.me.nonzero$COAD[gene] = "TRUE"
        } else {
          names.me.nonzero$COAD[gene] = "FALSE"
        }
      }
    }

    
    if (nrow(names.me.nonzero)!=0) {
      text(x = names.me.nonzero[,ccf.cols.name[i]],
           y = names.me.nonzero$y,
           labels = names.me.nonzero$geneName,
           offset = .5,
           srt = "45",
           col = as.vector(names.me.nonzero[,"Colors"]),
           font = ifelse(names.me.nonzero$COAD=="TRUE",2,3),
           cex = ifelse(names.me.nonzero$COAD=="TRUE",1.2,0.9))
    }
    next()
  }
  dev.off()
}

