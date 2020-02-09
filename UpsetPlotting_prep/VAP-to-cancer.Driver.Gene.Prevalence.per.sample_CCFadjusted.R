rm(list=ls())

library(stringr);library(dplyr);library(tidyr)

source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")
source("~/aarons_FAP_github_repository/vapReduce_cover,ccf,somatic.R")
source("~/DNAseq_WGS/scripts/find.sample.and.tumor.names.R")
source("~/aarons_FAP_github_repository/phylogeny/MaxParsimony_inspiredbyWillCrossTutorial_input.R")

load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
ccf.data = list("A001" = ptA001,"A002" = ptA002,"EP" = ptF,"JP" = ptG)

pt = "JP"
pts = c("EP","JP","A001","A002")
pts = c("A001","A002")
for (pt in pts) {
  #read in VAP output file
  # if (pt %in% c("A001","A002")) {
  #   vap.path = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
  # } else if(pt %in% c("EP","JP")){
  #   vap.path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt") #EP and JP original WGS
  # }
  # vap.path
  
  vap = ccf.data[[pt]]
  
  #read in Comprehensive Cancer Drive list
  pan.cancer.driver.path = path.expand("~/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv")
  drivers= read.csv(pan.cancer.driver.path,header = T,skip = 3,stringsAsFactors = F)
  
  # #create list of samples from VAP file
  # patient = "*" #collects both A001 and A002
  # names=find.sample.and.tumor.names(vap,patient)
  # sample.names = names$sample.names
  # normal.names=names$normal.names
  # tumor.names=names$tumor.names
  
  #Filter only the rows in the VAP output which have geneNames in the PanCancer Driver gene list AND are non-synonymous or stop-gain Muts
  vap.drivers.only = vap[which(vap$geneName %in% drivers$Gene),]
  #Filter NA and synonymous mutations from the dataset.
  #Filter for top 1% (CADD >= 20) of deleterious mutations
  vap.drivers.only = filter(vap.drivers.only,functionalClass!="synonymous SNV")%>%filter(CADD_phred>=20)%>%
    unite(mutID,c("chr","pos","id","ref","alt"),sep = ":")%>%
    select(mutID,geneName,AAChange)%>%
    separate(col = AAChange,into = c("AAChange",NA),sep = ",")
  
  #Add chr to chr column
  if(!grepl(pattern = "chr",x = vap.drivers.only$mutID[1])){
    vap.drivers.only$mutID = paste0("chr",vap.drivers.only$mutID)
  } else {}
  
  head(vap.drivers.only)
  
  # read in data from the CCF version of the VAP file
  vap.ccf.presence = maxpars_input(pt = pt)
  vap.ccf.presence$mutID = rownames(vap.ccf.presence)
  
  if(!grepl(pattern = "chr",x = vap.ccf.presence$mutID[1])){
    vap.ccf.presence$mutID = paste0("chr",vap.ccf.presence$mutID)
  } else {}
  
  head(vap.ccf.presence)
  
  #Column bind the Driver+NonSynonymous Mutations with the "Present/Absent" Matrix
  x = inner_join(x = vap.ccf.presence,y = vap.drivers.only,"mutID")
  x = data.frame(x,stringsAsFactors = F)
  str(x)
  
  for (r in 1:nrow(x)) {
    for(c in 1:grep(pattern = "Normal",x = colnames(x))){
      if (x[r,c]==1) {
        x[r,c]= x[r,"AAChange"]
      } else {}
    }
  }
  
  x
  
  write.table(x = x,
              file = paste0("~/Bulk_A001_A002/repeatMutations/mutect.snv.res...somatic.table.ccf_CADD20_NonSyn_",pt,".txt"),
              append = F,sep = "\t",row.names = F)
  
}

