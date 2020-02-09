# rm(list = ls())
library(dplyr);library(tidyr);library(stringr)

setwd("/Users/ahorning/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")

# pt = "EP"



dnds_input = function(pt){
  #For A001,A002
  if (c(pt) %in% c("A001", "A002")) {
    # CCF output
    path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
    # VAP output
    # path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table")
  } else if (c(pt) %in% c("EP","JP")){
    path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt")
  } else {}
  path
  
  # Because this VAP has germline column and more than 1 patient in it, filter out columns and rows
  vap =  read.delim(file = path,header = T,sep = "\t")
  # new.vap = vap
  # old.vap = vap
  
  # Filter out the germline mutations
  # index.germline = !is.na(vap$germline)
  # x = vap[!index.germline,]
  # 
  # #samples
  # colnames(x)
  # For filtering the rows for patient specific somatic mutations
  # numsamples = length(samples)
  # rindexTF = vector()
  # cindex = vector()
  # for (si in 1:length(samples)) {
  #   if (length(rindexTF) > 0) {
  #     rindexTF = rindexTF | grepl(paste(samples[si],"\\[",sep=""), d$somatic)
  #   } else {
  #     rindexTF = grepl(paste(samples[si],"\\[",sep=""), d$somatic)
  #   }
  #   cindex = append(cindex, match(paste(samples[si], "maf", sep=""), colnames(d)))
  # }
  # rindex = which(rindexTF)
  
  
  #creates maf table and sample name lists
  # looks similar to VAP output but with less columns and has Somatic column for Present and Absent mutations
  maf = vapReduce_cover_vcover_mutect(vap = vap,patient = pt,mafOrccf = "ccf")[[1]]
  
  # 
  somatic = select(maf,"chr" ,"pos","id","ref","alt",ends_with("_somatic"))%>%
    unite(col = "mutID",c("chr" ,"pos","id","ref","alt"),sep = ":")
  #Add functional class to table
  somatic$functionalClass = vap$functionalClass
  somatic$geneName = vap$geneName
  # Melt the table
  somatic = gather(data = somatic,key ="sample",value = "somatic",
                   grep(pattern = "mutID|functionalClass|geneName",colnames(somatic),invert = T))
  
  
  data = data.frame(mutID = somatic$mutID,
                    sample = str_remove(somatic$sample,"_somatic"),
                    # mafc = round(mafc$mafc * d$d,digits = 0),
                    functionalClass = somatic$functionalClass,
                    geneName = somatic$geneName,
                    presence = somatic$somatic)

  head(data)
  table(data$sample)
  #Filter out only present mutations.
  data = data[data$presence=="present",]
  table(data$sample)
  
  # Select and rename final columns for 
  # The functionalClass, geneName and presence columns will need to be reomved before dnds
  data_input = separate(data,col = mutID,into = c("chr" ,"pos","id","ref","alt"),sep = ":")%>%
    select("sampleID" = sample,chr,pos,ref, "mut" = alt,
           functionalClass, geneName, presence)
  table(data$sample)
  head(data_input)
  
  return(data_input)
}


# example of dnds input
#   sampleID chr      pos ref mut
# 1 Sample_1   1   871244   G   C
# 2 Sample_1   1  6648841   C   G


#remade booleans from Zheng
# These Heuristics are unnecessary because the "somatic" columns are the filters
# data$presence[data$mutect == "yes"] = "present"
# data$presence[data$mutect=="unknown" & data$d >= 10] = "absent"
# data$presence[data$mutect=="unknown" & data$d < 10 & data$mafc >= 1] = "present"
# 
# data$presence[data$mutect == "unknown" & data$d >= 10 & data$mafc >=3] = "present"
# data$presence[data$mutect == "unknown" & data$d >= 10 & data$mafc == 0] = "absent"
# data$presence[data$mutect == "unknown" & data$d <= 10 & data$mafc >= 1] = "present"
# data$presence[data$mutect == "unknown" & data$d >= 20 & data$mafc == 1] = "absent"

#Booleans from Zheng to detect "present" and "absent" mutations
# source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")

# #####done
# if(data$mutect == "yes"){
#   data$presence = "present"
# }
# # {allele = toString(var_allele[k])}
# 
# else if(mutect == "no") {
#   if(cover>=10)
#   {allele = toString(ref_allele[k])} #if total coverage of the reference allele is >10, then ref.
#   else{
#     if(vcover>=1) #else, if variant cover is greater than 1 (thus > 10% coverage)
#     {allele = toString(var_allele[k])} #variant allele
#     else
#     {allele = toString("N")} #else, unknown.
#   }
# }
# 
# else{ #if mutect is neither "yes" nor "no", and allele call is unclear, then...
#   if(cover>=10 & vcover>=3)
#   {allele = toString(var_allele[k])}
#   else if(cover>=10 & vcover==0)
#   {allele = toString(ref_allele[k])}
#   else if(cover<=10 & vcover>=1)
#   {allele = toString(var_allele[k])}
#   else if(cover>=20 & vcover==1)
#   {allele = toString(ref_allele[k])}
#   else
#   {allele = toString("N")}
# }
# 
