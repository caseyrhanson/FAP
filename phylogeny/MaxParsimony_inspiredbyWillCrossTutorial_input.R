
library(dplyr);library(tidyr);library(stringr)

setwd("/Users/ahorning/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")

maxpars_input = function(pt,vap){
  #For A001,A002
  if (pt %in% c("A001", "A002")) {
    path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
  } else if (pt %in% c("EP","JP")){
    path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt")
  }
  path

  # is vap added as an argument or not? if yes, use that table, if not, find it based on the path
  if (!missing(vap)) {
    vap = vap
  } else {
    vap =  read.delim(file = path,header = T,sep = "\t")
  }
  
  
  #creates maf table and sample name lists
  maf = vapReduce_cover_vcover_mutect(vap = vap,patient = pt,mafOrccf = "ccf")[[1]]
  #make column of just original maf.mutID
  maf.mutID = unite(maf,col = "mutID",c("chr" ,"pos","id","ref","alt"),sep = ":")%>%
    select(mutID)
  
  # mafc = select(maf,"chr" ,"pos","id","ref","alt",ends_with("_mafc"))%>%
  #   unite(col = "mutID",c("chr" ,"pos","id","ref","alt"),sep = ":")
  # mafc = gather(data = mafc,key ="sample",value = "mafc",grep(pattern = "_mafc",colnames(mafc)))
  
  
  d = select(maf,"chr" ,"pos","id","ref","alt",ends_with("d"))%>%
    unite(col = "mutID",c("chr" ,"pos","id","ref","alt"),sep = ":")
  d = gather(data = d,key ="sample",value = "d",grep(pattern = "mutID",colnames(d),invert = T))
  
  somatic = select(maf,"chr" ,"pos","id","ref","alt",ends_with("_somatic"))%>%
    unite(col = "mutID",c("chr" ,"pos","id","ref","alt"),sep = ":")
  somatic = gather(data = somatic,key ="sample",value = "somatic",
                   grep(pattern = "mutID",colnames(somatic),invert = T))
  
  
  data = data.frame(mutID = somatic$mutID,
                    sample = str_remove(somatic$sample,"_somatic"),
                    # mafc = round(mafc$mafc * d$d,digits = 0),
                    d = d$d,
                    somatic = somatic$somatic,
                    presence = rep("absent",length(d$d)),stringsAsFactors = F)
  
  # Heuristic Thresholds - ignore these and just use Somatic Column
  # data$presence[data$mutect == "yes"] = "present"
  # data$presence[data$mutect=="unknown" & data$d >= 10] = "absent"
  # data$presence[data$mutect=="unknown" & data$d < 10 & data$mafc >= 1] = "present"
  # 
  # data$presence[data$mutect == "unknown" & data$d >= 10 & data$mafc >=3] = "present"
  # data$presence[data$mutect == "unknown" & data$d >= 10 & data$mafc == 0] = "absent"
  # data$presence[data$mutect == "unknown" & data$d <= 10 & data$mafc >= 1] = "present"
  # data$presence[data$mutect == "unknown" & data$d >= 20 & data$mafc == 1] = "absent"
  
  #Filter out only mutations that are Present
  data = data[data$somatic=="present",]
  data$presence = data$somatic #Somatic mutations are considerd present
  
  # remake original ccf format. Fill in "spread" with "absent" because the tall format was 
  # filtered above for present. So when it is spread, it will have holes for "absent" mutations
  presence = select(data,mutID,sample, presence)%>%
    spread(key = sample,value = presence,fill = "absent")
  
  #to account for all the original CCF file mutations, align these filtered mutations  
  presence.maf.mutID = left_join(x = maf.mutID,y = presence,"mutID")
  
  # prepare the dataframe for binary.
  rownames(presence.maf.mutID) = presence.maf.mutID[,1]
  data.input = data.frame(presence.maf.mutID[,-1],stringsAsFactors = F)
  
  #Make binary input: 1 is present, 0 is absent
  data.input[data.input=="present"] = 1
  data.input[data.input=="absent"] = 0

  # the rows added back in the left_join were missing because the mutations were absent
  # there were NA so i called them "0" below
  data.input[is.na(data.input)] = 0 

      
  for (i in 1:dim(data.input)[2]) {
    data.input[,i] = as.numeric(as.character(data.input[,i]))
  }
  
  
  #Add a column for the Normal. All somatic mutations should be 0 (absent)
  data.input$Normal = 0
  
  str(data.input)
  return(data.input)
}


