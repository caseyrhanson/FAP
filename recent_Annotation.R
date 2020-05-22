recent_Annotation = function(go_online_clinical, go_online_pathology){
  library(dplyr);library(tidyr)
  
  # Paths
  clinical = "https://docs.google.com/spreadsheets/d/e/2PACX-1vRsC3CJUE-kFKdfrA5CB-K8V4896QyvVGONRk2B9BEqv6bwf3ovAtJ_W4nJmhfo18l61AOH7QztmYEO/pub?gid=0&single=true&output=csv"
  pathology.path = "https://docs.google.com/spreadsheets/d/e/2PACX-1vTZVugWQtSvlTlUMIbUhjT7dWspkiMu_EiellVkbL3XbO9mS6neobwd9tfBP_ZljpupFOUK4U7y7Rv9/pub?gid=1098597315&single=true&output=tsv"
  
  original.wd = getwd()
  setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/SampleInformation/")
  
  if (go_online_clinical=="yes" && go_online_pathology!="yes") {
    recent = read.csv(clinical, stringsAsFactors = F,header = T)
    write.table(x = recent,file = paste0("SampleTrackingSheet_Roxanne_",format(Sys.time(),"%b%d%Y"),".txt"),sep = "\t",row.names = F)
  
  } else if (go_online_clinical=="yes" && go_online_pathology=="yes") {
    
    #Clinical annotations
    recent = read.csv(clinical, stringsAsFactors = F,header = T)
    recent = separate(data = recent,col = "SampleName",into = "SampleName",sep = " ")
    
    
    #Pathology
    pathology = read.table(pathology.path, header = T,sep = "\t")
  
    # num.samples = length(unique(pathology$Sample))-2
    pathology = pathology[-1,]
    # pathology = pathology[c(1:num.samples),]
    
    colnames(pathology)
    
    summary(pathology)
    #Select and Rename only specific columns from Pathology Report
    pathology = dplyr::select(pathology,"SampleName" = "Sample",
                              "isPolyp"="Polyp",
                              "PolypType"="Polyp.Type",
                              "Dysplasia",
                              "HighDysplasia"="High","LowDysplasia"="Low",
                       "NeoCellsInTumor"  = "Percent.of.neoplastic.cell.nuclei.as.a.total.of.all.cell.nuclei.in.tumor.ONLY.area",
                       "nonNeoStromaInTotal"  = "Percent.of.non.neoplastic.stroma..by.area",   
                       "TumorInTotal" = "Percent.of.entire.tissue.involved.by.tumor",
                      "PercentNormalOverall" = "Percent.Normal.Overall",                                        
                      "PercentStromaOverall" = "Percent.Stroma.Overall",
                      "PercentCancerOverall" = "Percent.Cancer.Overall",
                      "PercentAdenomaOverall" = "Percent.Adenoma.Overall",
                      "PercentStromaInCancer" = "Percent.Stroma.in.Cancer",
                      "PercentStromaInAdenoma" = "Percent.Stroma.in.Adenoma",
                      "PercentNecrosisInTumor" = "Percent.Necrosis.in.Tumor")%>%
      separate(col = "SampleName",into = "SampleName",sep = " ")
    
    #convert specific columns from character to numeric
    idx.first.char2num = grep(pattern = "NeoCellsInTumor",x = colnames(pathology))
    cols.for.char.2.numeric = colnames(pathology)[idx.first.char2num:length(colnames(pathology))]
    
    for (cols in cols.for.char.2.numeric) {
      # str(pathology[,cols])
      pathology[,cols] = as.numeric(as.character(pathology[,cols]))
      # str(pathology[,cols])
    }
    str(pathology)
    pathology[pathology==""] = NA
    
    # pathology$NeoCellsInTumor = as.numeric(levels(pathology$NeoCellsInTumor))[pathology$NeoCellsInTumor]
    
    #Combine Clinical and Pathology information
    recent = full_join(x = recent,y = pathology,"SampleName")
    
    write.table(x = recent,file = paste0("SampleTrackingSheet_Roxanne_Pathology",format(Sys.time(),"%b%d%Y"),".txt"),sep = "\t",row.names = F)
    
  } else {
    #Just read in the latest Clinical Annotation and Pathology Files
    files = list.files(pattern = "SampleTrackingSheet",)
    recent = read.table(file = files[which.max(file.mtime(files))],
                        header = T,
                        sep = "\t",stringsAsFactors = F)
  }
  
  setwd(original.wd)
  getwd()
  return(recent)
} 



