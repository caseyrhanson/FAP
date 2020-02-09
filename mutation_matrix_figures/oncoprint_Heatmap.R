rm(list=ls())

library(dplyr);library(tidyr);library(ComplexHeatmap);library(stringr);library(circlize);library(RColorBrewer);library(yarrr)

setwd("~/Bulk_A001_A002/")

load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
# muts = read.table("mutect.snv.res.filtered.classified.founds.nopara.somatic.table",
#            header = T,sep = "\t",stringsAsFactors = F)
# 
# muts.FG = read.table("../mutect.snv.res.filtered.classified.founds.nopara.somatic.table",
#                      header = T,sep = "\t",stringsAsFactors = F)

pan.cancer.driver.path = path.expand("~/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv")
drivers = read.csv(pan.cancer.driver.path,header = T,skip = 3,stringsAsFactors = F)

vap.selectCol.filter = function(VAPoutput){
  # practice
  # VAPoutput = muts
  
  #Filter VAP file for specific mutations of interest
  muts.cols = unite(VAPoutput,mutID,c("chr",  "pos", "id",   "ref",  "alt"))%>%
    dplyr::select(mutID,somatic, geneName, functionalClass, geneLoc)%>%
    filter(geneLoc=="exonic")%>%
    arrange(geneName)%>%
    mutate(functionalClass = str_replace_all(functionalClass," ",""))%>%
    mutate(geneLoc = str_replace_all(geneLoc," ",""))
  
  #Find drivers specific to COADREAD
  drivers.small = drivers%>%
    dplyr::select(Gene,Cancer,Tumor.suppressor.or.oncogene.prediction..by.20.20..,Decision)
  
  
  #Find all the cancer driver genes in the VAP output folder
  muts.cols.drivers = muts.cols[muts.cols$geneName %in% drivers$Gene,]
  gene.drivers = sort(unique(muts.cols.drivers$geneName))
  
  #Tumor suppressors or oncogenes: classifying the mutations i collected
  tsg.or.onc.drivers = drivers.small[which(drivers.small$Gene %in% gene.drivers),]
  #specific to coadread
  tsg.or.onc.coadread.drivers = filter(tsg.or.onc.drivers,Cancer=="COADREAD")%>%
    dplyr::rename("Tsg.or.Onc" = "Tumor.suppressor.or.oncogene.prediction..by.20.20..")
  #just the genes with coadread-specific tsg or onc status
  tsg.or.onc.coadread.genes = filter(tsg.or.onc.drivers,Cancer=="COADREAD")[,1]
  #table of genes without COADREAD-specific driver status
  tsg.or.onc.noncoadread.drivers = tsg.or.onc.drivers[!(tsg.or.onc.drivers$Gene %in% tsg.or.onc.coadread.genes),]
  colnames(tsg.or.onc.noncoadread.drivers)[3] = "Tsg.or.Onc"
  tsg.or.onc.drivers = rbind(tsg.or.onc.coadread.drivers,tsg.or.onc.noncoadread.drivers)
  
  return(list(muts.cols.drivers,gene.drivers,tsg.or.onc.drivers))
}

#Put all Driver mutations together in same table
muts.cols.drivers = rbind(vap.selectCol.filter(VAPoutput = ptA001)[[1]],
                          vap.selectCol.filter(VAPoutput = ptA002)[[1]],
                          vap.selectCol.filter(VAPoutput = ptF)[[1]],
                          vap.selectCol.filter(VAPoutput = ptG)[[1]])

#Filter out synonymousSNVs
muts.cols.drivers = filter(muts.cols.drivers,functionalClass!="synonymousSNV")

#Collect all gene Drivers from the samples. Create Unique and sorted list
gene.drivers = sort(unique(c(as.vector(vap.selectCol.filter(VAPoutput = ptA001)[[2]]),
                             as.vector(vap.selectCol.filter(VAPoutput = ptA002)[[2]]),
                             as.vector(vap.selectCol.filter(VAPoutput = ptF)[[2]]),
                             as.vector(vap.selectCol.filter(VAPoutput = ptG)[[2]]))))

#only run when/if new genes are added to the list of mutations in the samples we have. 
# tsg.or.onc = rbind(vap.selectCol.filter(muts)[[3]], vap.selectCol.filter(muts.FG)[[3]])
# tsg.or.onc.allMuts = arrange(unique(tsg.or.onc),Gene)
# write.table(x = tsg.or.onc.allMuts,
#             file = "drivers_TsgOrOncogenic.txt",append = F,quote = F,
#             sep = "\t",row.names = F)
# write.table(x = muts.cols.drivers,
#             file = "drivers_mutationTypes.txt",append = F,quote = F,
#             sep = "\t",row.names = F)


#Make a list of the samples which are in the VAP output file
ccf.data = list(ptA001,ptA002,ptF,ptG)
ccf.sets = length(ccf.data)

samples = c()
for (set in 1:ccf.sets) {
  ccf = ccf.data[[set]]
  sample.index = which(grepl("ccf$",colnames(ccf)))
  samples = c(samples,str_remove_all(colnames(ccf)[sample.index],pattern = "ccf$"))
}

samples

#Make empty matrix prepared for entering all the mutation information
onc.matrix = matrix(data = "",nrow = length(gene.drivers),ncol = length(samples))
colnames(onc.matrix) = samples #has "." instead of "-"
rownames(onc.matrix) = gene.drivers


for (i in 1:nrow(muts.cols.drivers)) {
  sample.bracket = head(str_split(muts.cols.drivers$somatic[i],pattern = ",")[[1]],-1)
  for (sample in sample.bracket) {
    if (grepl("doubt",sample)) {
      next()
    } else {
      sample.name = str_split(sample,"\\[")[[1]][1]
      sample.name = str_replace_all(sample.name,"-",".")
    }
      gene.name = as.character(muts.cols.drivers$geneName[i])
      mut.type = muts.cols.drivers$functionalClass[i] #exonic mutation
      # mut.type = ifelse(is.na(mut.type),muts.cols.drivers$geneLoc[i],mut.type) #non-exonic mutation assigned
      
      if(onc.matrix[gene.name,sample.name]==""){
        onc.matrix[gene.name,sample.name] = mut.type 
      } else if(onc.matrix[gene.name,sample.name]!=""){
        mut.type1 = onc.matrix[gene.name,sample.name] #put previous mutation type in object
        onc.matrix[gene.name,sample.name] = paste0(mut.type,";",mut.type1) #paste each mutation type together
      }
       
    }
    
  }

#########
# Rearrange order of Columns so All the Patients are together and ordered by Number of Mutations
onc.matrix = onc.matrix[,names(sort(colSums(onc.matrix!=""),decreasing = T))]
onc.matrix.less = onc.matrix[,c(grep("A001",colnames(onc.matrix),value = T),
           grep("A002", colnames(onc.matrix),value = T),
           grep("EP",colnames(onc.matrix),value = T),
           grep("JP",colnames(onc.matrix),value = T))]

#removes genes/rows without mutations
onc.matrix.less = onc.matrix.less[rowSums(onc.matrix.less!="")>0,] 

#Plot Mutations
get_type_fun = function(x) strsplit(x, ";")[[1]] #Separating the types of mutations

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                              gp = gpar(fill = "#BDBDBD")),
  nonsynonymousSNV = function(x,y,w,h)  grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["nonsynonymousSNV"], col = NA)),
  stopgain = function(x,y,w,h)  grid.rect(x, y, w*0.9, h*0.4, 
                                          gp = gpar(fill = col["stopgain"], col = NA)),
  synonymousSNV = function(x,y,w,h) grid.rect(x, y, w*0.4, h*0.9, 
                                              gp = gpar(fill = col["synonymousSNV"], col = NA)),
  activating = function(x, y, w, h) grid.points(x, y, pch = 16),
  inactivating = function(x, y, w, h) grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5,
                                                    gp = gpar(lwd = 2))
)

# Mutation colors
col = c(nonsynonymousSNV = "red", stopgain = "blue",synonymousSNV = "purple")


#draw initial oncoplot
onc.plot = draw(oncoPrint(mat = onc.matrix.less,get_type = get_type_fun,
                          alter_fun = alter_fun,
                          col = col,show_column_names = T, column_order = colnames(onc.matrix.less)))

#read in the sample table with the sample names and its associated tissue types
source("~/aarons_FAP_github_repository/recent_Annotation.R")
# anno = recent_Annotation(go_online_clinical = "yes",go_online_pathology = "yes")
anno = recent_Annotation("no")
samples.tissue.types = data.frame(Samples = str_remove_all(string = anno$VAP_Names,pattern = "-"),
                                  Patient = anno$Patient,
                                  Stage = anno$Stage..Polyp..Normal..AdCa.,
                                  Size = sapply(str_split(anno$Size..mm.,"x",simplify = F),"[",1),
                                  Location = anno$Location,
                                  NeoCellsInTumor = anno$NeoCellsInTumor,
                                  TumorInTotal = anno$TumorInTotal,
                                  nonNeoStromaInTotal = anno$nonNeoStromaInTotal,
                                  PercentNormalOverall = anno$PercentNormalOverall,
                                  PercentStromaOverall = anno$PercentStromaOverall,
                                  PercentCancerOverall = anno$PercentCancerOverall,
                                  PercentAdenomaOverall = anno$PercentAdenomaOverall,
                                  PercentStromaInCancer = anno$PercentStromaInCancer,
                                  PercentStromaInAdenoma = anno$PercentStromaInAdenoma,
                                  PercentNecrosisInTumor = anno$PercentNecrosisInTumor,
                                  stringsAsFactors = F)
#Add Purity and Ploidy Values to Annotation
load("~/aarons_FAP_github_repository/ccf_estimation/titan_PloidyAndPurity_data.Rdata")
pp = as.data.frame(pp.matrix)
pp$Samples = rownames(pp)
samples.tissue.types = full_join(samples.tissue.types,pp,"Samples")

#If Tissue is Normal, remove Size information
samples.tissue.types$Size[samples.tissue.types$Stage=="Normal"] = "" 

##Determine order of columns of heatmap
#get column order
onc.plot.col.order = colnames(onc.matrix.less)
#reorganize clinical data by columns in heatmap
samples.tissue.types = samples.tissue.types[match(onc.plot.col.order,samples.tissue.types$Samples),] 


#Colors for continuous Size variable
#Makes size a numeric values
samples.tissue.types$Size = as.numeric(samples.tissue.types$Size)

#colors for continuous variables
col_size = colorRamp2(c(min(samples.tissue.types$Size,na.rm = T),
                        max(samples.tissue.types$Size,na.rm = T)),
                      c("lightgreen", "navy"))

#colors for discrete variables
samples.tissue.types$Patient = str_replace_all(samples.tissue.types$Patient,pattern = "EP",replacement = "F")
samples.tissue.types$Patient = str_replace_all(samples.tissue.types$Patient,pattern = "JP",replacement = "G")
col_patient = brewer.pal(n = length(unique(samples.tissue.types$Patient)), name = "Accent")
names(col_patient) = unique(samples.tissue.types$Patient)

col_stage = brewer.pal(n = length(unique(samples.tissue.types$Stage)),name = "Reds")
names(col_stage) = c("Normal","Polyp",  "AdCa") #unique(samples.tissue.types$Stage)

#Make color vector for Location
col_location = piratepal(palette = "southpark",length.out = length(unique(samples.tissue.types$Location)))
names(col_location) = c("Ascending","Transverse", "Descending", "Rectum")
# unique(samples.tissue.types$Location)

#Make Continuous Colors for Neoplastic Cells in Tumor
col_neocellsintumor = colorRamp2(c(min(samples.tissue.types$NeoCellsInTumor,na.rm = T),
                                   max(samples.tissue.types$NeoCellsInTumor,na.rm = T)),
                                 c("orange", "darkred"))

#Make Continuous Colors for Tumor in Total
col_TumorInTotal = colorRamp2(c(min(samples.tissue.types$TumorInTotal,na.rm = T),
                                max(samples.tissue.types$TumorInTotal,na.rm = T)),
                              c("lightblue", "black"))

#Make Continuous Colors for Neoplastic Cells in Tumor
col_neocellsintumor = circlize::colorRamp2(c(min(samples.tissue.types$NeoCellsInTumor,na.rm = T),
                                             max(samples.tissue.types$NeoCellsInTumor,na.rm = T)),
                                           c("orange", "darkred"))

#Make Continuous Colors for Tumor in Total
col_TumorInTotal = circlize::colorRamp2(c(min(samples.tissue.types$TumorInTotal,na.rm = T),
                                          max(samples.tissue.types$TumorInTotal,na.rm = T)),
                                        c("lightblue", "black"))


# Make Continuous Colors for PercentNormalOverall
# PercentNormalOverall light greeen dark green
col_PercentNormalOverall = circlize::colorRamp2(c(min(samples.tissue.types$PercentNormalOverall,na.rm = T),
                                                  max(samples.tissue.types$PercentNormalOverall,na.rm = T)),
                                                c("lightgreen", "darkgreen"))

# Make Continuous Colors for Percent Stroma Overall
# "PercentStromaOverall"  
col_PercentStromaOverall = circlize::colorRamp2(c(min(samples.tissue.types$PercentStromaOverall,na.rm = T),
                                                  max(samples.tissue.types$PercentStromaOverall,na.rm = T)),
                                                c("lightyellow1", "yellowgreen"))


# Make Continuous Colors for Percent Cancer Overall
# "PercentCancerOverall"   
col_PercentCancerOverall = circlize::colorRamp2(c(min(samples.tissue.types$PercentCancerOverall,na.rm = T),
                                                  max(samples.tissue.types$PercentCancerOverall,na.rm = T)),
                                                c("mediumpurple1", "mediumpurple4"))

# Make Continuous Colors for Percent Stroma in Cancer
# "PercentStromaInCancer"  
col_PercentStromaInCancer = circlize::colorRamp2(c(min(samples.tissue.types$PercentStromaInCancer,na.rm = T),
                                                   max(samples.tissue.types$PercentStromaInCancer,na.rm = T)),
                                                 c("mediumpurple4","yellowgreen"))

# Make Continuous Colors for Percent Necrosis in Tumor
col_PercentNecrosisInTumor = circlize::colorRamp2(c(min(samples.tissue.types$PercentNecrosisInTumor,na.rm = T),
                                                    1),
                                                  c("grey89","red4"))


# Make Continuous Colors for Percent Adenoma Overall
# "PercentAdenomaOverall" 
col_PercentAdenomaOverall = circlize::colorRamp2(c(min(samples.tissue.types$PercentAdenomaOverall,na.rm = T),
                                                   max(samples.tissue.types$PercentAdenomaOverall,na.rm = T)),
                                                 c("lightpink","deeppink3"))

# Make Continuous Colors for PercentStromaInAdenoma
# "PercentStromaInAdenoma"
col_PercentStromaInAdenoma = circlize::colorRamp2(c(min(samples.tissue.types$PercentStromaInAdenoma,na.rm = T),
                                                    max(samples.tissue.types$PercentStromaInAdenoma,na.rm = T)),
                                                  c("yellowgreen","deeppink3"))
########


# ha = HeatmapAnnotation(which = "column",show_annotation_name = T,
#                        df = samples.tissue.types[,2:length(colnames(samples.tissue.types))],
#                        col = list(Patient = col_patient,
#                                   Stage = col_stage,
#                                   Size = col_size,
#                                   Location = col_location,
#                                   NeoCellsInTumor = col_neocellsintumor,
#                                   TumorInTotal = col_TumorInTotal,
#                                   PercentNormalOverall = col_PercentNormalOverall,
#                                   PercentStromaOverall = col_PercentStromaOverall,
#                                   PercentCancerOverall = col_PercentCancerOverall,
#                                   PercentStromaInCancer = col_PercentStromaInCancer,
#                                   PercentNecrosisInTumor = col_PercentNecrosisInTumor),
#                        na_col = "white",
#                        gp = gpar(col = "black"))

ha = HeatmapAnnotation(which = "column",show_annotation_name = T,
                       df=samples.tissue.types[,c("Patient","Stage","Size","Location",
                                                  "NeoCellsInTumor","TumorInTotal","Purity", "Ploidy")],
                       col = list(Patient = col_patient,
                                  Stage = col_stage,
                                  Size = col_size,
                                  Location = col_location,
                                  NeoCellsInTumor = col_neocellsintumor,
                                  TumorInTotal = col_TumorInTotal),
                       na_col = "white",
                       gp = gpar(col = "black"))


#rename sample names: EP to F and JP to G
colnames(onc.matrix.less) = str_replace_all(string = colnames(onc.matrix.less),pattern = "EP",replacement = "F")
colnames(onc.matrix.less) = str_replace_all(string = colnames(onc.matrix.less),pattern = "JP",replacement = "G")

######### Create the Oncoplot with Annotation ############
onc.plot = oncoPrint(mat = onc.matrix.less,get_type = get_type_fun,
                          alter_fun = alter_fun,
                          col = col,show_column_names = T,column_order = colnames(onc.matrix.less),
                          pct_gp = gpar(cex = 0.8),column_names_gp = gpar(cex=0.6),row_names_gp = gpar(cex = 0.8),
                          bottom_annotation = ha,
                          width = unit(x = 12,units = "cm"))
onc.plot

#Save the plot
pdf(file = "A001_A002_F_G_oncoplot_annotations.pdf",width = 14, height = 20)
onc.plot
dev.off()

pdf(file = "A001_A002_F_G_oncoplot_heatmap.pdf",width = 14, height = 10)
onc.plot
dev.off()
