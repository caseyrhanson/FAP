rm(list=ls())

library(dplyr);library(tidyr);library(ComplexHeatmap);library(stringr);library(circlize);library("RColorBrewer");library(yarrr)

setwd("~/Bulk_A001_A002/")

muts = read.table("mutect.snv.res.filtered.classified.founds.nopara.somatic.table",
           header = T,sep = "\t",stringsAsFactors = F)

muts.FG = read.table("../mutect.snv.res.filtered.classified.founds.nopara.somatic.table",
                     header = T,sep = "\t",stringsAsFactors = F)

pan.cancer.driver.path = path.expand("~/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv")
drivers = read.csv(pan.cancer.driver.path,header = T,skip = 3,stringsAsFactors = F)

vap.selectCol.filter = function(VAPoutput){
  #Filter VAP file for specific mutations of interest
  muts.cols = select(VAPoutput,somatic, geneName, functionalClass, geneLoc)%>%
    filter(geneLoc=="exonic")%>%
    arrange(geneName)%>%
    mutate(functionalClass = str_replace_all(functionalClass," ",""))%>%
    mutate(geneLoc = str_replace_all(geneLoc," ",""))
  
  #Find all the cancer driver genes in the VAP output folder
  muts.cols.drivers = muts.cols[muts.cols$geneName %in% drivers$Gene,]
  gene.drivers = sort(unique(muts.cols.drivers$geneName))
  
  return(list(muts.cols.drivers,gene.drivers))
}

muts.cols.drivers = rbind(vap.selectCol.filter(muts)[[1]], vap.selectCol.filter(muts.FG)[[1]])
gene.drivers = sort(unique(c(vap.selectCol.filter(muts)[[2]],vap.selectCol.filter(muts.FG)[[2]])))

#Make a list of the samples which are in the VAP output file
sample.index = which(grepl("maf",colnames(muts)))
samples = str_remove_all(colnames(muts)[sample.index],pattern = "maf")
sample.index = which(grepl("maf",colnames(muts.FG)))
samples2 = str_remove_all(colnames(muts.FG)[sample.index],pattern = "maf")

samples = c(samples,samples2)

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
      gene.name = muts.cols.drivers$geneName[i]
      mut.type = muts.cols.drivers$functionalClass[i] #exonic mutation
      mut.type = ifelse(is.na(mut.type),muts.cols.drivers$geneLoc[i],mut.type) #non-exonic mutation assigned
      
      if(onc.matrix[gene.name,sample.name]==""){
        onc.matrix[gene.name,sample.name] = mut.type 
      } else if(onc.matrix[gene.name,sample.name]!=""){
        mut.type1 = onc.matrix[gene.name,sample.name] #put previous mutation type in object
        onc.matrix[gene.name,sample.name] = paste0(mut.type,";",mut.type1) #paste each mutation type together
      }
       
    }
    
  }

#########

onc.matrix = onc.matrix[,names(sort(colSums(onc.matrix!=""),decreasing = T))]
onc.matrix.less = onc.matrix[,c(grep("A001",colnames(onc.matrix),value = T),
           grep("A002", colnames(onc.matrix),value = T),
           grep("EP",colnames(onc.matrix),value = T),
           grep("JP",colnames(onc.matrix),value = T))]

onc.matrix.less = onc.matrix.less[rowSums(onc.matrix.less!="")>0,] #removes genes/rows without mutations

#Plot Mutations
get_type_fun = function(x) strsplit(x, ";")[[1]] #Separating the types of mutations

alter_fun = list(
  background = function(x, y, w, h) grid.rect(x, y, w, h, 
                                              gp = gpar(fill = "#BDBDBD")),
  nonsynonymousSNV = function(x,y,w,h)  grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["nonsynonymousSNV"], col = NA)),
  stopgain = function(x,y,w,h)  grid.rect(x, y, w*0.9, h*0.4, 
                                          gp = gpar(fill = col["stopgain"], col = NA)),
  synonymousSNV = function(x,y,w,h) grid.rect(x, y, w*0.4, h*0.9, gp = gpar(fill = col["synonymousSNV"], col = NA))
)

#mutation colors
col = c(nonsynonymousSNV = "red", stopgain = "blue",synonymousSNV = "purple")
# x = grep("pink",colors(),value = T)
# plot(1:length(x),1:length(x),col =x,pch = 17)

#draw initial oncoplot
onc.plot = draw(oncoPrint(mat = onc.matrix.less,get_type = get_type_fun,
                          alter_fun = alter_fun,
                          col = col,show_column_names = T,column_order = NULL))

#read in the sample table with the sample names and its associated tissue types
source("~/aarons_FAP_github_repository/recent_Annotation.R")
anno = recent_Annotation(go_online_clinical = "yes",go_online_pathology = "yes")
samples.tissue.types = data.frame(Samples = str_remove_all(string = anno$VAP_Names,pattern = "-"),
                                  Patient = anno$Patient,
                                  Stage = anno$Stage..Polyp..Normal..AdCa.,
                                  Size = sapply(str_split(anno$Size..mm.,"x",simplify = F),"[",1),
                                  Location = anno$Location,
                                  NeoCellsInTumor = anno$NeoCellsInTumor,
                                  TumorInTotal = anno$TumorInTotal,
                                  stringsAsFactors = F)

#If Tissue is Normal, remove Size information
samples.tissue.types$Size[samples.tissue.types$Stage=="Normal"] = "" 

##Determine order of columns of heatmap
#get column order
onc.plot.col.order = colnames(onc.matrix.less[,column_order(onc.plot)[[2]]]) 
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
names(col_stage) = c("","Normal","Polyp",  "AdCa") #unique(samples.tissue.types$stage)

#Make color vector for Location
piratepal(palette = "southpark",plot.result = T)
col_location = piratepal(palette = "southpark",length.out = length(unique(samples.tissue.types$Location)))
names(col_location) = c("","Ascending","Transverse", "Descending", "Rectum")

#Make Continuous Colors for Neoplastic Cells in Tumor
col_neocellsintumor = colorRamp2(c(min(samples.tissue.types$NeoCellsInTumor,na.rm = T),
                                   max(samples.tissue.types$NeoCellsInTumor,na.rm = T)),
                                 c("orange", "darkred"))

#Make Continuous Colors for Tumor in Total
col_TumorInTotal = colorRamp2(c(min(samples.tissue.types$TumorInTotal,na.rm = T),
                                max(samples.tissue.types$TumorInTotal,na.rm = T)),
                              c("lightblue", "black"))


ha = HeatmapAnnotation(which = "column",show_annotation_name = T,
                       df = samples.tissue.types[,2:length(colnames(samples.tissue.types))],
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
onc.plot = draw(oncoPrint(mat = onc.matrix.less,get_type = get_type_fun,
                          alter_fun = alter_fun,
                          col = col,show_column_names = T,column_order = NULL,
                          pct_gp = gpar(cex = 0.8),column_names_gp = gpar(cex=0.6),row_names_gp = gpar(cex = 0.8),
                          bottom_annotation = ha,
                          width = unit(x = 12,units = "cm")))

#Save the plot
svg(filename = "A001_A002_F_G_oncoplot.svg",width = 12)
onc.plot
dev.off()


# 
# ### Adding colors to the annotations
# col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
# ha = HeatmapAnnotation(foo = 1:10, col = list(foo = col_fun))
# 
# 
# #Make the annotations now that we have the order of the columns
# onc.annos = columnAnnotation(samples.tissue.types)



# alter_fun = list(
#   mut1 = function(x, y, w, h) 
#     grid.rect(x, y, w, h, gp = gpar(fill = "red", col = NA)),
#   mut2 = function(x, y, w, h) 
#     grid.rect(x, y, w, h, gp = gpar(fill = "blue", col = NA)),
#   mut3 = function(x, y, w, h) 
#     grid.rect(x, y, w, h, gp = gpar(fill = "yellow", col = NA)),
#   mut4 = function(x, y, w, h) 
#     grid.rect(x, y, w, h, gp = gpar(fill = "purple", col = NA)),
#   mut5 = function(x, y, w, h) 
#     grid.rect(x, y, w, h, gp = gpar(fill = NA, lwd = 2)),
#   mut6 = function(x, y, w, h) 
#     grid.points(x, y, pch = 16),
#   mut7 = function(x, y, w, h) 
#     grid.segments(x - w*0.5, y - h*0.5, x + w*0.5, y + h*0.5, gp = gpar(lwd = 2))
# )
# test_alter_fun(alter_fun)
# 
# 
# 
