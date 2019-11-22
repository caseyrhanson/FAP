rm(list=ls())

plot.related.colors = function(color){
  # plot.related.colors("green")
  x = colors()[grep(color,colors())]
  plot(1:length(x),1:length(x),col = x,pch = 15)
  text(1:length(x),1:length(x),labels = x,pos = 4,cex = 0.6)
  return(x)
}


library(dplyr);library(tidyr);library(stringr)
library(yarrr);library(RColorBrewer);library(circlize);library(ComplexHeatmap)
setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/SampleInformation/")

#publish to web just the first sheet with first sheet in CSV format
#roxanne's sample sheet is below
source("~/aarons_FAP_github_repository/recent_Annotation.R")
# recent = recent_Annotation(go_online_clinical = "yes",go_online_pathology = "yes")
recent = recent_Annotation("no")
colnames(recent)
head(recent)

assay.list = c("DNA_WGS","DNA_WES","WG_Bisulfite","RNAseq","bulk.ATAC","snRNAseq",
               "snATACseq","scRNAseq","HE.CODEX","Proteomics","Lipidomics","Metabolomics")
samplenames = colnames(recent)[c(1,8)] #to ensure duplicate samples are not included

#replace assay specific details with word "done" in "recent"
for (assay in assay.list) {
  assay.bool = recent[,assay]!="" #find all opposites of empty spaces 
  recent[,assay][assay.bool] = "done" #within each assay column, replace assay specifics with word "done"
  recent[,assay][!assay.bool] = "not done"
}

head(recent[,assay.list])

#select specifc columns
assaydone=recent[,c(samplenames,assay.list)]
head(assaydone)


#make all sample names de-duplicated
sample.rename.index = which(duplicated(assaydone$SampleName))
assaydone$SampleName[sample.rename.index] = str_c(assaydone$SampleName[sample.rename.index],assaydone$Physical.state.of.tissue[sample.rename.index],sep = "_")
head(assaydone)

#choose specific samples to look at
sample.choices = paste0("A001|A002|EP|JP|B001")
# sample.choices = paste0("B001")
samples.of.interest = which(grepl(pattern = sample.choices,x = assaydone$SampleName))
assaydone.pts = assaydone[samples.of.interest,]
head(assaydone.pts)


#Arrange rows by DNA_WGS, change rownames to sample names, Remove unnecessary columns and transpose
assaydone.pts = arrange(assaydone.pts,DNA_WES)
assaydone.pts = arrange(assaydone.pts,DNA_WGS)
rownames(assaydone.pts) = assaydone.pts$SampleName
assaydone.pts = select(assaydone.pts,-Physical.state.of.tissue, -SampleName)
head(assaydone.pts)


#Transpose for matrix formation
assaydone.pts.t = t(as.matrix(assaydone.pts))
head(assaydone.pts.t)

# #remove "not done" rows and cols
assays.done = apply(assaydone.pts.t, 1, function(r) any(r %in% c("done")))
assaydone.pts.t = assaydone.pts.t[assays.done,]

samples.done = apply(assaydone.pts.t, 2, function(cols) any(cols %in% c("done"))) 
assaydone.pts.t = assaydone.pts.t[,samples.done]
head(assaydone.pts.t)

#Arrange the columns so that A001 samples comes first and JP comes last
assaydone.pts.t = assaydone.pts.t[,c(grep("B001",colnames(assaydone.pts.t),value = T),
                                grep("A001",colnames(assaydone.pts.t),value = T),
                                grep("A002", colnames(assaydone.pts.t),value = T),
                                grep("EP",colnames(assaydone.pts.t),value = T),
                                grep("JP",colnames(assaydone.pts.t),value = T))]
head(assaydone.pts.t)

#make heatmap. Rough Draft Heatmap
fig.heatmap = ComplexHeatmap::Heatmap(matrix = assaydone.pts.t,col = c("done"="red","not done" = "grey"),
                        cluster_rows = F,cluster_columns = F,
                        show_row_names = T,show_column_names = T,
                        column_names_gp = gpar(cex = 1),width = 4,
                        rect_gp = gpar(col = "white"))
fig.heatmap
# getwd()
# svg(filename = "AssaysPerformed.svg",width = 8)
# fig.heatmap
# dev.off()

########### Add Clinical and other Annotations
#Change the order of the rows to match that of the heatmap
clinical.anno = recent[,c("SampleName","Patient","Stage..Polyp..Normal..AdCa.","Location",
                          "NeoCellsInTumor","nonNeoStromaInTotal","TumorInTotal",
                          "PercentNormalOverall","PercentStromaOverall","PercentCancerOverall",
                          "PercentAdenomaOverall","PercentStromaInCancer","PercentStromaInAdenoma",
                          "PercentNecrosisInTumor")]
#Only collect clinical info on samples which have assays done
clinical.anno = clinical.anno[clinical.anno$SampleName %in% colnames(assaydone.pts.t),]
#place in same order at heatmap
clinical.anno = clinical.anno[match(colnames(assaydone.pts.t),clinical.anno$SampleName),] 

#rename some column names
clinical.anno = select(clinical.anno, SampleName,
                       Patient,
                       Stage = "Stage..Polyp..Normal..AdCa.",
                       Location,
                       NeoplasticInTumor = "NeoCellsInTumor",
                       TumorInTotal,
                       PercentNormalOverall,
                       PercentStromaOverall,
                       PercentCancerOverall,
                       PercentStromaInCancer,
                       PercentNecrosisInTumor,
                       PercentAdenomaOverall,
                       PercentStromaInAdenoma)
str(clinical.anno)

# Some figures don't want to include the Normal B001 sample. But some do.
# This bool is true when the annotation table has B001 in it.
bool.B001 = any(grepl(pattern = "B001",x = clinical.anno$Patient))

#Make Color Vector for Patients. 
#Change the names of the patients: EP to F and JP to G
clinical.anno$Patient = str_replace_all(clinical.anno$Patient,pattern = "EP",replacement = "F")
clinical.anno$Patient = str_replace_all(clinical.anno$Patient,pattern = "JP",replacement = "G")
col_patient = brewer.pal(n = length(unique(clinical.anno$Patient)), name = "Accent")
names(col_patient) = unique(clinical.anno$Patient)

#Make Color Vector for Stage
unique(clinical.anno$Stage) 
col_stage = brewer.pal(n = length(unique(clinical.anno$Stage)),name = "Reds")

if(bool.B001){
   names(col_stage) = c("","Normal (no FAP)","Normal","Polyp",  "AdCa")
 } else {
   names(col_stage) = c("","Normal","Polyp",  "AdCa")
 }
#unique(samples.tissue.types$stage)
#unique(samples.tissue.types$stage)

#Make Color Vector for Location
unique(clinical.anno$Location)
col_location = piratepal(palette = "info2",
                         length.out = length(unique(clinical.anno$Location)))
if (bool.B001) {
  names(col_location) = c("","Ascending","Transverse", "Descending", "Rectum","Sigmoid","Ileum","Jejunum","Mid Jejunum", "Duodenum")
} else {
  names(col_location) = c("","Ascending","Transverse", "Descending", "Rectum")
}

#Make Continuous Colors for Neoplastic Cells in Tumor
col_neocellsintumor = circlize::colorRamp2(c(min(recent$NeoCellsInTumor,na.rm = T),
                        max(recent$NeoCellsInTumor,na.rm = T)),
                      c("orange", "darkred"))

#Make Continuous Colors for Tumor in Total
col_TumorInTotal = circlize::colorRamp2(c(min(recent$TumorInTotal,na.rm = T),
                                   max(recent$TumorInTotal,na.rm = T)),
                                 c("lightblue", "black"))


# Make Continuous Colors for PercentNormalOverall
# PercentNormalOverall light greeen dark green
col_PercentNormalOverall = circlize::colorRamp2(c(min(recent$PercentNormalOverall,na.rm = T),
                                          max(recent$PercentNormalOverall,na.rm = T)),
                                        c("lightgreen", "darkgreen"))

# Make Continuous Colors for Percent Stroma Overall
# "PercentStromaOverall"  
col_PercentStromaOverall = circlize::colorRamp2(c(min(recent$PercentStromaOverall,na.rm = T),
                                                  max(recent$PercentStromaOverall,na.rm = T)),
                                                c("lightyellow1", "yellowgreen"))


# Make Continuous Colors for Percent Cancer Overall
# "PercentCancerOverall"   
col_PercentCancerOverall = circlize::colorRamp2(c(min(recent$PercentCancerOverall,na.rm = T),
                                                  max(recent$PercentCancerOverall,na.rm = T)),
                                                c("mediumpurple1", "mediumpurple4"))

# Make Continuous Colors for Percent Stroma in Cancer
# "PercentStromaInCancer"  
col_PercentStromaInCancer = circlize::colorRamp2(c(min(recent$PercentStromaInCancer,na.rm = T),
                                                  max(recent$PercentStromaInCancer,na.rm = T)),
                                                c("mediumpurple4","yellowgreen"))

# Make Continuous Colors for Percent Necrosis in Tumor
"PercentNecrosisInTumor" 
col_PercentNecrosisInTumor = circlize::colorRamp2(c(min(recent$PercentNecrosisInTumor,na.rm = T),
                                                   1),
                                                 c("grey89","red4"))


# Make Continuous Colors for Percent Adenoma Overall
# "PercentAdenomaOverall" 
col_PercentAdenomaOverall = circlize::colorRamp2(c(min(recent$PercentAdenomaOverall,na.rm = T),
                                                   max(recent$PercentAdenomaOverall,na.rm = T)),
                                                 c("lightpink","deeppink3"))

# Make Continuous Colors for PercentStromaInAdenoma
# "PercentStromaInAdenoma"
col_PercentStromaInAdenoma = circlize::colorRamp2(c(min(recent$PercentStromaInAdenoma,na.rm = T),
                                                   max(recent$PercentStromaInAdenoma,na.rm = T)),
                                                 c("yellowgreen","deeppink3"))



#make an object for the annotations
colnames(clinical.anno)
anno_df = data.frame(clinical.anno[,2:11])
colnames(anno_df)

ha = HeatmapAnnotation(df = anno_df,show_annotation_name = T,col = list(Patient = col_patient,
                                       Stage = col_stage,
                                       Location = col_location,
                                       NeoplasticInTumor = col_neocellsintumor,
                                       TumorInTotal = col_TumorInTotal,
                                       PercentNormalOverall = col_PercentNormalOverall,
                                       PercentStromaOverall = col_PercentStromaOverall,
                                       PercentCancerOverall = col_PercentCancerOverall,
                                       PercentStromaInCancer = col_PercentStromaInCancer,
                                       PercentNecrosisInTumor = col_PercentNecrosisInTumor),
                       na_col = "white", gp = gpar(col = "black"))


# ha = HeatmapAnnotation(which = "column",show_annotation_name = T,
#                        df = anno_df,
#                        col = list(Patient = col_patient,
#                                   Stage = col_stage,
#                                   Location = col_location,
#                                   NeoplasticInTumor = col_neocellsintumor, 
#                                   TumorInTotal = col_TumorInTotal),
#                        na_col = "white",
#                        gp = gpar(col = "black"))


# Make oncoplot with Heatmap data. Oncoplot is special type of ComplexHeatmap
assaydone.pts.t[assaydone.pts.t=="not done"] = "" #remove "not done" from matrix. Cleans up the matrix
get_type_fun = function(x) strsplit(x, ";")[[1]]
col = c(done = "red", notdone = "grey")
alter_fun_list = list(background = function(x,y,w,h) {grid.rect(x,y,w,h, gp = gpar(fill = "grey", col = "white"))},
                      done = function(x, y, w, h) {grid.rect(x, y, w, h, gp = gpar(fill = col["done"], col = NA))},
                      notdone = function(x, y, w, h) {grid.rect(x, y, w, h, gp = gpar(fill = col["notdone"], col = NA))})
#change EP to F and JP to G
colnames(assaydone.pts.t) = str_replace_all(colnames(assaydone.pts.t),pattern = "EP",replacement = "F")
colnames(assaydone.pts.t) = str_replace_all(colnames(assaydone.pts.t),pattern = "JP",replacement = "G")

#Make and print the Heatmap
assaysdone.heatmap = oncoPrint(mat = assaydone.pts.t,get_type = get_type_fun,alter_fun = alter_fun_list,col = col,column_order = NULL,
          show_column_names = T,column_names_gp = gpar(cex = 0.8),bottom_annotation = ha,
          row_title = "Assays Performed",row_title_side = "right",
          heatmap_legend_param = list(title = "Assay Performed"))

assaysdone.heatmap
getwd()

svg(filename = paste0("AssaysPerformed_annotated",format(Sys.time(),"%b%d%Y"),"_annotations.svg"),width = 14, height = 20)
assaysdone.heatmap
dev.off()

svg(filename = paste0("AssaysPerformed_annotated",format(Sys.time(),"%b%d%Y"),"_heatmap.svg"),width = 14, height = 10)
assaysdone.heatmap
dev.off()

