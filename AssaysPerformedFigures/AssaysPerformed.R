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
library(scales);library(ggsci)

setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/SampleInformation/")

#publish to web just the first sheet with first sheet in CSV format
#roxanne's sample sheet is below
source("~/aarons_FAP_github_repository/recent_Annotation.R")
# recent = recent_Annotation(go_online_clinical = "yes",go_online_pathology = "yes")
recent = recent_Annotation("no")
colnames(recent)
head(recent)
# #Just for RNA plan purposes. Needed to remove these because they weren't done yet
# recent$RNAseq[recent$RNAseq=="Batch 3: Jan 2020"] = ""
assay.list = c("DNA_WGS","DNA_WES","WG_Bisulfite","RNAseq","bulk.ATAC","snRNAseq",
               "snATACseq","scRNAseq","HE.CODEX","Proteomics","Lipidomics","Metabolomics")

#replace assay specific details with word "done" in "recent"
for (assay in assay.list) {
  assay.bool = recent[,assay]!="" #find all opposites of empty spaces 
  recent[,assay][assay.bool] = "done" #within each assay column, replace assay specifics with word "done"
  recent[,assay][!assay.bool] = "not done"
}

head(recent[,assay.list])

#select specifc columns
#cols: SampleName, Physical state of tissue, Assays of interest
assaydone=recent[,c("SampleName",assay.list)]
head(assaydone)

#choose specific samples to look at
sample.choices = paste0("A001|A002|EP|JP|B001")
# sample.choices = paste0("B00")
samples.of.interest = which(grepl(pattern = sample.choices,x = assaydone$SampleName))
assaydone.pts = assaydone[samples.of.interest,]
head(assaydone.pts)

# Merge Samples names and all Assays done on them
# Takes care of duplicate rows for when the same sample is preserved different ways
samples = unique(assaydone.pts$SampleName)
m = matrix(data = "0",nrow = length(samples),ncol = length(colnames(assaydone.pts)))
colnames(m) = colnames(assaydone.pts)
rownames(m) = samples

i = "B001-A-001"
i = "A001-C-004"
for (i in samples) {
  sample.done.notdone = assaydone.pts[assaydone.pts$SampleName==i,-1]
  bool = apply(sample.done.notdone,2,function(x) !grepl("not done",x) & grepl("done",x))
  if (!is.vector(bool)) {
    m[i,] = c("SampleName" = i,colSums(bool))
  } else {
    bool.0.1 = ifelse(bool,1,0)
    bool.0.1 = c("SampleName" = i,bool.0.1)
    m[i,] = bool.0.1
  }
}

m[m=="1"] = "done"
m[m=="0"] = "not done"
assaydone.pts = data.frame(m)


#Arrange rows by DNA_WGS, change rownames to sample names, 
# Remove unnecessary columns and transpose
assaydone.pts = arrange(assaydone.pts,DNA_WES)
assaydone.pts = arrange(assaydone.pts,DNA_WGS)
rownames(assaydone.pts) = assaydone.pts$SampleName
# assaydone.pts = select(assaydone.pts,-Physical.state.of.tissue, -SampleName)
head(assaydone.pts)


#Transpose for matrix formation
assaydone.pts.t = t(as.matrix(assaydone.pts))
head(assaydone.pts.t)

# #remove "not done" rows and cols
assays.done = apply(assaydone.pts.t, 1, function(r) any(r %in% c("done"))) #check rows for presence of "done"
assaydone.pts.t = assaydone.pts.t[assays.done,]
#check cols for presence of "done"
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

# only for tesing
# assaydone.pts.t = assaydone.pts.t[-1,sort(colnames(assaydone.pts.t))]

#make heatmap. Rough Draft Heatmap
fig.heatmap = ComplexHeatmap::Heatmap(matrix = assaydone.pts.t,col = c("done"="red","not done" = "grey"),
                        cluster_rows = F,cluster_columns = F,
                        show_row_names = T,show_column_names = T,
                        column_names_gp = gpar(cex = 1),width = 4,
                        rect_gp = gpar(col = "white"))
fig.heatmap
# getwd()
pdf(file = "~/Desktop/AssaysPerformed.pdf",width = 10,height = 5)
fig.heatmap
dev.off()
# svg(filename = "AssaysPerformed.svg",width = 8)
# fig.heatmap
# dev.off()

########### Add Clinical and other Annotations
#Change the order of the rows to match that of the heatmap
clinical.anno = recent[,c("SampleName","Patient","Stage..Polyp..Normal..AdCa.",
                          "Location","Dysplasia", "NeoCellsInTumor","nonNeoStromaInTotal","TumorInTotal",
                          "PercentNormalOverall","PercentStromaOverall","PercentCancerOverall",
                          "PercentAdenomaOverall","PercentStromaInCancer","PercentStromaInAdenoma",
                          "PercentNecrosisInTumor","CurrentState")]
#Only collect clinical info on samples which have assays done
clinical.anno = clinical.anno[clinical.anno$SampleName %in% colnames(assaydone.pts.t),]
#place in same order at heatmap
clinical.anno = clinical.anno[match(colnames(assaydone.pts.t),clinical.anno$SampleName),] 

#rename some column names
clinical.anno = dplyr::select(clinical.anno, SampleName,
                       Patient,
                       Stage = "Stage..Polyp..Normal..AdCa.",
                       Location,Dysplasia,
                       NeoplasticInTumor = "NeoCellsInTumor",
                       TumorInTotal,
                       PercentNormalOverall,
                       PercentStromaOverall,
                       PercentCancerOverall,
                       PercentStromaInCancer,
                       PercentNecrosisInTumor,
                       PercentAdenomaOverall,
                       PercentStromaInAdenoma,
                       StateOfSamples = CurrentState)
str(clinical.anno)
summary(clinical.anno)
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
clinical.anno$Stage = factor(clinical.anno$Stage,levels = c("Blood","Normal", "Polyp",  "AdCa"))
col_stage = brewer.pal(n = length(unique(clinical.anno$Stage)),name = "Reds")
names(col_stage) = c("Blood","Normal","Polyp",  "AdCa")


#Make Color Vector for Location
unique(clinical.anno$Location)
col_location = piratepal(palette = "info2",
                         length.out = length(unique(clinical.anno$Location)))
if (bool.B001) {
  names(col_location) = c("","Ascending","Transverse", "Descending", "Rectum","Sigmoid","Ileum","Jejunum","Mid-jejunum", "Duodenum")
} else {
  names(col_location) = c("","Ascending","Transverse", "Descending", "Rectum")
}

###
# Make color Vector for Dysplasia: Yes and no
unique(clinical.anno$Dysplasia)
clinical.anno$Dysplasia[is.na(clinical.anno$Dysplasia)] = ""
col_dysplasia = c("white", "rosybrown1", "lightgrey")
names(col_dysplasia) = unique(clinical.anno$Dysplasia)

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
# "PercentNecrosisInTumor" 
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

#Make Color Vector for State of tissue
unique(clinical.anno$StateOfSamples)
table(clinical.anno$StateOfSamples)
x = clinical.anno$StateOfSamples
x[x!="Depleted" & x!="Unfound" & x!= ""] = "Available"
table(x)
clinical.anno$StateOfSamples = x
col_state = piratepal(palette = "southpark",length.out = length(unique(clinical.anno$StateOfSamples)))
names(col_state) = unique(clinical.anno$StateOfSamples)


#make an object for the annotations
colnames(clinical.anno)
# anno_df = data.frame(clinical.anno[,2:length(colnames(clinical.anno))])
anno_df = data.frame(clinical.anno[,c(2:7,15)])
colnames(anno_df)

#With less Annotations
ha = HeatmapAnnotation(df = anno_df,show_annotation_name = T,
                       col = list(Patient = col_patient,
                                  Stage = col_stage,
                                  Location = col_location,
                                  Dysplasia = col_dysplasia,
                                  NeoplasticInTumor = col_neocellsintumor,
                                  TumorInTotal = col_TumorInTotal,
                                  StateOfSamples = col_state),
                       na_col = "white", gp = gpar(col = "black"))

# ha = HeatmapAnnotation(df = anno_df,show_annotation_name = T,
#                        col = list(Patient = col_patient,
#                                   Stage = col_stage,
#                                   Location = col_location,
#                                   NeoplasticInTumor = col_neocellsintumor,
#                                   TumorInTotal = col_TumorInTotal,
#                                   PercentNormalOverall = col_PercentNormalOverall,
#                                   PercentStromaOverall = col_PercentStromaOverall,
#                                   PercentCancerOverall = col_PercentCancerOverall,
#                                   PercentStromaInCancer = col_PercentStromaInCancer,
#                                   PercentNecrosisInTumor = col_PercentNecrosisInTumor,
#                                   StateOfSamples = col_state),
#                        na_col = "white", gp = gpar(col = "black"))


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

# Text Vector of number of assays completed
num.assays.completed = apply(assaydone.pts.t,1,function(x) sum(x=="done"))

#Make and print the Heatmap
assaysdone.heatmap = oncoPrint(mat = assaydone.pts.t,
                               get_type = get_type_fun,alter_fun = alter_fun_list,
                               col = col,column_order = colnames(assaydone.pts.t),
          show_column_names = T,column_names_gp = gpar(cex = 0.8),
          bottom_annotation = ha,show_pct = F,
          row_title = "Assays Performed",row_title_side = "right",
          heatmap_legend_param = list(title = "Assay Performed"),
          right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(),
                                           values = anno_text(num.assays.completed)))

assaysdone.heatmap
getwd()

pdf(file = paste0("AssaysPerformed_annotated",format(Sys.time(),"%b%d%Y"),
                  "_annotations.pdf"),
    width = 14, height = 10)
assaysdone.heatmap
dev.off()

pdf(file = paste0("AssaysPerformed_annotated",format(Sys.time(),"%b%d%Y"),
                  "_heatmap.pdf"),
    width = 14, height = 10)

dev.off()

########### Assays per patient ##########
data = t(assaydone.pts.t)
rownames(data)
pts = c("A001","A002","F","G","B001")

#create matrix for per patient number of assays done
m = matrix(0,nrow = length(colnames(data)),ncol = length(pts))
colnames(m) = pts
rownames(m) = colnames(data)

for (pt in pts) {
  data.reduced = data[grep(paste0("^",pt),rownames(data)),]
  assay.nums = apply(data.reduced,2,function(x) sum(grepl("done",x)))
  m[,pt] = assay.nums
}

m

assays.to.include = c("DNA_WGS","RNAseq", "HE.CODEX", "Proteomics", "Lipidomics","snRNAseq", "snATACseq","scRNAseq")
m = m[assays.to.include,]

assay.colors = c(pal_tron("legacy",alpha = 1)(7),pal_npg()(dim(m)[1]-7))
names(assay.colors) = rownames(m)
assay.colors["Lipidomics"] = "mediumpurple3"


# Print out figure
pdf("NumberAssaysperPatient.pdf",width = 15,height = 7)
layout.matrix = matrix(c(1,1,1,2),nrow = 1)
layout(layout.matrix)
bp = barplot(m,beside = T,main = "Number of Assays done per Patient",
             col = assay.colors[rownames(m)],cex.axis = 2,cex.names = 2)
text(x = bp,y = as.numeric(m)+1,
     labels = str_remove_all(as.numeric(m),"^0"),
     cex = 1, srt = 30)
plot(1,1,"n",axes = F,xlab = "", ylab = "")
legend("left",legend = rownames(m),fill = assay.colors[rownames(m)])
dev.off()

###### Assays on X axis and Patient on Y
pdf("AssaysPerformed_perPt_barplot.pdf",height = 5,width = 10)
layout(1)
par(mar = c(6, 4.1, 3.5, 2.1))
pt.assays = t(m)
pt.assays.sum = pt.assays
for (row in 2:nrow(pt.assays)) {
  pt.assays.sum[row,] = pt.assays.sum[row-1,]+pt.assays[row,]
}

x = barplot(pt.assays,legend.text = T,las = 2,col = col_patient,
        main = "Total Assays performed for Every Patient")
for (assay in 1:ncol(pt.assays)) {
  text(x = x[assay],y = pt.assays.sum[,assay]-2,
       labels = str_replace(pt.assays[,assay],"0",""),xpd = NA)
}

barheight = pt.assays.sum[nrow(pt.assays.sum),]
text(x,barheight+3,barheight,xpd = NA,font = 2)

dev.off()


# Write Table
getwd()
write.table(assaydone.pts.t,file = paste0("assaysdone_matrix_",format(Sys.time(),"%b%d%Y"),".txt"),append = F,sep = "\t")
