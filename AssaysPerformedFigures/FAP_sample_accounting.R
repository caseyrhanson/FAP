rm(list=ls())

library(dplyr);library(tidyr);library(GenomicRanges);library(reshape2);library(stringr)

setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/SampleInformation/")
# 
# x = tbl_df(read.table(file = "A001_A002_2018NovDec.txt",header = T,sep = "\t"))
# colnames(x)
# head(x)
#publish to web just the first sheet with first sheet in CSV format
source("~/aarons_FAP_github_repository/recent_Annotation.R")
recent = recent_Annotation(go_online_clinical = "yes",go_online_pathology = "yes")
# recent = recent_Annotation("no")
# recent = read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRd0tND0K2rcfq_C6pte4_W42vkuwsemkv_A_jnrm6mm5N72nxxw5Bezumy898XmtzEQlcxfrYW9j2E/pub?gid=0&single=true&output=csv")
colnames(recent)
summary(recent)
head(recent)

less = select(recent,"SampleName", PT = "Patient","Collection",Tissue_Type = "Stage..Polyp..Normal..AdCa.",
              Anatomical_Location = "Location","Size..mm.",
              "Notes", "VAP_Names")%>%separate(col = "Size..mm.",into = "Size..mm.",sep = "x")
#Add a character size column based on first dimension of 3 measurements of size
less$Size..mm. = as.numeric(as.character(less$Size..mm.))
less$Size = rep(NA,nrow(less))
less$Size[less$Size..mm.<5] = "S"
less$Size[less$Size..mm.>=5 & less$Size..mm. < 10] = "M"
less$Size[less$Size..mm.>=10] = "L"

# less = select(recent,"sampleName", "Sample_Type",PT="Patient",
#        "Sample_Number", "BigSize",Tissue_Type="Stage", Anatomical_Location="Location_of_Sampled_Specimen", "Size_.mm.",Size="Size_SML.S..5.5.M..10.L.10.",
#        Storage_Method="PreservationMethod", "Grade","Phenotype", "Weight","Collection", "PatientSample_label",
#        "blank", "Sections", "DNA_WGS","DNA_WES", "WG_Bisulfite", "RNAseq","ATAC", 
#        "single.cell.dissociation", "scRNA","scATAC")
x = less
colnames(less)
head(less)
less$Collection

unique(x$PT)
x$PT = str_replace(string = x$PT,pattern = "EP",replacement = "F")
x$PT = str_replace(string = x$PT,pattern = "JP",replacement = "G")
x$PT = str_replace(string = x$PT,pattern = "JC",replacement = "E")
x$PT = str_replace(string = x$PT,pattern = "NBM",replacement = "C")
x$PT = str_replace(string = x$PT,pattern = "ES",replacement = "D")
unique(x$Collection)
unique(x$Tissue_Type)
unique(x$Phenotype)

unique(x$Anatomical_Location)
id.Rect = grep(pattern = "rec",ignore.case = T, x = x$Anatomical_Location) #if the Anatomical location has the word "Rect" just rename it "Rectum"
x$Anatomical_Location[id.Rect]="Rectum"
id.Sig = grep(pattern = "Sig",x = x$Anatomical_Location) #if the Anatomical location has the word "Sig" just rename it "Sigmoid"
x$Anatomical_Location[id.Sig]="Sigmoid"
id.Desc = grep(pattern = "desc",x = x$Anatomical_Location) #if the Anatomical location has the word "Sig" just rename it "Sigmoid"
x$Anatomical_Location[id.Desc]="DESC"
id.Desc = grep(pattern = "ecum",x = x$Anatomical_Location) #if the Anatomical location has the word "Sig" just rename it "Sigmoid"
x$Anatomical_Location[id.Desc]="Cecum"
id.Desc = grep(pattern = "junum",x = x$Anatomical_Location) #if the Anatomical location has the word "Sig" just rename it "Sigmoid"
x$Anatomical_Location[id.Desc]="Jejunum"

unique(x$Anatomical_Location)
unique(x$Size)
unique(x$Storage_Method) #may need to adjust this column. maybe make a new one

#make a table with these columns:
c("PT", "Collection", "Tissue_Type", "Phenotype", "Anatomical_Location", "Size", "Storage_Method")

unique(x$PT)
x$PT=factor(x$PT,levels=c("F", "G", "C",  "D", "E", "A001", "A002", "A003","A005", "A006", "A008","A009","A010","A011","A012","B001","B002"))
# x$PT=factor(x$PT,levels=c("EP", "JP", "NBM",  "ES", "JC", "A001", "A002", "A003","B001","A005", "A006", "A008"))

x$Anatomical_Location = factor(x$Anatomical_Location, levels = c("Duodenum","Jejunum","Ileum",
                                                                 "Cecum","Ascending","Transverse",
                                                                 "Descending","Sigmoid", "Rectum",
                                                                 "", "Heart"))

x$Collection= factor(x$Collection, levels = c("Endoscopy", "Colectomy","Post-Operation Check Up","Autopsy"))
x$Size = factor(x$Size, levels = c("S","M","L","NA"))
x$Tissue_Type = factor(x$Tissue_Type, levels = c("Normal","unsure adenoma","Polyp", "AdCa",  "", "Normal (no FAP)"))
x$Phenotype = factor(x$Phenotype, levels = c("Normal","ND","Sessile","Stalk"))

sum(is.na(x$Size))

x = select(x,PT, Tissue_Type, Anatomical_Location,Size, Collection)
sum = summary(x, maxsum = 15); print(sum)

write.table(x = sum,file =  paste0(format(Sys.Date(),"%b%d%Y"),"_summary_uniq.txt"),
            append = F,sep = "\t",
            col.names = T,row.names = F,quote = T)

# x_grouped = group_by(x,PT,Collection,Anatomical_Location,Size, Tissue_Type)
# x_grouped = group_by(x,PT,Collection,Anatomical_Location,Size, Tissue_Type,Phenotype)
# 
# y = summarise(x_grouped, n = n())
# 
# print(y,n = nrow(y))
# 
# write.table(x = y,file = paste0(format(Sys.Date(),"%b%d%Y"),"_summary_uniq.txt"),append = F,sep = "\t",row.names = F)
            
