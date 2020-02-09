rm(list=ls())

library(dplyr);library(tidyr);library(stringr);library(yarrr);library(RColorBrewer);library(ggrepel)
setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/SampleInformation/")

#publish to web just the first sheet with first sheet in CSV format
#roxanne's sample sheet is below
source("~/aarons_FAP_github_repository/recent_Annotation.R")
recent = recent_Annotation(go_online_clinical = "no")
colnames(recent)

recent = recent%>%select(SampleName,Location,Patient, 
                Size = "Size..mm.",
                "NeoCellsInTumor",
                "nonNeoStromaInTotal",
                "TumorInTotal")
recent$Patient = str_replace_all(string = recent$Patient,pattern = "EP",replacement = "F")
recent$Patient = str_replace_all(string = recent$Patient,pattern = "JP",replacement = "G")

sum(!is.na(recent$NeoCellsInTumor))
#remove NA pathology rows
recent = filter(recent,NeoCellsInTumor!="")
#Use only first value of size and make it numeric
recent$Size = sapply(str_split(recent$Size,"x",simplify = F),"[",1)
recent$Size = as.numeric(recent$Size)

str(recent$Size)
recent$NeoCellsInTumor = as.numeric(recent$NeoCellsInTumor)
str(recent$NeoCellsInTumor)



  
  
# Patient Colors
col_patient = yarrr::piratepal(palette = "basel",length.out = length(unique(recent$Patient)))
names(col_patient) = unique(recent$Patient)

svg(filename = "size_v_NeoCellsInTumor_names.svg")
cor(x = recent$Size,y = recent$NeoCellsInTumor, use = "complete.obs")
plot(x = recent$Size,y = recent$NeoCellsInTumor,bg = col_patient[recent$Patient],
     pch = 21, cex = 1.5,ylab = "Percent (%) Neoplastic Cells in Tumor Area",xlab = "Size of Lesion (mm)")
text(x = recent$Size,y = recent$NeoCellsInTumor + sample(-20:20,1),labels = recent$SampleName,srt = 45,xpd = NA)
legend("topright",legend = names(col_patient)[1:4],pt.bg = col_patient,pch = 21)
summary(lm(recent$Size ~ recent$NeoCellsInTumor))
dev.off()

recent %>%
ggplot(aes(x = Size, y = NeoCellsInTumor,col = Patient)) +
  geom_point() +
  geom_text_repel(aes(label=SampleName))
ggsave("size_v_NeoCellsInTumor_names.pdf")

svg(filename = "size_v_TumorInTotal_names.svg")
cor(x = recent$Size,y = recent$TumorInTotal, use = "complete.obs")
plot(x = recent$Size,y = recent$TumorInTotal,bg = col_patient[recent$Patient],
     pch = 21, cex = 1.5,ylab = "Percent (%) Tumor in Total Area", xlab = "Size of Lesion (mm)")
text(x = recent$Size,y = recent$TumorInTotal, labels = recent$SampleName,srt = 90,xpd = NA)
legend("topright",legend = names(col_patient)[1:4],pt.bg = col_patient,pch = 21)
summary(lm(recent$Size ~ recent$TumorInTotal))
dev.off()

recent %>%
  ggplot(aes(x = Size, y = TumorInTotal,col = Patient)) +
  geom_point() +
  geom_text_repel(aes(label=SampleName))
ggsave("size_v_TumorInTotal_names.pdf")


svg(filename = "size_v_nonNeoStromaInTotal_names.svg")
cor(x = recent$Size,y = recent$nonNeoStromaInTotal, use = "complete.obs")
plot(x = recent$Size,y = recent$nonNeoStromaInTotal,bg = col_patient[recent$Patient],
     pch = 21, cex = 1.5,ylab = "Percent (%) Non-Neoplastic Stromal Cells in Total Area", xlab = "Size of Lesion (mm)")
text(x = recent$Size,y = recent$nonNeoStromaInTotal,labels = recent$SampleName,srt = 80, xpd = NA)
legend("topright",legend = names(col_patient)[1:4],pt.bg = col_patient,pch = 21)
summary(lm(recent$Size ~ recent$nonNeoStromaInTotal))
dev.off()
getwd()


recent %>%
  ggplot(aes(x = Size, y = nonNeoStromaInTotal,col = Patient)) +
  geom_point() +
  geom_text_repel(aes(label=SampleName))
ggsave("size_v_nonNeoStromaInTotal_names.pdf")
