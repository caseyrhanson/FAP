rm(list=ls())

library(dplyr);library(tidyr);library(stringr)

pt = "A001"
jaccard = read.csv(file = "~/Bulk_A001_A002/treeomics/A001/A001_12_jsc-matrix.csv",
         header = T,sep = ",",skip = 1)
str(jaccard)
jaccard

jaccard = jaccard%>%gather(key = SampleID,value = Jac,-Sample)
jaccard = jaccard[grep(pattern = "blood",x = jaccard$Sample,invert = T),]
jaccard = jaccard[grep(pattern = "blood",x = jaccard$SampleID,invert = T),]

#Remove A001C and 00 from the sample names
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "A001C")
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "00")
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "^0")



#Remove A001C and 00 from the sample names
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "A001C")
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "00")
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "^0")

jaccard = jaccard%>%unite(col = "sample_to_sample",c("Sample","SampleID"),sep = "_")

write.table(jaccard,file = "~/Bulk_A001_A002/schematic_colons/A001jaccard.txt",append = F,sep = "\t",row.names = F)

head(jaccard)
## read in the physical distance data
euc = read.table("A001_ColonSchematic_Polypslocations_distances.txt",header = T)
head(euc)
euc = euc%>%unite(col = "sample_to_sample",c("Sample","SampleID"),sep = "_")

jac.euc 
x = left_join(jaccard,euc,"sample_to_sample")%>%filter(Jac!=1)%>%filter(!is.na(Euc_dist))
duplicated(x$sample_to_sample)

pdf(file = "A001jaccard_v_Euc.pdf",width = 8,height = 6)
plot(x$Jac,x$Euc_dist,
     ylab = "Euclidean distance (0 = self)",
     xlab = "Jaccard Statistics (1 = self)")
abline(lm(x$Euc_dist ~ x$Jac))
title(main = paste0(pt,": Is Euclidean Distance Correlated with Jaccard Similarity"))
test = cor.test(x$Jac,x$Euc_dist)
p = round(test$p.value,digits = 4)
r = round(test$estimate,digits = 2)
text(0.14,50,
     paste0("R is ",r,"and R^2 is ",round(r*r,digits = 2)," with p-value ",p), cex = 0.8)
dev.off()

############################
pt = "A002"
jaccard = read.csv(file = "~/Bulk_A001_A002/treeomics/A002/A002_13_jsc-matrix.csv",
                   header = T,sep = ",",skip = 1)
str(jaccard)

jaccard = jaccard%>%gather(key = SampleID,value = Jac,-Sample)
jaccard = jaccard[grep(pattern = "blood",x = jaccard$Sample,invert = T),]
jaccard = jaccard[grep(pattern = "blood",x = jaccard$SampleID,invert = T),]

#Remove A002C and 00 from the sample names
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "A002C")
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "00")
jaccard$Sample = str_remove_all(string = jaccard$Sample,pattern = "^0")



#Remove A002C and 00 from the sample names
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "A002C")
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "00")
jaccard$SampleID = str_remove_all(string = jaccard$SampleID,pattern = "^0")

jaccard = jaccard%>%unite(col = "sample_to_sample",c("Sample","SampleID"),sep = "_")


write.table(jaccard,file = "~/Bulk_A001_A002/schematic_colons/A002jaccard.txt",append = F,sep = "\t",row.names = F)

## read in the physical distance data
euc = read.table("A002_ColonSchematic_Polypslocations_distances.txt",header = T)
head(euc)
euc = euc%>%unite(col = "sample_to_sample",c("Sample","SampleID"),sep = "_")

jac.euc 
x = left_join(jaccard,euc,"sample_to_sample")%>%filter(Jac!=1)%>%filter(!is.na(Euc_dist))
duplicated(x$sample_to_sample)

pdf(file = "A002jaccard_v_Euc.pdf",width = 8,height = 6)
plot(x$Jac,x$Euc_dist,
     ylab = "Euclidean distance (0 = self)",
     xlab = "Jaccard Statistics (1 = self)")
abline(lm(x$Euc_dist ~ x$Jac))
title(main = paste0(pt,": Is Euclidean Distance Correlated with Jaccard Similarity"))
test = cor.test(x$Jac,x$Euc_dist)

p = round(test$p.value,digits = 4)
r = round(test$estimate,digits = 2)
text(0.14,35,
     paste0("R is ",r,"and R^2 is ",round(r*r,digits = 2)," with p-value ",p), cex = 0.8)
dev.off()


