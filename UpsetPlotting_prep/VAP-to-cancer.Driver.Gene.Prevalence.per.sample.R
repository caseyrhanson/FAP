#This code creates mutation tables for each patient.
# Final object created is: tabl.of.muts
# JP3                                 JP31  JP38                                     JP34B
# 1 <NA> KRAS:NM_004985:exon2:c.G35A:p.G12D <NA> APC:NM_001127511:exon14:c.G3876T:p.K1292N
# 2 <NA>                               <NA> <NA>                                      <NA>
# 3 <NA>                               <NA> <NA>                                      <NA>

# makes this sort of table for every patient.


rm(list=ls())

library(stringr);library(dplyr);library(tidyr)

source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")
source("~/DNAseq_WGS/scripts/find.sample.and.tumor.names.R")

load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
ccf.data = list("A001" = ptA001, "A002"= ptA002, "EP" = ptF, "JP" = ptG)


# #read in VAP output file
# vap.path = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt"
# # vap.path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table" #EP and JP original WGS
# vap = read.table(vap.path,header = T,sep = "\t",stringsAsFactors = F)
pt = "EP"
pts = c("A001", "A002", "EP", "JP")
for (pt in pts) {
  vap = ccf.data[[pt]]
  
  
  #read in Comprehensive Cancer Drive list
  pan.cancer.driver.path = path.expand("~/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv")
  drivers= read.csv(pan.cancer.driver.path,header = T,skip = 3,stringsAsFactors = F)
  
  #create list of samples from VAP file
  patient = "*" #collects both A001 and A002
  names=find.sample.and.tumor.names(vap,patient,maforccf = "ccf")
  sample.names = names$sample.names
  normal.names=names$normal.names
  tumor.names=names$tumor.names
  
  #Filter only the rows in the VAP output which have geneNames in the PanCancer Driver gene list AND are non-synonymous or stop-gain Muts
  vap.drivers.only = vap[which(vap$geneName %in% drivers$Gene),]
  #Filter NA and synonymous mutations from the dataset.
  #Dont Filter for top 1% (CADD >= 20) of deleterious mutations
  vap.drivers.only = filter(vap.drivers.only,functionalClass!="synonymous SNV")#%>%filter(CADD_phred>=20) 
  
  #create a mutation ID column for the specific mutations
  mutID = unite(vap.drivers.only,col = mutID,c("chr", "pos", "ref", "alt","geneName","CADD_phred"),sep = ":")[,1]
  vap.drivers.only = cbind(vap.drivers.only,mutID) #add a column to the VAP table which is for the distinct mutations
  
  number.of.driver.genes.found.in.vap = length(sort(unique(vap.drivers.only$geneName)))
  print(paste0("There are ",number.of.driver.genes.found.in.vap," cancer driver genes with non-synonymous or stopgain mutations in this VAP output file"))
  number.of.driver.genes.mutations.found.in.vap = length(sort(unique(vap.drivers.only$mutID)))
  print(paste0("There are ",number.of.driver.genes.mutations.found.in.vap," distinct non-synonymous or stopgain mutations in this VAP output file"))
  
  #Creates an empty "Sample vector" for each sample name
  for (sample in sample.names) {
    sample = str_replace_all(string = sample,pattern = "-",replacement = "\\.") #for EP and JP VAP files
    assign(sample,character(0))
  }
  
  # Change class of column to Character
  vap.drivers.only$somatic = as.character(vap.drivers.only$somatic)
  
  #Create a vector of cancer driver somatic ("good","goodSub" only) mutations per sample name
  for (i in 1:nrow(vap.drivers.only)) {
    vap.drivers.only$somatic[i] = str_replace_all(string = vap.drivers.only$somatic[i],pattern = '\\[good]' ,replacement = "")
    vap.drivers.only$somatic[i] = str_replace_all(string = vap.drivers.only$somatic[i],pattern = "\\[goodSub]" ,replacement = "")
    
    #collect the samples with mutations
    samples.with.muts = str_split(vap.drivers.only$somatic[i],pattern = ",")[[1]] #split up the sample names into a vector
    samples.with.muts=samples.with.muts[1:length(samples.with.muts)-1] #the above command leaves an empty "" so this command removes it
    samples.with.muts=samples.with.muts[!grepl(pattern = "doubt",x = samples.with.muts)] #some of the somatic mutations are "doubt", so this removes those
    
    #collect the associated amino acid change
    aminoacidchange.multi = as.character(vap.drivers.only$AAChange[i])
    aminoacidchange = str_split(string = aminoacidchange.multi,pattern = ",")[[1]][1]
    
    #creates character vectors for each sample. Fills with somatic mutations
    for (sample in samples.with.muts) {
      sample = str_replace_all(string = sample,pattern = "-",replacement = "\\.") #for EP and JP VAP files
      temp = get(sample)
      # temp = c(temp,vap.drivers.only$geneName[i]) #for counting up genes with mutations (may have repeats if there are more than 1 mutation per gene)
      #temp = c(temp,as.character(vap.drivers.only$AAChange[i])) #for collecting distinct mutations ID. Change to AAChange ##################
      temp = c(temp,aminoacidchange) #for collecting distinct mutations ID. Change to AAChange ##################
      temp = sort(temp)
      assign(x = sample,value = temp)
    }
  }
  
  #Collect number of driver mutations in each sample
  number.of.mutations.per.sample = character(0)
  for (sample in sample.names) {
    sample = str_replace(string = sample,pattern = "-",replacement = "\\.") #for EP and JP VAP files
    number.of.mutations.per.sample = c(number.of.mutations.per.sample,length(get(sample)))
  }
  
  ##### Create a table of all the Somatic Cancer Driver genes in each sample
  
  #highest number of driver mutations per sample. Need this to make all of the columns the same length
  max.number.of.mutations = max(as.numeric(number.of.mutations.per.sample))
  
  #To get the table started
  temp = c(get(sample.names[1]),rep(NA,max.number.of.mutations-length(get(sample.names[1])))) #adjust the length of the column by adding "NA"s
  tabl.of.muts = data.frame(temp,fix.empty.names = T) #create beginning of mutation table
  
  #bind all of the other samples together into a table
  for (sample.num in 2:length(sample.names)){
    temp = c(get(sample.names[sample.num]),rep(NA,max.number.of.mutations-length(get(sample.names[sample.num]))))
    tabl.of.muts = cbind(tabl.of.muts,temp)
  }
  
  plot(sort(table(A001C007),decreasing = T))
  unique(sort(A001C007))
  
  colnames(tabl.of.muts) = sample.names #change the column names to sample names
  tabl.of.muts
  apply(tabl.of.muts, 2, FUN = function(x) sum(grepl("APC",x)))
  apply(tabl.of.muts, 2, FUN = function(x) sum(grepl("PTPRD",x)))
  
  # #For A001 and A002
  if (pt %in% c("A001","A002")) {
    write.table(x = tabl.of.muts,paste0("~/Bulk_A001_A002/repeatMutations/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.SomaticDrivers_SNVmutations",pt,".csv"),
                append = F,sep = ",",row.names = F)
  } else if (pt %in% c("EP","JP")) {
    write.table(x = tabl.of.muts,paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/repeatMutations/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.EPandJP.SomaticDrivers_SNVmutations",pt,".csv"),
                append = F,sep = ",",row.names = F)
  }
}
