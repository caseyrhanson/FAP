rm(list=ls())
library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2);library(GGally);library(vioplot)
set.seed(426)



################################Functions#########################################

full_join_remove_repeats = function(data,data.ccf,by) {
  library(dplyr);library(tidyr);library(stringr)
  #merge data tables by mut.ID
  x = full_join(x = data,y = data.ccf,by = by)
  #remove .x from duplicate column names
  colnames(x) = str_remove(string = colnames(x),pattern = ".x")
  #removes repeat columns which ended in ".y"
  joined.table = x[,which(str_detect(string = colnames(x),pattern = ".y",negate = T))]
  return(joined.table)
}

create_vector_of_sample_names = function(pt,data,blood.normal.sample.string,ccf.or.maf){
  # pt = "A001"
  # data = data
  # ccf.or.maf = "maf$"
  # blood.normal.sample.string = "blood"
  # Create samples vector of samples from CCF Table
  ccf.index = str_detect(string = colnames(data),pattern = ccf.or.maf)
  samples.ccf = colnames(data)[ccf.index]
  samples = str_remove(string = samples.ccf,pattern = ccf.or.maf)
  #remove bp label if there is one
  sample.bp = str_detect(string = samples,pattern = "bp",negate = T)
  samples = samples[sample.bp]
  #remove blood sample
  samples.blood = str_detect(string = samples,pattern = blood.normal.sample.string)
  samples = samples[!samples.blood]
  #collect samples names of specific patient
  samples.pt = str_detect(string = samples,pattern = pt)
  samples = samples[samples.pt]
  return(samples)
}

get_cancer_drivers = function(){
  x = read.csv(file = "~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/CancerDriverGenes/PanCanDrivers_Cell2018.csv",
             header = T,sep = ",",stringsAsFactors = F,skip = 3)
  listofgenes = x$Gene
  return(listofgenes)
}

# CCFs above below 60%
add.percent.above.below = function(percent,vector,side,color){
  #percent = 0.6
  #vector = truncal.ccfs
  #side = "right"
  #color = "blue"
  abline(h = percent,col = color,lty = 2)
  
  total.muts = length(vector)
  high.muts = sum(vector>=percent)
  low.muts = sum(vector<percent)
  
  high.muts.percent = round(high.muts/total.muts*100)
  low.muts.percent = round(low.muts/total.muts*100)
  
  side = ifelse(side=="right",1.4,
             ifelse(side=="left",0.6))
  
  text(x = side,y = .8,labels = paste0(high.muts," (",high.muts.percent,"%)"),cex = .8,col = color)
  text(x = side,y = .35,labels = paste0(low.muts," (",low.muts.percent,"%)"),cex = .8,col = color)
}
##################################################################################

#which patient do we want to analyze first
pt = "A001";setwd("~/Bulk_A001_A002/")

# # path = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
# path = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt")
# #read in MAF data table
# data = read.table(file = path,header = T,sep = "\t",stringsAsFactors = F)
# #make a new united column: mut.ID
# data = unite(data,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
# any(duplicated(data$mut.ID)) #FALSE only if all the mut.ID are distinct

#read in CCF data table
path.ccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
data.ccf = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
#make a new united column: mut.ID
data.ccf = unite(data.ccf,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
any(duplicated(data.ccf$mut.ID)) #FALSE only if all the mut.ID are distinct

#Read in treeomics data for specific patient; treeomics has the "truncal" mutations
#read in treeomics data specifically performed on "ccf adjusted"
variants.path = list.files(path = paste0("~/Bulk_A001_A002/treeomics/",pt,"_ccf_adjusted/"),pattern = paste0("_variants.csv"),full.names = T)
treeomics = read.table(file = variants.path, header = T,sep = ",")
head(colnames(treeomics))
treeomics = unite(treeomics,col = "mut.ID",c("Chromosome","StartPosition","RefAllele","AltAllele"),sep = ":")
any(duplicated(treeomics$mut.ID))

#Find mutations from Treeomics output which are "Trunk"
truncal.muts = treeomics$mut.ID[treeomics$Phylogeny=="Trunk"] #mutations found in treeomics to be "trunkal"
truncal.data.ccf = data.ccf[data.ccf$mut.ID %in% truncal.muts,] #find the above mutations in the original CCF-adjusted file
# truncal.data = data[data$mut.ID %in% truncal.muts,]

#create vector of sample names from Treeomics File
colnames(treeomics)
samples = create_vector_of_sample_names(pt = "A001",
                                        data = treeomics,
                                        blood.normal.sample.string = "blood",
                                        ccf.or.maf = "VAF_")
#plot the trunkal mutations CCFs
sample = 1;layout(1)
plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)
for (sample in 1:length(samples)) {
  plot(1,1,type = "n", xlim = c(0.5,1.5), ylim = c(0,1),
       main = paste0(samples[sample], " \n Trunk Mutation CCF"),
       ylab = "CCF",xlab = paste0(samples[sample]),xaxt = "n")
  truncal.ccfs = truncal.data.ccf[,paste0(samples[sample],"ccf")]
  myjitter = jitter(x = rep(1,length(truncal.ccfs)),amount = .25,factor = 1)
  points(x = myjitter,y = truncal.ccfs,pch = 21,bg = "green",col = "grey")
  
  #Add line and percent values above and below to the right of the chart
  add.percent.above.below(0.6,truncal.ccfs,side = "right",color = "blue")
  add.percent.above.below(0.5,truncal.ccfs,side = "left",color = "magenta")
}

# Plot the Trunkal Mutations MAFs*2
sample = 1
plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)
for (sample in 1:length(samples)) {
  plot(1,1,type = "n", xlim = c(0.5,1.5), ylim = c(0,1),
       main = paste0(samples[sample], " \n Trunk Mutation CCF (MAF*2)"),
       ylab = "MAF * 2",xlab = paste0(samples[sample]),
       xaxt = "n") #
  truncal.mafs = truncal.data[,paste0(samples[sample],"maf")]
  myjitter = jitter(x = rep(1,length(truncal.mafs)),amount = .25,factor = 1)
  #Multiply the MAFs by 2. If it is larger than 1, bring value down to 1
  truncal.mafs=truncal.mafs*2
  truncal.mafs[truncal.mafs>=1] = 1
  points(x = myjitter,y = truncal.mafs,pch = 21,bg = "green",col = "grey")
  abline(h = 0.6,col = "blue",lty = 2)
  
  total.muts = length(truncal.mafs)
  high.muts = sum(truncal.mafs>=0.6)
  low.muts = sum(truncal.mafs<0.6)
  
  high.muts.percent = round(high.muts/total.muts*100)
  low.muts.percent = round(low.muts/total.muts*100)
  
  text(x = 1.4,y = .8,labels = paste0(high.muts," (",high.muts.percent,"%)"),cex = .8)
  text(x = 1.4,y = .35,labels = paste0(low.muts," (",low.muts.percent,"%)"),cex = .8)
}


################  All Phylogenetic-types of mutations
################ Trunkal, private, shared ###############

#get vector of sample names
samples = create_vector_of_sample_names(pt = "A001",
                                        data = treeomics,
                                        blood.normal.sample.string = "blood",
                                        ccf.or.maf = "VAF_")

#Find mutations from Treeomics output which are "Trunk"
truncal.muts = treeomics$mut.ID[treeomics$Phylogeny=="Trunk"]
private.muts = treeomics$mut.ID[treeomics$Phylogeny=="Private"]
shared.muts = treeomics$mut.ID[treeomics$Phylogeny=="Shared"]

#add phylogenetic types of muts to data
data$Phylogeny = rep("Not Determined",nrow(data))
data$Phylogeny[data$mut.ID %in% truncal.muts] = "Trunk"
data$Phylogeny[data$mut.ID %in% private.muts] = "Private"
data$Phylogeny[data$mut.ID %in% shared.muts] = "Shared"
str(data$Phylogeny)

#change the phylogeny mutation types
data$Phylogeny = factor(data$Phylogeny,levels = c("Trunk","Shared","Private","Not Determined"))

table(data$Phylogeny)


############# All Mutations in Treeomics ############
svg(filename = paste0("treeomics/",pt,"_treeomics_AllMutations_Trunk_Shared_Private.svg"),width = 15,height = 10)
plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)

for (i in 1:length(samples)) {
  #adjust maf: multiply by 2 and anything larger than 1 is equal to 1
  sample.maf.index = str_detect(string = colnames(data), pattern = paste0(samples[i],"maf"))
  maf_x_2 = data[,sample.maf.index] *2
  maf_x_2[maf_x_2 >= 1] = 1
  
  boxplot(maf_x_2 ~ data$Phylogeny, ylab = "CCF (MAF *2)",
          main = paste0(samples[i])) #make basic box plot
  abline(h = 0.6,col = "blue",lty = 2) # add line at 60% CCF
  
  # Statistics on Trunk to figure
  type = c("Trunk","Shared", "Private", "Not Determined")
  type.plot.loc = c(1.5, 2.5, 3.5, 4.5)
  names(type.plot.loc) = type
  
  for (mut.type in type) {
    total.type.muts = length(maf_x_2[data$Phylogeny == mut.type])
    high.trunk.muts = sum(maf_x_2[data$Phylogeny == mut.type] >=0.6)
    low.trunk.muts = sum(maf_x_2[data$Phylogeny == mut.type] < 0.6)
    
    high.muts.percent = round(high.trunk.muts/total.type.muts*100)
    low.muts.percent = round(low.trunk.muts/total.type.muts*100)
    
    text(x = type.plot.loc[mut.type],
         y = .8,labels = paste0(high.trunk.muts," (",high.muts.percent,"%)"),
         cex = .8,font = 2,xpd = "NA") #adds high values
    text(x = type.plot.loc[mut.type],
         y = .4,labels = paste0(low.trunk.muts," (",low.muts.percent,"%)"),
         cex = .8,font = 2, xpd = "NA") #adds low values
  }
}

dev.off()


############# Only Exonic Mutations in Treeomics ############
data.exon = data[data$geneLoc == "exonic",]

svg(filename = paste0("treeomics/",pt,"_treeomics_ExonicMutations_Trunk_Shared_Private.svg"),width = 15,height = 10)
plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)

for (i in 1:length(samples)) {
  #adjust maf: multiply by 2 and anything larger than 1 is equal to 1
  sample.maf.index = str_detect(string = colnames(data.exon), pattern = paste0(samples[i],"maf"))
  maf_x_2 = data.exon[,sample.maf.index] *2
  maf_x_2[maf_x_2 >= 1] = 1
  
  boxplot(maf_x_2 ~ data.exon$Phylogeny, ylab = "CCF (MAF *2)",
          main = paste0(samples[i])) #make basic box plot
  abline(h = 0.6,col = "blue",lty = 2) # add line at 60% CCF
  
  # Statistics on Trunk to figure
  type = c("Trunk","Shared", "Private", "Not Determined")
  type.plot.loc = c(1.5, 2.5, 3.5, 4.5)
  names(type.plot.loc) = type
  
  for (mut.type in type) {
    total.type.muts = length(maf_x_2[data.exon$Phylogeny == mut.type])
    high.trunk.muts = sum(maf_x_2[data.exon$Phylogeny == mut.type] >=0.6)
    low.trunk.muts = sum(maf_x_2[data.exon$Phylogeny == mut.type] < 0.6)
    
    high.muts.percent = round(high.trunk.muts/total.type.muts*100)
    low.muts.percent = round(low.trunk.muts/total.type.muts*100)
    
    text(x = type.plot.loc[mut.type],
         y = .8,labels = paste0(high.trunk.muts," (",high.muts.percent,"%)"),
         cex = .8,font = 2,xpd = "NA") #adds high values
    text(x = type.plot.loc[mut.type],
         y = .4,labels = paste0(low.trunk.muts," (",low.muts.percent,"%)"),
         cex = .8,font = 2, xpd = "NA") #adds low values
  }
}

dev.off()



mafs = data[,str_detect(colnames(data),"maf")]
head(mafs)
mafs = mafs*2
head(mafs)
mafs[mafs>=1] = 1
datawoMAFs = data[,str_detect(colnames(data),"maf",negate = T)]

data2 = cbind(datawoMAFs,mafs)
drivers = get_cancer_drivers() #get list of drivers
data2$cancerdriver = rep("Not Driver",nrow(data2)) # create new column in data2
data2$cancerdriver[data2$geneName %in% drivers] = "Driver" #add string "Driver" to genes which are known drivers

data2.exonic = data2[data2$geneLoc=="exonic",]

i=1
for (i in 1:length(samples)) {
  plot = ggplot(data = data2.exonic,mapping = aes(x = Phylogeny))+
    geom_boxplot(aes_string(y = paste0(samples[i],"maf")),outlier.shape = NA)+
    ggtitle(paste0(samples[i]," Exonic Mutations \n and Labelled Cancer Driver Genes")) +
    theme(plot.title = element_text(size = 30, face = "bold"))+
    labs(y=paste0(samples[i],"CCF (MAF * 2)"), x = "Phylogenetic Mutation Type")+
    geom_jitter(aes_string(y = paste0(samples[i],"maf"),colour = "cancerdriver"))+
    ggrepel::geom_label_repel(data = data2.exonic[which(data2.exonic[,paste0(samples[i],"maf")] > 0 & data2.exonic$cancerdriver=="Driver"),],
                              aes_string(y = paste0(samples[i],"maf"),label = "geneName"),size = 5)
  
  plot
  ggsave(filename = paste0("treeomics/",samples[i],"_treeomics_ExonicMutations_Trunk_Shared_Private.pdf"))
}


