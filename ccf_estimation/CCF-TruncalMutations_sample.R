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
  # pt = "EP"
  # data = treeomics
  # ccf.or.maf = "VAF_"
  # blood.normal.sample.string = "NL"
  # Create samples vector of samples from CCF Table
  ccf.index = str_detect(string = colnames(data),pattern = ccf.or.maf)
  samples.ccf = colnames(data)[ccf.index]
  samples = str_remove(string = samples.ccf,pattern = ccf.or.maf)
  #remove bp label if there is one
  sample.bp = str_detect(string = samples,pattern = "bp",negate = T)
  samples = samples[sample.bp]
  #remove blood or normal sample
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
  
  #becaus numeric vector now has "nas" in it, remove them and get the counts
  nonzeroccfs = vector[!is.na(vector)]
  total.muts = length(nonzeroccfs)
  high.muts = sum(vector>=percent,na.rm = T)
  low.muts = sum(vector<percent,na.rm = T)
  
  high.muts.percent = round(high.muts/total.muts*100)
  low.muts.percent = round(low.muts/total.muts*100)
  
  side = ifelse(side=="right",1.4,
             ifelse(side=="left",0.6))
  
  text(x = side,y = .8,labels = paste0(high.muts," (",high.muts.percent,"%)"),cex = 1,col = color)
  text(x = side,y = .35,labels = paste0(low.muts," (",low.muts.percent,"%)"),cex = 1,col = color)
}
##################################################################################

pt = "A001"
pt = "JP"
pts = c("A001","A002","EP","JP")
pts = c("EP","JP")
if (pt %in% c("A001","A002")) {
  #which patient do we want to analyze first
  setwd("~/Bulk_A001_A002/")
  
  #read in CCF data table
  path.ccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
  data.ccf = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
  #make a new united column: mut.ID
  data.ccf = unite(data.ccf,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
  any(duplicated(data.ccf$mut.ID)) #FALSE only if all the mut.ID are distinct
  
  #Read in treeomics data for specific patient; treeomics has the "truncal" mutations
  #read in treeomics data specifically performed on "ccf adjusted"
  variants.path = list.files(path = paste0("~/Bulk_A001_A002/treeomics/",pt,"_ccf_adjusted"),pattern = paste0("_variants.csv"),full.names = T)
  treeomics = read.table(file = variants.path, header = T,sep = ",")
  head(colnames(treeomics))
  treeomics = unite(treeomics,col = "mut.ID",c("Chromosome","StartPosition","RefAllele","AltAllele"),sep = ":")
  any(duplicated(treeomics$mut.ID))
} else if (pt %in% c("EP","JP")) {
  #which patient do we want to analyze first
  setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")
  
  #read in CCF data table
  path.ccf = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt")
  data.ccf = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
  #make a new united column: mut.ID
  data.ccf = unite(data.ccf,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
  any(duplicated(data.ccf$mut.ID)) #FALSE only if all the mut.ID are distinct
  
  #Read in treeomics data for specific patient; treeomics has the "truncal" mutations
  #read in treeomics data specifically performed on "ccf adjusted"
  variants.path = list.files(path = paste0("treeomics_WGS_output/",pt,"_WGS_ccf"),pattern = paste0("_variants.csv"),full.names = T)
  treeomics = read.table(file = variants.path, header = T,sep = ",")
  head(colnames(treeomics))
  treeomics = unite(treeomics,col = "mut.ID",c("Chromosome","StartPosition","RefAllele","AltAllele"),sep = ":")
  any(duplicated(treeomics$mut.ID))
  if(!grepl(pattern = "chr",x = treeomics$mut.ID[1])){
    treeomics$mut.ID = paste0("chr",treeomics$mut.ID)
  } else {}
}

# mutations.not.in.treeomics = data.ccf[!(data.ccf$mut.ID %in% treeomics$mut.ID),]

#Find mutations from Treeomics output which are "Trunk"
truncal.muts = treeomics$mut.ID[treeomics$Phylogeny=="Trunk"] #mutations found in treeomics to be "trunkal"
truncal.data.ccf = data.ccf[data.ccf$mut.ID %in% truncal.muts,] #find the above mutations in the original CCF-adjusted file
# truncal.data = data[data$mut.ID %in% truncal.muts,]

#create vector of sample names from Treeomics File
colnames(treeomics)
if(pt %in% c("A001","A002")){
  samples = create_vector_of_sample_names(pt = pt,
                                          data = treeomics,
                                          blood.normal.sample.string = "blood",
                                          ccf.or.maf = "VAF_")
} else if(pt %in% c("EP", "JP")){
  samples = create_vector_of_sample_names(pt = pt,
                                          data = data.ccf,
                                          blood.normal.sample.string = "NL",
                                          ccf.or.maf = "ccf$")
  #The EP samples have periods in weird spaces and i removed them for the Treeomics files
  # samples = str_remove_all(string = samples,pattern = "\\.")
}




################  All Phylogenetic-types of mutations
################ Trunkal, private, shared ###############

# #get vector of sample names
# samples = create_vector_of_sample_names(pt = "A001",
#                                         data = treeomics,
#                                         blood.normal.sample.string = "blood",
#                                         ccf.or.maf = "VAF_")

samples

#Find mutations from Treeomics output which are "Trunk"
truncal.muts = treeomics$mut.ID[treeomics$Phylogeny=="Trunk"]
private.muts = treeomics$mut.ID[treeomics$Phylogeny=="Private"]
shared.muts = treeomics$mut.ID[treeomics$Phylogeny=="Shared"]

#add phylogenetic types of muts to data
data.ccf$Phylogeny = rep("Not Determined",nrow(data.ccf))
data.ccf$Phylogeny[data.ccf$mut.ID %in% truncal.muts] = "Trunk"
data.ccf$Phylogeny[data.ccf$mut.ID %in% private.muts] = "Private"
data.ccf$Phylogeny[data.ccf$mut.ID %in% shared.muts] = "Shared"
str(data.ccf$Phylogeny)
table(data.ccf$Phylogeny)

#Adjust the levels of the phylogeny mutation types
data.ccf$Phylogeny = factor(data.ccf$Phylogeny,levels = c("Trunk","Shared","Private","Not Determined"))
table(data.ccf$Phylogeny)





############ ggplots for All Mutations_ Trunk, shared, Private
data.ccf.original = data.ccf #for practice purposes
# data.ccf = data.ccf.original

drivers = get_cancer_drivers() #get list of drivers
data.ccf$cancerdriver = rep("Not Driver",nrow(data.ccf)) # create new column in data2

# #add string "Driver" to genes which are known drivers
# data.ccf$cancerdriver[data.ccf$geneName %in% drivers] = "Driver" 

#relabel what is Driver: Only Known Drivers with >= 20 CADD score
data.ccf$CADD_phred[is.na(data.ccf$CADD_phred)] = 0 #change NAs to "0"
data.ccf$cancerdriver[data.ccf$geneName %in% drivers & data.ccf$CADD_phred >= 20] = "Driver"

# samples = create_vector_of_sample_names(pt = "A001",
#                                         data = treeomics,
#                                         blood.normal.sample.string = "blood",
#                                         ccf.or.maf = "VAF_")

samples

#For every sample CCF column, replace every "0" with an "NA"
sample.col.index = colnames(data.ccf) %in% paste0(samples,"ccf")
bool.0.ccf = data.ccf[,sample.col.index]==0
data.ccf[sample.col.index][bool.0.ccf] = NA

#Is the mutation present: does it have enough "altc" and "d"
sample.d.index = colnames(data.ccf) %in% paste0(samples,"d") #get all depths
d = data.ccf[,sample.d.index]
sample.altc.index = colnames(data.ccf) %in% paste0(samples,"altc") #get all altcs
altc = data.ccf[,sample.altc.index]
dim(d)==dim(altc) #check the dimensions are the same
mut.present = matrix(nrow = dim(d)[1],ncol = dim(d)[2]) #make new matrix
colnames(mut.present) = paste0(samples,"PresentOrNot")

#use the 2 logicals from Zheng to determine mutations present-ness
mut.present[d>=10 & altc>=3 | d<=10 & altc >=1] = "Present"
mut.present[is.na(mut.present)] = "Not Present"
mut.present = as.data.frame(mut.present)


#Bring the mut.present columns together withthe data.ccf table
data.ccf = cbind(data.ccf,mut.present)
layout(1)
layout.show(1)
i=1
mutlocation = "all"
for(mutlocation in c("all","exonic")){
  ##set up mutlocation
  if (mutlocation=="exonic") {
    data.ccf = data.ccf[data.ccf$geneLoc==mutlocation,]
  }else{
    data.ccf=data.ccf
  }
  #Make plots for each set of mutlocations
  for (i in 1:length(samples)) {
    plot = ggplot(data = data.ccf,mapping = aes(x = Phylogeny))+
      geom_violin(aes_string(y = paste0(samples[i],"ccf")),outlier.shape = NA)+
      ggtitle(paste0(samples[i]," ",mutlocation," Mutations \n and Labelled Cancer Driver Genes CADD >= 20")) +
      theme(plot.title = element_text(size = 20, face = "bold"))+
      labs(y=paste0(samples[i],"CCF"), x = "Phylogenetic Mutation Type")+
      geom_jitter(aes_string(y = paste0(samples[i],"ccf"),colour = "cancerdriver",
                             size = paste0(samples[i],"PresentOrNot")))+
      scale_size_discrete(range = c(0.2,1))+
      ggrepel::geom_label_repel(data = data.ccf[which(data.ccf[,paste0(samples[i],"ccf")] > 0 & data.ccf$cancerdriver=="Driver"),],
                                aes_string(y = paste0(samples[i],"ccf"),label = "geneName"),size = 7)
    
    plot
    
    if(pt %in% c("A001", "A002")){
      ggsave(filename = paste0("treeomics/",samples[i],"_treeomics_",mutlocation,"CADD>= 20 Mutations_Trunk_Shared_Private.pdf"),
             width = 12,height = 10)
    } else if(pt %in% c("EP","JP")){
      ggsave(filename = paste0("treeomics_WGS_output/",pt,"_WGS_ccf/",samples[i],"_treeomics_",mutlocation,"CADD>=20_Mutations_Trunk_Shared_Private.pdf"),
             width = 12,height = 10)
    }
  }
}


###########
# Using base R
#plot the trunkal mutations CCFs
if(pt %in% c("A001","A002")){
  svg(filename = paste0("treeomics/",pt,"ccf_adjusted_treeomics_Trunk.svg"),width = 15,height = 10)
} else if(pt %in% c("EP", "JP")){
  svg(filename = paste0("treeomics_WGS_output/",pt,"_WGS_ccf/",pt,"_ccf_adjusted_treeomics_Trunk.svg"),width = 15,height = 10)
}
dev.off()
sample = 1;layout(1)
plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)
for (sample in 1:length(samples)) {
  plot(1,1,type = "n", xlim = c(0.5,1.5), ylim = c(0,1),
       main = paste0(samples[sample], " \n Trunk Mutation CCF"),
       ylab = "CCF",xlab = paste0(samples[sample]),xaxt = "n")
  #CCFs
  truncal.ccfs = truncal.data.ccf[,paste0(samples[sample],"ccf")]
  is.na(truncal.ccfs) = truncal.ccfs==0
  myjitter = jitter(x = rep(1,length(truncal.ccfs)),amount = .25,factor = 1)
  
  #Make 2 samples specific vectors for Depths and Altc
  color = rep("color",nrow(truncal.data.ccf))
  altc = truncal.data.ccf[,paste0(samples[sample],"altc")]
  d = truncal.data.ccf[,paste0(samples[sample],"d")]
  
  #determine Mutation depth validity (Got 2 logicals from Zheng for Phylip prep)
  color[d>=10 & altc>=3 | d<=10 & altc >=1] = "green" #Likely mutation
  color[color=="color"] = "grey" #Unlikely Mutation
  
  #Add points to graph. Color is determined above based on depth and altc
  points(x = myjitter,y = truncal.ccfs,pch = 21,bg = color,col = "darkgrey")
  
  #Add line and percent values above and below to the right of the chart
  add.percent.above.below(0.6,truncal.ccfs,side = "right",color = "blue")
  add.percent.above.below(0.5,truncal.ccfs,side = "left",color = "magenta")
}
dev.off()

############# All Mutations in Treeomics ############
#plot the trunkal mutations CCFs
if(pt %in% c("A001","A002")){
  svg(filename = paste0("treeomics/",pt,"_treeomics_AllMutations_Trunk_Shared_Private.svg"),width = 15,height = 10)
} else if(pt %in% c("EP", "JP")){
  svg(filename = paste0("treeomics_WGS_output/",pt,"_WGS_ccf/",pt,"_treeomics_AllMutations_Trunk_Shared_Private.svg"),width = 15,height = 10)
}

plot.matrix = matrix(data = 1:12,nrow = 3,ncol = 4,byrow = T)
layout(mat = plot.matrix)
layout.show(12)
i=1
for (i in 1:length(samples)) {
  #Choose sample and associated ccf values
  sample.ccf.index = str_detect(string = colnames(data.ccf), pattern = paste0(samples[i],"ccf$"))
  
  #Make all "0" values into NA so they are not plotted.
  ccf = data.ccf[,sample.ccf.index]
  is.na(ccf) = ccf==0 #make all "0" value ccfs into NA
  
  #make boxplot
  boxplot(ccf ~ data.ccf$Phylogeny, ylab = "CCF",
          main = paste0(samples[i])) #make basic box plot
  abline(h = 0.5,col = "blue",lty = 2) # add line at 60% CCF
  
  
  # Where to place stats on boxplot
  type = c("Trunk","Shared", "Private", "Not Determined")
  type.plot.loc = c(1.2, 2.2, 3.2, 4.2)
  names(type.plot.loc) = type
  
  #side
  side="left"
  mut.type=type[2]
  for (mut.type in type) {
    
    mut.ccf = ccf[data.ccf$Phylogeny==mut.type]
    
    is.na(mut.ccf) = mut.ccf==0 #make all "0" value ccfs into NA
    nonzeroccfs = mut.ccf[!is.na(mut.ccf)]
    total.muts = length(nonzeroccfs) #number of nonzero ccfs
    high.muts = sum(mut.ccf>=0.5,na.rm = T) #number of high ccfs
    low.muts = sum(mut.ccf<0.5,na.rm = T)
    
    high.muts.percent = round(high.muts/total.muts*100)
    low.muts.percent = round(low.muts/total.muts*100)
    
    side = ifelse(side=="right",1.4,
                  ifelse(side=="left",0.6,NA))
    
    text(x = type.plot.loc[mut.type],
         y = .8,labels = paste0(high.muts,"\n (",high.muts.percent,"%)"),
         cex = .8,font = 2,xpd = "NA",srt=45) #adds high values
    text(x = type.plot.loc[mut.type],
         y = .4,labels = paste0(low.muts,"\n (",low.muts.percent,"%)"),
         cex = .8,font = 2, xpd = "NA",srt=45) #adds low values
  }
}

dev.off()



# ######### KRAS mutation in A001
# #which patient do we want to analyze first
# pt = "A001";setwd("~/Bulk_A001_A002/")
# #read in CCF data table
# path.ccf = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
# data.ccf = read.table(file = path.ccf,header = T,sep = "\t",stringsAsFactors = F)
# #make a new united column: mut.ID
# data.ccf = unite(data.ccf,col = "mut.ID",c("chr","pos","ref","alt"),sep = ":")
# any(duplicated(data.ccf$mut.ID)) #FALSE only if all the mut.ID are distinct
# 
# kras = data.ccf[data.ccf$geneName=="KRAS",]
# 
# x = select(kras,1:2,contains("A001C007"),contains("A001C107"))
# 

