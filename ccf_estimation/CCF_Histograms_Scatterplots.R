# Read in libraries and set seed

rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR") #already installed on Aaron's macbookpro
# Load the GenVisR package
library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2);library(GGally)
set.seed(426)

#read in the sample table with the sample names and its associated tissue types
setwd("~/aarons_FAP_github_repository/")
source("recent_Annotation.R")
recent = recent_Annotation("no")
recentWGS = recent[recent$DNA_WGS!="",]
recentWGS$samples = str_remove_all(string = recentWGS[,1],pattern = "-")


patients = c("A001","A002");setwd("~/Bulk_A001_A002/");path = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf."

# patients = c("EP","JP");setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/"); path = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100."


for (pt in patients) {
        # file<-read.table(file = paste0(path,pt,".txt"), header = T, sep = "\t")
        file<-read.table(file = paste0(path,pt,".txt"), header = T, sep = "\t")
        #Select only statistically interesting columns. Rename the columns so that are compatible with sample information
        file_2<-select(file, chr, pos, id, ref, alt, geneName, geneLoc, functionalClass, somatic, rep,ends_with(match = "ccf"),-mergeCCF,-contains("bpccf"))
                # filter(!(rep==1 & sc==1))%>% #remove low grade mutations-repeats and poor mapping?
                # filter(!(rep==1 & geneLoc!="exonic")) #only allow mutations in repeat regions in exons.
        
        colnames(file_2)<-str_replace_all(colnames(file_2),pattern = "ccf",replacement = "")
        
        #select only variants per polyp which are known to be somatic
        file_varreads_melt<-gather(file_2,"SampleName", "ccf",11:20)%>%
                filter(grepl(pt,SampleName))%>%
                filter(str_detect(somatic, paste0(SampleName,"\\","[")))%>% #filters for only somatic mutations per sample
                filter(ccf>=0.1)
          if (pt=="EP") {
            file_varreads_melt$SampleName = str_replace(file_varreads_melt$SampleName,paste0(pt,"."),"F_")
          } else if (pt=="JP") {
            file_varreads_melt$SampleName = str_replace(file_varreads_melt$SampleName,pt,"G_")
          }
        
        #Rename samples to include Stage and Size (mm)
        samples = unique(file_varreads_melt$SampleName)
        for (sample in samples) {
          sample.index = file_varreads_melt$SampleName %in% sample
          stage = recentWGS$Stage..Polyp..Normal..AdCa.[recentWGS$samples==sample]
          size = recentWGS$Size..mm.[recentWGS$samples==sample]
          file_varreads_melt$SampleName[sample.index] = paste0(sample,"_",stage,"_",size)
        }

        ###### plot all mutations
        ccf_hist <- ggplot(file_varreads_melt,aes(x=ccf, fill = functionalClass)) +
                geom_histogram() +
                coord_cartesian(ylim = c(500,0),xlim = c(0,1.10)) +
                facet_wrap( ~ SampleName,nrow = 4,scales = "free_y")
        ccf_hist
        ggsave(filename = paste0(pt,"_ccf_distributions_func.pdf"),plot = ccf_hist)
        
        ccf_hist <- ggplot(file_varreads_melt,aes(x=ccf)) +
          geom_histogram() +
          coord_cartesian(ylim = c(500,0),xlim = c(0,1.10)) +
          facet_wrap( ~ SampleName,nrow = 4,scales = "free_y")
        ccf_hist
        ggsave(filename = paste0(pt,"_ccf_distributions.pdf"),plot = ccf_hist)
        
        ##### plot only exonic mutations
        file_varreads_melt_exonic<-filter(file_varreads_melt,geneLoc=="exonic")
        ccf_hist <- ggplot(file_varreads_melt_exonic,aes(x=ccf, fill = functionalClass)) +
          geom_histogram() +
          facet_wrap( ~ SampleName,nrow = 4,scales = "free_y") +
          coord_cartesian(ylim = c(20,0))
          
        ccf_hist
        ggsave(filename = paste0(pt,"_ccf_distributions_exonic.pdf"),plot = ccf_hist)
}

##scatter plots of CCFs with GGally
#https://drsimonj.svbtle.com/pretty-scatter-plots-with-ggplot2
for (pt in patients) {
  file<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
  file<-select(file,"chr","pos","id","ref","alt",ends_with("ccf"),-mergeCCF,functionalClass,-contains("bpccf"))
  
  #change the EP and JP values
  if (pt=="EP") {
    colnames(file) = str_replace(colnames(file),paste0(pt,"."),"F_")
  } else if (pt=="JP") {
    colnames(file) = str_replace(colnames(file),pt,"G_")
  }
  
  ## http://ggobi.github.io/ggally/#custom_functions

  #Rename samples to include Stage and Size (mm)
  samples = unique(file_varreads_melt$SampleName)
  for (sample in samples) {
    sample.index = file_varreads_melt$SampleName %in% sample
    stage = recentWGS$Stage..Polyp..Normal..AdCa.[recentWGS$samples==sample]
    size = recentWGS$Size..mm.[recentWGS$samples==sample]
    file_varreads_melt$SampleName[sample.index] = paste0(sample,"_",stage,"_",size)
  }
  
  
## i know this one works  
  my_bin <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
    ggplot(data = data, mapping = mapping) +
      geom_bin2d() +
      scale_fill_gradient(low = low, high = high) +
      theme(strip.placement = "outside")
  }
  
  ## All mutations pairwise comparison
  pm <- ggpairs(data = file,columns = c(7:ncol(file)-1),
                upper = "blank",
                diag = NULL,
                lower = list(continuous = my_bin),switch = "both")
  #  pm
  ggsave(filename = paste0(pt,"_pairwise_CCF.pdf"),plot = pm,width = 20,height = 20)
  
  #exonic filter
  file_exonic<-filter(file,!is.na(functionalClass))%>%
    select(-functionalClass)
  pm <- ggpairs(data = file_exonic,columns = c(6:ncol(file_exonic)),
                upper = "blank",
                diag = NULL,
                lower = list(continuous = my_bin),switch = "both")
#  pm
  ggsave(filename = paste0(pt,"_pairwise_CCF_exonic.pdf"),plot = pm,width = 20,height = 20)

  ###chr7 filter
  file_chr7<-filter(file,chr=="chr7")%>%
    select(-functionalClass)
  pm <- ggpairs(data = file_chr7,columns = c(6:ncol(file_chr7)),
                upper = "blank",
                diag = NULL,
                lower = list(continuous = my_bin),switch = "both")
  #  pm
  ggsave(filename = paste0(pt,"_pairwise_CCF_chr7.pdf"),plot = pm,width = 20,height = 20)
  
  
  }


#################### remake of CCF distribution plots

#read in the sample table with the sample names and its associated tissue types
setwd("~/aarons_FAP_github_repository/")
source("recent_Annotation.R")
recent = recent_Annotation("no")
recent = data.frame(recent,stringsAsFactors = F)

pts = c("A001","A002"); setwd("~/Bulk_A001_A002/");pt = "A001"

  file<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
  # file<-select(file,"chr","pos","id","ref","alt",ends_with("ccf"),-mergeCCF,functionalClass,-contains("bpccf"))
  
  #change the EP and JP values
  if (pt=="EP") {
    colnames(file) = str_replace(colnames(file),paste0(pt,"."),"F_")
  } else if (pt=="JP") {
    colnames(file) = str_replace(colnames(file),pt,"G_")
  }
  
  ## http://ggobi.github.io/ggally/#custom_functions
  
  #Rename samples to include Stage and Size (mm)
  samples = str_remove(string = colnames(file)[str_detect(string = colnames(file),pattern = "ccf")],pattern = "ccf")
  recent$SampleName = str_remove_all(string = recent$SampleName,pattern = "-")
  for (sample in samples) {
    sample.index = recent$SampleName %in% sample
    stage = recent$Stage..Polyp..Normal..AdCa.[recent$SampleName == sample]
    size = recent$Size..mm.[recent$SampleName == sample]
    recent$SampleName[recent$SampleName == sample] = paste0(sample,"_",stage,"_",size)
  }

polyp1 = "A001C004"
polyp2 = "A001C007"

polyp1ccf = file[,colnames(file) %in% paste0(polyp1,"ccf")]
polyp2ccf = file[,colnames(file) %in% paste0(polyp2,"ccf")]

plot(polyp2ccf,polyp1ccf,xlab = polyp2, ylab = polyp1)
  
  
  
# https://www.reed.edu/data-at-reed/resources/R/loops_with_ggplot2.html
# for making an iterative for loop for plotting
  
  