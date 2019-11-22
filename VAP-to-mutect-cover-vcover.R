rm(list=ls())
#For A001,A002
setwd("/Users/ahorning/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")

#For EP and JP WGS data
# setwd("/Users/ahorning/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/")
#########
# Consider something like this to collect all of the mutations per sample
# x = data.frame(matrix(c(1,2,3,10),byrow = T,nrow = 2),stringsAsFactors = F)
# y = data.frame(matrix(c("a","b","c","z"),byrow = T,nrow = 2),stringsAsFactors = F)
# y[x<=1 | x==10 | y=="b"] = "aaron"
#######

source("~/aarons_FAP_github_repository/vapReduce_cover,vcover,mutect.R")

##UchiCase2
#maf <- read.table("../../UchiCase2_MuTectSNV_Indel_Coding_Filtered.txt",header=T)
    #maf <- maf[which(maf$CtoT_filter == "Pass"),]

#For A001,A002
  # vap = read.table(file = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified",
  # header = T,sep = "\t")

  
# For EP and JP WGS
# vap = read.table(file = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt",
#                  header = T,sep = "\t")

######### for practice
#vap = head(vap,100)

pts = c("A001", "A002")
pts = "A001"
pt = "A002"
# pts = c("EP","JP")
for (pt in pts) {
  
  path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
  vap =  read.table(file = path,
                    header = T,sep = "\t")
  
  #creates maf table and sample name lists
  maf = vapReduce_cover_vcover_mutect(vap,pt,"ccf")
  message(paste("Patient",pt,"table chosen"))
  
  
  #number of mutations
  mnumber = nrow(maf$phylip.input.maf.file)
  message(paste(mnumber, "of mutations to sift through. phew thats a lot!"))
  
  #Creates an empty "Sample vector" for each sample name
  for (sample in maf$patient.samplenames) {
    sample = str_replace_all(string = sample,pattern = "-",replacement = "\\.") #for EP and JP VAP files
    assign(sample,as.vector(rep("0",mnumber)))
  }
  message(paste("made",length(maf$patient.samplenames)," sample vectors for patient",pt))
  
  #Because indels and SNVs are included, it is useful to use a single symbol in the FASTA files
  #That is why we are creating this makeshift set of basepair sequenced to represent our maf file
  var_allele = c()
  ref_allele = c()
  for (i in 1:mnumber){
    pair <- sample(c("T","C","G","A"),2)
    ref_allele = c(ref_allele,pair[1])
    var_allele = c(var_allele,pair[2])
  }
  message("made the makeshift ref and var alleles. so creative Zheng!")
  
  #Create little tables with the specific values of interest: 1 table of the depth, 1 table of VAF (maf), 1 table of the mutect calls ("yes","unknown")
  cover.columns = maf$phylip.input.maf.file%>%select(ends_with("d"))%>%select(-id)
  cover.columns[is.na(cover.columns)]<-0
  mutect.columns = maf$phylip.input.maf.file%>%select(ends_with("mutect"))
  freq.columns = maf$phylip.input.maf.file%>%select(contains("maf"))
  freq.columns[is.na(freq.columns)]<-0
  
  #Based on the individual columns values, determine whether or not there was a mutation, and add either "ref" or "alt" to the makeshift 
  # sample vector
  for (sample in 1:length(maf$patient.samplenames)) {
    message(paste("make sample specific columns. Up to bat is: ",maf$patient.samplenames[sample]))
    print(colnames(mutect.columns)[sample])
    mutect.col = mutect.columns[,sample]
    cover.col = cover.columns[,sample]
    freq.col = freq.columns[,sample]
    
    message(paste("Determine if the mutations for sample",maf$patient.samplenames[sample],"is real or nah"))
    for (k in 1:mnumber) {
      mutect = mutect.col[k]
      cover = cover.col[k]
      vcover = round(cover.col[k]*freq.col[k])
      
      if(mutect == "yes")
      {allele = toString(var_allele[k])}
      else if(mutect == "no") {
        if(cover>=10)
        {allele = toString(ref_allele[k])} #if total coverage of the reference allele is >10, then ref.
        else{
          if(vcover>=1) #else, if variant cover is greater than 1 (thus > 10% coverage)
          {allele = toString(var_allele[k])} #variant allele
          else
          {allele = toString("N")} #else, unknown.
        }
      }
      
      else{ #if mutect is neither "yes" nor "no", and allele call is unclear, then...
        if(cover>=10 & vcover>=3)
        {allele = toString(var_allele[k])}
        else if(cover>=10 & vcover==0)
        {allele = toString(ref_allele[k])}
        else if(cover<=10 & vcover>=1)
        {allele = toString(var_allele[k])}
        else if(cover>=20 & vcover==1)
        {allele = toString(ref_allele[k])}
        else
        {allele = toString("N")}
      }
      
      temp = get(paste0(maf$patient.samplenames[sample]))
      temp[k] = allele
      assign(paste0(maf$patient.samplenames[sample]),temp)
    }
    print(head(temp,30))
    message("on to the next sample")
  }
  
  ########################## Save the fasta files by patient
  message(paste("Creating and saving the fasta file for patient",pt))
  # sink(paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.phylipinput.",pt,".fasta"))
  sink(paste0("phylip/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt",".fasta"))
  if (pt=="EP" | pt=="JP") {
    cat(">Normal")
    cat("\n")
    cat(paste(ref_allele,collapse=""))
    cat("\n")
    for (sample in 1:length(maf$patient.samplenames)) {
      temp = get(paste0(maf$patient.samplenames[sample]))
      str_flatten(temp)
      
      cat(paste0(">",maf$patient.samplenames[sample])) #name of the sample
      cat("\n") #creates a next row in the file
      cat(str_flatten(temp))
      cat("\n")
    }
  }else{
    cat(">Normal")
    cat("\n")
    cat(paste(ref_allele,collapse=""))
    cat("\n")
    for (sample in 1:length(maf$patient.samplenames)) {
      temp = get(paste0(maf$patient.samplenames[sample]))
      str_flatten(temp)
      
      cat(paste0(">",maf$patient.samplenames[sample])) #name of the sample
      cat("\n") #creates a next row in the file
      cat(str_flatten(temp))
      cat("\n")
    }
  }  
  
  sink()
  message(paste("done saving the fasta file for patient",pt))
}

save.image("Pre_phylip_input.RData")


# May be useful in the future.... ####################
# cat(">Normal")
# cat("\n")
# cat(paste(ref_allele,collapse=""))
# cat("\n")

