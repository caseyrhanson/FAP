#https://bioconductor.org/packages/3.7/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf

rm(list=ls())

library(dplyr);library(tidyr);library(stringr);library(GenomicRanges)


setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")
# Pyclone Input #######

#Collect Amino Acid Change data for each of the Mutations in the data set
vap = read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table",header = T,sep = "\t",stringsAsFactors = F)
mutID.AA = unite(data = vap, col = "mutation_id",chr,pos,geneName,ref,alt,sep = ":")%>%
  select(mutation_id,AAChange)%>%
  separate(AAChange,"AAchange_1",",",extra = "drop")%>%
  unite(mutation_id_2,AAchange_1,mutation_id,remove = F,sep = ",")%>%
  mutate(mutation_id_2 = str_replace(mutation_id_2,"^,",replacement = ""))%>%
  select(mutation_id,"mutation_id_AAchange" = mutation_id_2,AAchange_1)
head(mutID.AA)

######## Read in the SNV data and reshape it so only ref and alt counts are included
pts=c("A001","A002")
pt="A001"
for (pt in pts) {
  # sampAB<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt"), header = T, sep = "\t")
  sampAB<-read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
  
  #select ref and alt columns and melt the table. Modify sample names too.
  #some files have this weird column called bprefc and bpaltc.
  if (sum(grepl("bprefc",colnames(sampAB)))==0) {
    snvs<-select(sampAB, chr, pos, id, ref, alt,geneName, ends_with("refc"), ends_with("altc"))%>%
      gather(key = sampleName,value = altc,ends_with("altc"))%>%
      gather(key = SampleName_2, value = refc, ends_with("refc"))%>%
      select(-SampleName_2)%>%
      filter(!(altc == 0))%>%
      separate(col = sampleName,into = "sampleName",sep = "altc")%>%
      mutate(sampleName_2=str_replace(string = sampleName,pattern = "\\.",replacement = "-"))%>%
      select(-sampleName)
  } else {
    snvs<-select(sampAB, chr, pos, id, ref, alt,geneName, ends_with("refc"), ends_with("altc"),-bprefc,-bpaltc)%>%
      gather(key = sampleName,value = altc,ends_with("altc"))%>%
      gather(key = SampleName_2, value = refc, ends_with("refc"))%>%
      select(-SampleName_2)%>%
      filter(!(altc == 0))%>%
      separate(col = sampleName,into = "sampleName",sep = "altc")%>%
      mutate(sampleName_2=str_replace(string = sampleName,pattern = "\\.",replacement = "-"))%>%
      select(-sampleName)
  }
  
  colnames(snvs)[colnames(snvs)=="sampleName_2"] = "sampleName"
  
  samples<-sort(unique(snvs$sampleName)) #list of sample names
  
  #i=1 
  for (i in 1:length(samples)) {
    snvsA = snvs[snvs$sampleName==samples[i],]
    
    
    ######### Read in CNV data
    cnvFiles<-list()
    for (sample in list.files("titan/segments/")) {
      if (!grepl(pattern = "_titan",x = sample) && grepl(pt,sample)) {
        print(sample)
        cnvFiles<-append(cnvFiles,sample)
      }
    }
    cnvFiles = sort(unlist(cnvFiles))
    #j=grep(pattern = samples[i],x = cnvFiles) #the sample being collected for SNVs is collected for CNVs too
    cnvA<-read.table(file = paste0("titan/segments/",cnvFiles[i]), header = T, sep = "\t")
    cnvA = mutate(cnvA,"sampleName_cnvfile"=cnvFiles[i])
    
    # add "chr" before every chromosome number
    cnvSeqNames = cnvA$chrom
    if (!grepl("chr",cnvA$chrom[1])){
      cnvSeqNames = paste("chr",cnvA$chrom,sep="")
    }
    snvSeqNames = snvsA$chr
    if (!grepl("chr",snvsA$chr[1])){
      snvSeqNames = paste("chr",snvsA$chr,sep="")
    }
    
    cnvRangeA_M = GRanges(seqnames = cnvSeqNames,
                          ranges = IRanges(cnvA$loc.start, end=cnvA$loc.end),
                          strand=rep('+',dim(cnvA)[1]),
                          normal_cn = cnvA$copynumber,
                          minor_cn = cnvA$minor_cn,
                          major_cn = cnvA$major_cn,
                          sample_cnv = cnvA$sampleName_cnvfile)
    
    snvRange_M = GRanges(seqnames = snvSeqNames,
                         ranges = IRanges(snvsA$pos, end=snvsA$pos),
                         strand=rep('+',dim(snvsA)[1]),
                         mcols = snvsA[,c(3:length(colnames(snvsA)))])

    x = mergeByOverlaps(query = snvRange_M,subject = cnvRangeA_M)
    
    
    y = as.data.frame(x)%>%
      select(chr = snvRange_M.seqnames,
             pos = snvRange_M.start,
             geneName = snvRange_M.mcols.geneName,
             ref = snvRange_M.mcols.ref,
             alt = snvRange_M.mcols.alt,
             ref_counts = snvRange_M.mcols.refc,
             var_counts = snvRange_M.mcols.altc,
             normal_cn,minor_cn,major_cn,
             sample_snv = mcols.sampleName,
             sample_cnv = cnvRangeA_M.sample_cnv)%>%
      unite(col = "mutation_id",chr,pos,geneName,ref,alt,sep = ":")
    
    # for Including Amino acid changes in the output
    # x = full_join(y,mutID.AA,"mutation_id")

    
    dir.create(path = "pyclone/",showWarnings = F)
    write.table(y,file = paste0("pyclone/",samples[i],".tsv"),
                sep = "\t",
                append = F,
                row.names = F)
  }
}

