vapReduce_cover_vcover_mutect = function(vap, patient,mafOrccf){
  # Make a new table with sample_mutect, sample_cover sample_freq
  # Mutect is boolean (yes or no), cover is depth (d), freq is maf (percentage).
  
  library(dplyr);library(tidyr);library(stringr)
  
  # Make vap file patient specific
  # Choose only first columns ("chr", "pos", "id", "ref", "alt") and  columns related to specific patient
  index.patient = which(grepl(patient,colnames(vap)))
  vap.reduced.patientspecific = vap[,c(1:6,index.patient)]  
  
  vap.names = colnames(vap.reduced.patientspecific)
  
  
  
  if (mafOrccf=="maf") {
    
    #List out the Sample, Tumor and Normal sample names
    sample.names = vap.names[which(grepl("maf$",vap.names))]
    sample.names = str_replace(sample.names,pattern = "maf",replacement = "")
    
    normal.names = sample.names[grep("blood",sample.names)]
    tumor.names = sample.names[grep("blood",sample.names,invert = T)]
    
    #Select specific columns to keep.
    ##Keep depth and maf in the reduced table
    vars = c("CADD_phred", "Polyphen2_HVAR_pred","mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad")
    vap.reduced = vap.reduced.patientspecific%>%select(chr, pos, id, ref, alt, ends_with("maf"),ends_with("d"),-one_of(vars = vars))
    index.maf.cols = which(grepl("maf",colnames(vap.reduced)))
    
    #determine which maf values are not "0", then label them "yes". Label "0" as "unknown"
    maf = vap.reduced[,index.maf.cols]
    maf[maf!=0] = "yes"
    maf[maf==0] = "unknown"
    maf[is.na(maf)] = "unknown"
    
    ##Create new columns that say _mutect
    mutect.names = str_c(sample.names,"_mutect")
    
    ##Create new dataframe with all the samples having a new column and each value is "yes"
    df = maf
    colnames(df) = mutect.names
    
    
    #merge the two dataframes together
    vap.reduced.and.mutect = cbind(vap.reduced,df)
    
    
    
    
    return(list(phylip.input.maf.file = vap.reduced.and.mutect,
                patient.samplenames = sample.names,
                normals = normal.names,
                tumor = tumor.names))
  }else if(mafOrccf=="ccf"){
    
    
    #List out the Sample, Tumor and Normal sample names
    sample.names = vap.names[which(grepl("ccf$",vap.names))]
    sample.names = str_replace(sample.names,pattern = "ccf",replacement = "")
    
    normal.names = sample.names[grep("blood",sample.names)]
    tumor.names = sample.names[grep("blood",sample.names,invert = T)]
    
    #Select specific columns to keep.
    ##Keep depth and maf in the reduced table
    vars = c("CADD_phred", "Polyphen2_HVAR_pred","mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad")
    vap.reduced = vap.reduced.patientspecific%>%select(chr, pos, id, ref, alt, ends_with("mafc"),
                                                       ends_with("d",ignore.case = F),-one_of(vars = vars))
    
    
    #determine which maf values are not "0", then label them "yes". Label "0" as "unknown"
    index.mafc.cols = which(grepl("mafc",colnames(vap.reduced))) #find the mafc columns
    mafc = vap.reduced[,index.mafc.cols] #separate out the mafc columns
    mafc[mafc!=0] = "yes"
    mafc[mafc==0] = "unknown"
    mafc[is.na(mafc)] = "unknown"
    
    ##Create new columns that say _mutect
    mutect.names = str_c(sample.names,"_mutect")
    
    ##Create new dataframe with all the samples having a new column and each value is "yes"
    mutect = mafc
    colnames(mutect) = mutect.names
    
    
    #collect the "depth" columns for each of the samples we have
    d.cols = grep("d$",colnames(vap.reduced),value = T)  
    d.cols.rem.d = str_sub(string = d.cols,end = -2)
    index.d.cols = d.cols[which(d.cols.rem.d %in% sample.names)]
    d = vap.reduced[,index.d.cols]
    
    
    #merge the two dataframes together
    vap.reduced.and.mutect = cbind(vap.reduced[,1:5],
                                   vap.reduced[,index.mafc.cols],
                                   vap.reduced[,index.d.cols],
                                   mutect)
    
    return(list(phylip.input.maf.file = vap.reduced.and.mutect,
                patient.samplenames = sample.names,
                normals = normal.names,
                tumor = tumor.names))
  }
}
