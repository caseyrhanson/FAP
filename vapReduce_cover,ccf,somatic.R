vapReduce_cover_ccf_somatic = function(data, patient,mafOrccf){
  # Make a new table with sample_mutect, sample_cover sample_freq
  # Mutect is boolean (yes or no), cover is depth (d), freq is maf (percentage).
  
  library(dplyr);library(tidyr);library(stringr)
  
  # Make vap file patient specific
  # Choose only first columns ("chr", "pos", "id", "ref", "alt") and  columns related to specific patient
  index.patient = which(grepl(patient,colnames(data)))
  index.mutID = which(colnames(data) %in% c("chr", "pos","id" ,  "ref", "alt"))
  index.geneName = which(colnames(data) %in% c("geneName"))
  data.reduced.patientspecific = cbind(data[,c(index.mutID,index.geneName,index.patient)],somatic = data[,"somatic"])
  data.names = colnames(data.reduced.patientspecific)
  
  
  

    
    
    #List out the Sample, Tumor and Normal sample names
    sample.names = data.names[which(grepl("maf$",data.names))]
    sample.names = str_replace(sample.names,pattern = "maf",replacement = "")
    
    normal.names = sample.names[grep("blood",sample.names)]
    tumor.names = sample.names[grep("blood|EP.107.NL|JP_Dec_NL.1",sample.names,invert = T)]
    
    #Select specific columns to keep.
    ##Keep depth and maf in the reduced table
    vars = c("CADD_phred", "Polyphen2_HVAR_pred","mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad")
    data.reduced = data.reduced.patientspecific%>%
      select(chr, pos, id, ref, alt, geneName,
             ends_with("ccf"), ends_with("d",ignore.case = FALSE),
             -one_of(vars = vars), somatic)
    
    head(data.reduced)
    
    # somatic mutations chosen from Somatic Column created new Somatic matrix
    somatic.mat = matrix("",nrow = nrow(data.reduced),ncol = length(sample.names))
    colnames(somatic.mat) = sample.names
    
    #split up the somatic column by the "comma"
    somatic.row = str_split(string = data.reduced$somatic,pattern = ",")
    #removes the additional space after the last comma
    somatic.row = lapply(somatic.row,function(x) head(x,-1))
    #remove Samples that have "doubt" in them
    somatic.row = lapply(somatic.row, function(x) x[str_detect(string = x,pattern = "doubt", negate = T)])
    #remove the additional words in between the brackets next to the sample name
    somatic.row.good.goodSub = lapply(somatic.row,function(x) str_remove_all(string = x,pattern = "\\[[a-z]*[A-Z]*[a-z]*\\]"))
    
    if (patient %in% c("EP","JP")) {
      somatic.row.good.goodSub = lapply(somatic.row.good.goodSub, function(x) str_replace_all(string = x,pattern = "\\-",replacement = "\\."))
    } else {NULL}
    
    #create somatic vector for For Loop
    # somatic = vap.reduced$somatic
    
    for (row in 1:nrow(data.reduced)) {
      #For every somatic row, if the sample is present in the Somatic column, assign "present" to the cell
      somatic.mat[row,sample.names %in% somatic.row.good.goodSub[[row]]] = "present"
      print(row)
    }
    
    # create table with Present or Absent Somatic mutations
    head(somatic.mat)
    somatic.mat[somatic.mat != "present"] = "absent"
    head(somatic.mat)
    colnames(somatic.mat) = paste0(colnames(somatic.mat),"_somatic")
    
    
    #determine which columns have ccf values
    index.ccf.cols = which(grepl("ccf",colnames(data.reduced))) #find the mafc columns
    ccf = data.reduced[,index.ccf.cols] #separate out the mafc columns
    colnames(ccf) = str_replace_all(string = colnames(ccf),pattern = "ccf",replacement = "_ccf")
    head(ccf)
    
    #collect the "depth" columns for each of the samples we have
    d.cols = grep("d$",colnames(data.reduced),value = T)  
    d.cols.rem.d = str_sub(string = d.cols,end = -2)
    index.d.cols = d.cols[which(d.cols.rem.d %in% sample.names)]
    d = data.reduced[,index.d.cols]
    head(d)
    
    # Create ccf-adjusted Alternate Read Counts
    # Create empty numeric matrix
    ccf.adjusted.altc = matrix(data = 0,nrow = nrow(d),ncol = length(tumor.names))
    colnames(ccf.adjusted.altc) = paste0(tumor.names,"_ccf_adjusted_altc")
    head(ccf.adjusted.altc)
    # For every "Tumor" Sample, gather the Depth from d, and the ccf from ccf
    # then multiply them together and store the value in ccf.adjusted.altc
    i=1
    for (i in 1:length(tumor.names)) {
      temp.d = d[,paste0(tumor.names[i],"d")] #depth of "Tumor" Sample
      # head(temp.d,30)
      temp.ccf = ccf[,paste0(tumor.names[i],"_ccf")]
      # the CCF Adjusted Altc is Half the CCF (of the VAF in diploid regions) times the Depth.
      temp.ccf.adjusted.altc = round(temp.d * (temp.ccf/2))
      ccf.adjusted.altc[,paste0(tumor.names[i],"_ccf_adjusted_altc")] = temp.ccf.adjusted.altc
      # head(ccf.adjusted.altc[,1],30)
    }
    head(ccf.adjusted.altc)
    
    #merge the two dataframes together
    #data.reduced columns of chr, pos, id, ref, alt and geneName [1:6]
    data.reduced.and.somatic = cbind(data.reduced[,1:6],
                                   ccf,
                                   d,
                                   ccf.adjusted.altc,
                                   somatic.mat)
    head(data.reduced.and.somatic)
    
    return(list(data.reduced.and.somatic,
                patient.samplenames = sample.names,
                normals = normal.names,
                tumor = tumor.names))
}
