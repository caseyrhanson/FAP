###for changing Ruping's output file into something more useable

library(dplyr)
library(tidyr)
data = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt"
ccf.threshold = .1
mutsigFilePrep = function(data,ccf.threshold){
  file = read.table(file = data,header = T,sep = "\t")
  str(file)
  colnames(file)
  table(file$chr)
  
  # file_select = dplyr::select(file, chr, pos, id, ref, alt, ends_with("maf"), geneLoc, functionalClass, somatic, germline, founds)
  file_select = dplyr::select(file, chr, pos, id, ref, alt,
                              ends_with("ccf"), geneLoc, functionalClass, somatic)%>%
                dplyr::select(-starts_with("bp"))
  table(file_select$chr)
  colnames(file_select)
  #for mutational signature analysis
  #make a file(".*_MutSigInput_groups_filtered.txt") with only : "Chr Pos Ref Var SampleName"
  # file_mutsig<-select(file, Chr=chr, Pos=pos, Ref=ref, Var=alt, ends_with("maf"))%>%
  #   gather("SampleName", "MeanAlleleFreq", 5:34)%>%
  #   filter(MeanAlleleFreq > 0)%>%
  #   separate(SampleName, sep="maf",c("SampleName", "maf"))%>%
  #   select(-maf,-MeanAlleleFreq)
  file_mutsig = dplyr::select(file, Chr=chr, Pos=pos, Ref=ref, Var=alt,
                              ends_with("ccf"),-starts_with("bp"),-starts_with("merge"))
  
  idx.ccf.samples = which(grepl(pattern = "ccf",x = colnames(file_mutsig)))
  first = idx.ccf.samples[1]
  last = idx.ccf.samples[length(idx.ccf.samples)]
  
  file_mutsig = tidyr::gather(file_mutsig,"SampleName", "ccf", first:last)%>%
    filter(ccf > ccf.threshold)%>%
    separate(SampleName, sep ="ccf",c("SampleName", "ccf"))%>%
    dplyr::select(-ccf)
  
  str(file_mutsig)
  table(file_mutsig$Chr)
  distinct()
  return(file_mutsig)
  # write.table(x = file_mutsig,file = "FAP_snv_MutSigInput_groups_filtered.txt",append = F,sep = "\t",col.names = T,row.names = F)
}






