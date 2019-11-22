###for changing Ruping's output file into something more useable

library(dplyr)
library(tidyr)
file<-read.delim(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt",header = T,sep = "\t")
str(file)
colnames(file)
table(file$chr)

file_select<-select(file, chr, pos, id, ref, alt, ends_with("maf"), geneLoc, functionalClass, somatic, germline, founds)
table(file_select$chr)

#for mutational signature analysis
#make a file(".*_MutSigInput_groups_filtered.txt") with only : "Chr Pos Ref Var SampleName"
file_mutsig<-select(file, Chr=chr, Pos=pos, Ref=ref, Var=alt, ends_with("maf"))%>%
        gather("SampleName", "MeanAlleleFreq", 5:34)%>%
        filter(MeanAlleleFreq > 0)%>%
        separate(SampleName, sep="maf",c("SampleName", "maf"))%>%
        select(-maf,-MeanAlleleFreq)
str(file_mutsig)
table(file_mutsig$Chr)
distinct()
write.table(x = file_mutsig,file = "FAP_snv_MutSigInput_groups_filtered.txt",append = F,sep = "\t",col.names = T,row.names = F)





