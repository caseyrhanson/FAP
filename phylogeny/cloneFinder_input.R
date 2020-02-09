# Example input
# SNVID	Wild	Mut	R2:ref	R2:alt	R3:ref	R3:alt	R4:ref	R4:alt	
# S1	T	A	70	30	70	30	50	50	
# S2	C	A	70	30	70	30	50	50
# S3	A	T	70	30	70	30	50	50
# S4	A	T	70	30	70	30	50	50
# S5	C	T	70	30	70	30	40  60CPCOLS <- c("#1f78b4", "#33a02c", "#e31a1c", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")

rm(list=ls())

library(stringr);library(GenomicRanges);library(IRanges);library(dplyr);library(tidyr)

setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")

############ Prep for Treeomics

pts=c("A001","A002")
pt="A002"
file.type = "ccf"


sampAB = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"), header = T, sep = "\t")
sampAB$Wild = sampAB$ref #make new column just renamed
sampAB$Mut =  sampAB$alt #make new column just renamed
#select refc and altd columns and unite the SNVID column
snvs<-dplyr::select(sampAB,
             Chromosome = chr,
             Position = pos,id,
             Wild, Mut,
             ref, alt,
             ends_with("refc"), ends_with("altc"),-CADD_phred,-Polyphen2_HVAR_pred,-ends_with("ccfSD"),-starts_with("bp"))%>%
  unite(col = "SNVID",c("Chromosome","Position", "id", "ref", "alt" ),sep = ":")

#rename and reorder the sample columns, collect only samples with both refc and altc columns types
snvids.index = which(!grepl(pattern = pt,x = colnames(snvs))) #find the columns which dont have the pt name in it
pt.index = which(grepl(pattern = pt,x = colnames(snvs))) #find the columns with sample names in it
pt.samples = sort(colnames(snvs)[pt.index]) #sorted list of the patient columns: A001C004altc A001C005altc ... A001C219ref
pt.samples = str_remove_all(string = pt.samples,pattern = "refc") #remove refc
pt.samples = str_remove_all(string = pt.samples,pattern = "altc") # remove altc
pt.samples = pt.samples[duplicated(pt.samples)] #find paired refc and altc columns
pt.samples.refc = paste0(pt.samples,"refc") #create a list of sample names with refc after it
pt.samples.altc = paste0(pt.samples,"altc") #create a list of sample names with altc after it
pt.samples = c(pt.samples.refc,pt.samples.altc) # creaes new list of samples with refc first then altc


snvs = cbind(snvs[,snvids.index], snvs[,pt.samples]) #put all the SNVID and refc and altc columns together
snvids.index = which(!grepl(pattern = pt,x = colnames(snvs)))
pt.index = which(grepl(pattern = pt,x = colnames(snvs)))
pt.samples = colnames(snvs)[pt.index] #collect all sample names 
pt.samples = str_replace_all(string = pt.samples,pattern = "refc",replacement = ":ref") #change the column labels
pt.samples = str_replace_all(string = pt.samples,pattern = "altc",replacement = ":alt") #change the column labels
colnames(snvs)[pt.index] = pt.samples

head(snvs)

# To keep track of the mutations, re-add the columns with the ID infro
snvs$SNVID_2 = snvs$SNVID # Remake the mutID column
snvs.small = snvs%>%separate(col = SNVID_2,into = c("Chromosome","Position", "id", "ref", "alt"),sep = ":")
#for the snv/CNV filtering we need this column to be numeric
snvs.small$Position = as.numeric(as.character(snvs.small$Position)) 
head(snvs.small)
dim(snvs.small)


################# Filter out the SNVs within CNVs, so we only collect HET snps for the clonefinder
# Make a Granges object for the SNVs for filtering
snvs.ranges = GRanges(seqnames = snvs.small$Chromosome,
                      ranges = IRanges(start = snvs.small$Position,width = 1),
                      strand = NULL,mcols = snvs.small)

# collect all the segment files for the pt
pt.segments = list.files(path = "titan/segments",pattern = pt,full.names = T)
if(pt == "A002"){
  #For A002, only consider samles 102 and 202 for their CNV aberrations.
  # Their plots are the cleanest and I think the most reliable visually
  pt.segments = pt.segments[grep(pattern = "102|202",pt.segments)]
}

datalist = list()
for (seg in 1:length(pt.segments)) {
  file = read.delim(file = pt.segments[seg])
  datalist[[seg]] = file
}
pt.all.segments = do.call(rbind,datalist) #merge all the segments together
head(pt.all.segments)

#filter for HET regions only and select only location columns and arrange in order
pt.all.segments.het = pt.all.segments%>% 
  filter(LOHcall=="HET")%>%
  dplyr::select(chrom,loc.start,loc.end)%>%
  arrange(loc.start)%>%
  arrange(chrom)
head(pt.all.segments.het)

# Create segments GRanges object
segments.het.ranges = GRanges(seqnames = pt.all.segments.het$chrom,
                              ranges = IRanges(start= pt.all.segments.het$loc.start,
                                               end= pt.all.segments.het$loc.end,
                                               names = paste0("HET",1:nrow(pt.all.segments.het))),
                              strand = NULL)
head(segments.het.ranges)
#Filter for Non-HET regions only and select only location columns and arrange in order
pt.all.segments.nonhet = pt.all.segments%>% 
  filter(LOHcall!="HET")%>%
  dplyr::select(chrom,loc.start,loc.end)%>%
  arrange(loc.start)%>%
  arrange(chrom)
head(pt.all.segments.nonhet)

# Create segments GRanges object for Non-HET segments
segments.nonhet.ranges = GRanges(seqnames = pt.all.segments.nonhet$chrom,
                              ranges = IRanges(start= pt.all.segments.nonhet$loc.start,
                                               end= pt.all.segments.nonhet$loc.end,
                                               names = paste0("nonHET",1:nrow(pt.all.segments.nonhet))),
                              strand = NULL)
segments.nonhet.ranges
#Find segments which are HET in every sample. Should be no non-HET samples 
het.outside.nonhet = subsetByOverlaps(segments.het.ranges,segments.nonhet.ranges,
                                      invert = T) #invert ensures we select het segments outside non-het segments
segments.het.ranges = het.outside.nonhet

# x = segments.het.ranges[1,]
# y = subsetByOverlaps(snvs.ranges,x)

# Find Overlapping SNVs within HET segments
snvs.in.hets = subsetByOverlaps(snvs.ranges, segments.het.ranges)
snvs.in.hets
length(snvs.in.hets)
snvs = as.data.frame(snvs.in.hets) #turn the GRanges object into a dataframe for further downstream prep.


snvs = snvs[,grepl(pattern = "mcols",x = colnames(snvs))] #the "mcol" containing columns are the ones we want to work with
remove = c("mcols.Chromosome", "mcols.Position","mcols.id","mcols.ref","mcols.alt") #remove these columns
snvs = snvs%>%
  dplyr::select(-one_of(remove))
colnames(snvs)
colnames(snvs) = str_remove_all(string = colnames(snvs),pattern = "mcols.") # Remove "mcols" from colnames and replace the "." with a ":"
colnames(snvs) = str_replace_all(string = colnames(snvs),pattern = "\\.",replacement = ":") # Replace the "." with a ":"

head(snvs)

#sort the columns so that ":ref" follows ":alt"
pt.samples = sort(colnames(snvs)[grepl(pattern = pt,x = colnames(snvs))]) #sort the column names by pt sample
pt.samples = str_replace_all(pt.samples,pattern = "ref",replacement = "x")
pt.samples = str_replace_all(pt.samples,pattern = "alt", replacement = "ref")
pt.samples = str_replace_all(pt.samples,pattern = "x", replacement = "alt")
snvs.SNVID.Wild.Mut = snvs[,!grepl(pattern = pt,x = colnames(snvs))] #the SNVID, Wild, Mut columns
snvs.samples = snvs[,pt.samples] # the ref and alt columns

snvs = cbind(snvs.SNVID.Wild.Mut,snvs.samples) #column bind all SNVID, Wild, Mut and ref and alt columns
str(snvs)
write.table(x = snvs,file = paste0("clonefinder/",pt,"_HETsegmentswoNonHets_w_mutIDasSNVID.txt"),
            append = F,sep = "\t",row.names = F,quote = F)

# # python clonefinder.py snv path/Input.txt
path.to.input.file = paste0("~",getwd(),
                            "/",list.files(path = "clonefinder",
                                           pattern = paste0(pt,"_HETsegmentswoNonHets_w_mutIDasSNVID.txt"),
                                           full.names = T))
paste0("python clonefinder.py snv ",path.to.input.file)

################################ Stop here and run CloneFinder ##############

#### Or filter out the low read count samples

# #Make new SNVID column with S Number
# snvs$mutID = snvs$SNVID
# snvs$SNVID = paste0("S",1:nrow(snvs))

head(snvs)

# all sample names in this table
samples = str_remove_all(string = colnames(snvs)[str_detect(string = colnames(snvs),pattern = "ref")], pattern = ":ref")
samples
############################# Filter out low depth ##############
d = matrix(data = 0,nrow = nrow(snvs),ncol = length(samples)) #make matrix with dimensions: rows of snvs, columns of samples
colnames(d) = paste0(samples,"_d")

for (i in 1:nrow(snvs)) {
  for (j in 1:length(samples)) {
    depth = snvs[i,paste0(samples[j],":ref")] + snvs[i,paste0(samples[j],":alt")]
    d[i,paste0(samples[j],"_d")] = depth
  }
}

summary(d) #determined that a total depth of 20 was around the 1st quartile for each sample. Will remove SNVs if ANY sample has less than 20 read depth
cut = d < 20 #TRUE means d is less than 20, should be cut
cut[1:10,]
rows.to.cut = apply(cut,1,any) #Find SNV rows in which any depth is less than 20
table(rows.to.cut)

snvs = snvs[!rows.to.cut,]

################ Filter low Alternate Read counts

colnames(snvs)
samples.alt = paste0(samples,":alt")
a = snvs[,samples.alt]
a.low = a <= 3 #which alt reads are less than 3
cut.low.alt = apply(a.low,1,all) #which rows have all values less than 3.
sum(cut.low.alt)

snvs = snvs[!cut.low.alt,]
head(snvs)

write.table(x = snvs,
            file = paste0("clonefinder/",pt,"_HETsegmentsonly_20dCut_3altCut_HETsegmentswoNonHets.txt"),
            append = F,sep = "\t",row.names = F,quote = F)


####### Change sample names to Rs
# samples = unique(str_remove_all(string = pt.samples,pattern = ":(ref|alt)$"))
# pt.Rs = paste0("R",1:length(samples)) #create an R# value for every A001CXXX sample
# for (i in 1:length(samples)) {
#   for (n in 1:length(colnames(snvs))) {
#     if (grepl(pattern = samples[i],x = colnames(snvs)[n])) {
#       colnames(snvs)[n] = str_replace(string = colnames(snvs)[n],pattern = samples[i],replacement = pt.Rs[i])
#     } else {NA}
#   }
# }
# 
# colnames(snvs)
# 

#remove snvid column with S numbers. Replace with mutID column
snvs = snvs%>%select(-SNVID)%>%rename_("SNVID" = "mutID")
not.snvid.ids = which(!grepl(pattern = "SNVID",x = colnames(snvs)))
snvid.ids = which(grepl(pattern = "SNVID",x = colnames(snvs)))
snvs = cbind(SNVID = snvs[,snvid.ids],snvs[,not.snvid.ids])
         
write.table(x = snvs,file = paste0("clonefinder/",pt,".txt"),
            append = F,sep = "\t",row.names = F,quote = F)
# RtoSample = cbind(pt.Rs,samples)
# write.table(x = RtoSample,file = "clonefinder/Rs-to-Samplesl.txt",append = F,quote = F,sep = "\t")
# # python clonefinder.py snv path/Input.txt
# path.to.input.file = paste0(getwd(),"/",list.files(path = "clonefinder",pattern = "samplesRenamed",full.names = T))
# paste0("python clonefinder.py snv ",path.to.input.file)
