https://github.com/im3sanger/dndscv
rm(list=ls())
library(dndscv)# library(devtools); install_github("im3sanger/dndscv")

source("~/aarons_FAP_github_repository/dnds/dnds_input.R")

# data("dataset_simbreast", package="dndscv") # Example dataset
# head(mutations)

#read in mutation file so it looks like:
#   sampleID chr      pos ref mut
# 1 Sample_1   1   871244   G   C
# 2 Sample_1   1  6648841   C   G

############ All Patients together.... Start here ###########
# For gathering all patients together into 1 list/table
F = dnds_input(pt = "EP")
table(F$sampleID) #makes sure the samples are included after running the data rearranging program
table(F$functionalClass)
G = dnds_input(pt = "JP")
table(G$sampleID)
A001 = dnds_input(pt = "A001")
table(A001$sampleID)
A002 = dnds_input(pt = "A002")
table(A002$sampleID)

mutations.wAnnotations = rbind(F,G,A001,A002)
table(mutations.wAnnotations$sampleID)

keep.index = !(mutations.wAnnotations$sampleID %in% c("JPAdenoCa", "EP.AdenoCa"))
mutations.wAnnotations = mutations.wAnnotations[keep.index,]
table(mutations.wAnnotations$sampleID)

#add chr to columns that dont have it already
chr.index = which(!grepl(pattern = "chr",x = mutations.wAnnotations$chr))
mutations.wAnnotations$chr[chr.index] = paste0("chr",mutations.wAnnotations$chr[chr.index])
table(mutations.wAnnotations$chr)

#removes additional annotation columns
mutations = select(mutations.wAnnotations,-functionalClass,-geneName,-presence)

################### #For a single patient sample.... start here  ##########
pt = "A002"
mutations = dnds_input(pt = pt)
table(mutations$chr)
table(mutations$sampleID)
#add chr to columns that dont have it already
chr.index = which(!grepl(pattern = "chr",x = mutations$chr))
mutations$chr[chr.index] = paste0("chr",mutations$chr[chr.index])

#remove the JP and EP adenocarcinoma samples
keep.index = !(mutations$sampleID %in% c("JPAdenoCa", "EP.AdenoCa"))
mutations = mutations[keep.index,]

mutations = select(mutations,-functionalClass,-geneName,-presence)
head(mutations)
table(mutations$sampleID)

# mutations = mutations[mutations$chr=="3" & mutations$pos>segment[1] & mutations$pos<segment[2], ] # Restricting the mutations to those inside the segment
# mutations$pos = mutations$pos-segment[1]+1 # Correcting the position of the mutations to their position in the reference fasta file used

# load("~/aarons_FAP_github_repository/dnds/RefCDS_human_GRCh38.p12.rda")

#Run dN/dS analysis on the mutations
dndsout = dndscv(mutations, refdb="~/aarons_FAP_github_repository/dnds/RefCDS_human_GRCh38.p12.rda", cv=NULL)
sel_cv = dndsout$sel_cv
sel_cv = sel_cv[sel_cv$pallsubs_cv <= 0.1,]
sel_cv = sel_cv[order(sel_cv$pallsubs_cv),]

head(sel_cv)

write.table(x = sel_cv,file = paste0("~/Bulk_A001_A002/dnds/","dndsout_",pt,".txt"),
            append = FALSE,sep = "\t",row.names = FALSE)
# write.table(x = sel_cv,file = paste0("~/Bulk_A001_A002/dnds/","dndsout_","EPJPA001A002",".txt"),
#             append = F,sep = "\t",row.names = F)

head(dndsout$sel_cv)
print(dndsout$globaldnds)
print(dndsout$nbreg$theta)
head(dndsout$annotmuts)

# This is shown as an example but these results based on a few genes should not be trusted
# The output of the dndscv function is a list of objects. 
# For an analysis of exome or genome data, the most relevant output will often be the result
# of neutrality tests at gene level. P-values for substitutions are obtained by Likelihood-Ratio
# Tests as described in (Martincorena et al, 2017) and q-values are obtained by 
# Benjamini-Hodgberg’s multiple testing correction. The table also includes information on the
# number of substitutions of each class observed in each gene, as well as maximum-likelihood 
# estimates (MLEs) of the dN/dS ratios for each gene, for missense (wmis), nonsense (wnon), 
# essential splice site mutations (wspl) and indels (wind). 
# The global q-value integrating all mutation types are available in the qglobal_cv and 
# qallsubs_cv columns for analyses with and without indels, respectively.


#local dn/ds: more classical test
signif_genes_localmodel = as.vector(dndsout$sel_loc$gene_name[dndsout$sel_loc$pall_loc<0.1])
print(signif_genes_localmodel)
head(dndsout$sel_loc[dndsout$sel_loc$pall_loc<=0.1,])
dndsout.loc = dndsout$sel_loc[dndsout$sel_loc$pall_loc<=0.1,]

write.table(dndsout.loc,
            file = paste0("~/Bulk_A001_A002/dnds/","dndsout.loc_",pt,".txt"),
            append = FALSE,sep = "\t",row.names = FALSE)

# write.table(dndsout.loc,
#             file = paste0("~/Bulk_A001_A002/dnds/","dndsout.loc_","EPJPA001A002",".txt"),
#             append = F,sep = "\t",row.names = F)

