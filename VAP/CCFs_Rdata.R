rm(list=ls())

load("~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")
library(tidyr);library(dplyr)

# For gathering other columns from preCCF adjustment VAP table
path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table"
EPJP = read.delim(path,header = T,sep = "\t")
class(EPJP$AAChange)
EPJP.reduced = EPJP%>%
  unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")%>%
  select(mutID,AAChange)%>%
  mutate(mutID = paste0("chr",mutID))


# Patient F
# Add AAChange Column
path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.EP.txt"
ptF = read.delim(file = path,  header = T, sep = "\t")
ptF.mutID = ptF%>%unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")
ptF = left_join(ptF.mutID,EPJP.reduced,"mutID")%>%
  separate(col = mutID,into = c("chr","pos","id","ref","alt"),sep = ":")
head(ptF)

# Patient G
path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.JP.txt"
ptG = read.delim(file = path,  header = T, sep = "\t")
ptG.mutID = ptG%>%
  unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")
ptG = left_join(ptG.mutID,EPJP.reduced,"mutID")%>%
  separate(col = mutID,into = c("chr","pos","id","ref","alt"),sep = ":")
head(ptG)

# For Gathering other columns (ie AAChange) from the preCCF adjustment column
path = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table"
A1A2 = read.delim(path,header = T,sep = "\t")
A1A2.reduced = A1A2%>%
  unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")%>%
  select(mutID,AAChange)%>%
  mutate(mutID = paste0("chr",mutID))

# Patient A001
path = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt"
ptA001 = read.delim(file = path, header = T, sep = "\t")
ptA001.mutID = ptA001%>%
  unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")%>%
  mutate(mutID = paste0("chr",mutID))
ptA001 = left_join(ptA001.mutID,A1A2.reduced,"mutID")%>%
  separate(col = mutID,into = c("chr","pos","id","ref","alt"),sep = ":")
head(ptA001)

# Patient A002
path = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A002.txt"
ptA002 = read.delim(file = path,header = T,sep = "\t")
ptA002.mutID = ptA002%>%
  unite(col = "mutID","chr", "pos", "id","ref","alt",sep = ":")%>%
  mutate(mutID = paste0("chr",mutID))

ptA002 = left_join(ptA002.mutID,A1A2.reduced,"mutID")%>%
  separate(col = mutID,into = c("chr","pos","id","ref","alt"),sep = ":")
head(ptA002$AAChange)

save(ptA001, ptA002, ptF, ptG, file = "~/aarons_FAP_github_repository/VAP/FAPccfDatasets.Rdata")


