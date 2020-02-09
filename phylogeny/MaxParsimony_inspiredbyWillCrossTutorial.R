rm(list=ls())
library(ape)
library(phangorn)

#Gather data from individual samples

pt = "A001"
source("~/aarons_FAP_github_repository/phylogeny/MaxParsimony_inspiredbyWillCrossTutorial_input.R")
dataSub = maxpars_input(pt = pt)
head(dataSub)
#remove EP.AdenoCa for now because the data from that sample is bad
dataSub = dataSub[,!grepl(pattern = "EP.AdenoCa",x = colnames(dataSub))]
colnames(dataSub) = str_replace_all(string = colnames(dataSub),pattern = "EP",replacement = "F")
#remove JP.AdenoCa for now because the data from that sample is bad
dataSub = dataSub[,!grepl(pattern = "JPAdenoCa",x = colnames(dataSub))]
colnames(dataSub) = str_replace_all(string = colnames(dataSub),pattern = "JP",replacement = "G")

str(dataSub)
head(dataSub)
#make distance matrix
njInput <- dist(as.matrix(t(dataSub)), method = "euc")

#perform nj analysis
njPhy <- nj(njInput)

#set normal as root
njPhy <- root(njPhy, outgroup = "Normal")

#plot data
# layout(1)
# layout.show()
yarrr::piratepal(palette = "info2")
plot.phylo(njPhy, use.edge.length = TRUE,show.node.label = T,
           edge.color = yarrr::piratepal(palette = "info2"))
title("Neighbor Joining Tree")

###### for Maximum Parsimony tree building #######

#convert data to phyDat form
convData <- as.phyDat(as.matrix(t(dataSub)), type = "USER", levels=c(0,1),site.pattern = F)
#perform parsimony analysis
prachSave <- pratchet(data = convData, method = "fitch",
                      start = NULL, maxit = 1000000, all = TRUE, trace = 1)
#Root the tree
if (class(prachSave)=="multiPhylo") {
  prachSave <- root.phylo(prachSave[[1]], outgroup = "Normal")
} else {
  prachSave <- root.phylo(prachSave, outgroup = "Normal")
}

# #get parsimony score for NJ tree
# parsScore <- parsimony(tree = njPhy, data = convData)
# #get HI value: homoplasy index
# HIindex <- 1 - (nrow(dataIn) / parsScore)

#Build tree. Adjust edge lengths using acctran
tree.acctran = acctran(tree = prachSave,data = convData)
tree.acctran.parscore = parsimony(tree.acctran,convData)
# The Consistency Index isthe Minimum number of steps given your dataset divided by the actual number on your tree (observed steps).
# (so CI=1 means you have no homoplasy)
ci = CI(tree.acctran,convData)
# Homoplasy index which is 1-CI (so HI=0 means no homoplasy
# You might want that your index increases with increasinghomoplasy)
HIindex = 1-ci    #From Will: 1 - (nrow(dataIn)/tree.acctran.parscore)


#Plot Max Parsimony Tree without Bootstrapping
plot.phylo(tree.acctran,main = paste0(pt,": Max Parsimony Score ",tree.acctran.parscore,
                                      " with pratchet HI Index ",
                                      round(HIindex,3)),
     use.edge.length = TRUE,xpd = NA)
edgelabels(text = round(tree.acctran$edge.length,digits = 0),frame = "none",adj = c(0,-.5),cex = .8)
add.scale.bar(lcol = "purple")

write.tree(phy = root.phylo(tree.acctran,"Normal"),file = paste0("~/Bulk_A001_A002/phylogeny/",pt,"MaxPars.nwk"))


########## With Bootstrapping
set.seed(123)
BStrees <- bootstrap.phyDat(x = convData,FUN = pratchet, bs = 100)
BStrees = root.multiPhylo(BStrees,outgroup = "Normal")
BStrees
class(BStrees)
class(tree.acctran)
# The Consistency Index isthe Minimum number of steps given your dataset divided by the actual number on your tree (observed steps).
# (so CI=1 means you have no homoplasy)
ci = CI(tree.acctran,convData)
# Homoplasy index which is 1-CI (so HI=0 means no homoplasy
# You might want that your index increases with increasinghomoplasy)
HIindex = 1-ci    #From Will: 1 - (nrow(dataIn)/tree.acctran.parscore)

tree.acctran.parscore = parsimony(tree.acctran,convData)

svg(file = paste0("~/Bulk_A001_A002/phylogeny/",pt,"MaxPars_wBootStrap.svg"))
x = plotBS(tree.acctran, BStrees, "phylogram",p = 0,bs.col = "darkred",bs.adj = c(0,1),
       main = paste0(pt,": Max Parsimony Score ",tree.acctran.parscore,
                     " with pratchet HI Index ",round(HIindex,3)),
       xpd = NA)
plot.phylo(tree.acctran,
           main = paste0(pt,": Max Parsimony Score ",tree.acctran.parscore,
                                      " with pratchet HI Index ",round(HIindex,3)))
nodelabels(x$node.label,frame = "none",adj = -.5,cex = .7,col = "darkred")
add.scale.bar(lcol = "purple")
edgelabels(text = round(tree.acctran$edge.length,digits = 0),
           frame = "none",cex = .8,adj = c(0,-0.5))
dev.off()


# Assignment of somatic changes in WGS to the phylogenetic tree
# Substitutions
# # Substitutions were called as present or absent in each organoid 
# as described above. To assign these mutations to the tree, each branch 
# of the tree was considered in turn. If a mutation was called in all the 
# organoids that were descendants of a given branch, and in no organoids 
# that were not descendants of the branch, mutations were assigned to that 
# branch. Ignoring private mutations, which necessarily fit any tree, 
# 97.7% of shared mutations fitted the tree structure from patient 1 
# perfectly, 89.7% fitted the tree from patient 2 perfectly, and 88.1% fitted 
# the tree from patient 3 perfectly. The lower concordance with the tree for 
# patients 2 and 3 reflects the increased copy number changes that have occurred in 
# these phylogenies. Examination of the copy number state at loci where there were 
# discordant mutations showed that the majority could be explained by deletions
# of those mutations in a subclone. Substitutions that did not fit the tree
# perfectly were therefore assigned to the most recent common ancestor of
# the samples in which they were called.
