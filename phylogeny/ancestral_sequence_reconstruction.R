rm(list=ls())
library(phangorn)

data(Laurasiatherian)

dm <- dist.ml(Laurasiatherian)
class(dm)

tree <- NJ(dm)

plot(tree)

x = pratchet(as.phyDat(data("Laurasiatherian")),start = tree)
x$edge
x$edge.length
plot(x)
fit <- pml(tree, Laurasiatherian)
anc.ml <- ancestral.pml(fit, type = "ml")
anc.p <- ancestral.pars(tree, Laurasiatherian)
plotAnc(tree = tree,data = anc.p,1)
## Not run:
require(seqLogo)
seqLogo( t(subset(anc.ml, 48, 1:20)[[1]]), ic.scale=FALSE)
seqLogo( t(subset(anc.p, 48, 1:20)[[1]]), ic.scale=FALSE)
## End(Not run)
# plot the first site pattern
plotAnc(tree, anc.ml, 1)
# plot the third character
plotAnc(tree, anc.ml, attr(anc.ml, "index")[3])

anc.ml$Possum
tree$edge
tree$tip.label
tree$Nnode
anc.ml$`48`
anc.p$`48`


vignette("Ancestral")




source("~/aarons_FAP_github_repository/phylogeny/MaxParsimony_inspiredbyWillCrossTutorial_input.R")
data = maxpars_input(pt = "A001")
convData <- as.phyDat(as.matrix(t(data)), type = "USER", levels=c(0,1))
#perform parsimony analysis
prachSave <- pratchet(data = convData, method = "fitch",
                      start = NULL, maxit = 1000000, all = TRUE, trace = 1)
prachSave <- root.phylo(prachSave, outgroup = "Normal")

#Build tree. Adjust edge lengths using acctran
tree.acctran = acctran(tree = prachSave,data = convData)

asr = ancestral.pars(tree.acctran,convData)

plot(tree.acctran)
tree.acctran$edge
attr(asr,"index")
getMRCA(phy = tree.acctran,tip = c(7,3))
tree.acctran$tip.label
edgelabels()
nodelabels()
