##Mutational signatures from COSMIC website
##Use SNVs rescued from pileup files
##Use individual samples
#This is very helpful
# vignette("SomaticSignatures-vignette")
# https://cancer.sanger.ac.uk/cosmic/signatures

# Great Tutorial here:
# http://bioconductor.org/packages/release/bioc/vignettes/SomaticSignatures/inst/doc/SomaticSignatures-vignette.html

## input the file with Mutation Signature analysis

##all tumors
library(dplyr)
library(NMF)
library(SomaticSignatures)
library(ggplot2)
library(grDevices)
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# library(BSgenome.Hsapiens.UCSC.hg19)

# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

# library(GenomicRanges)
# library(VariantAnnotation)

mutsignature <- function(mutation,mutsig) {
  ##mutation: data frame with columns Chr, Pos, Ref, Var, SampleName
  ##mutsig: mutational signature matrix (96 x n_signatures) or vector of integers for numbers of de novo signatures

  data(signatures21, package = "SomaticSignatures")

  if ( ! grepl("chr",mutation$Chr[1])) {
    mutation$Chr <- paste("chr",mutation$Chr,sep="")
  }
  chrs <- paste("chr",c(1:22,"X","Y"),sep="")
  mutation <- mutation[mutation$Chr %in% chrs,]
  mutation_vr <- VRanges(seqnames=factor(mutation$Chr,levels=chrs),
			 ranges=IRanges(start=as.numeric(as.character(mutation$Pos)),width=1),
			 ref=mutation$Ref,
			 alt=mutation$Var,
			 sampleNames=mutation$SampleName)
  mutation_motif <- mutationContext(vr = mutation_vr, BSgenome.Hsapiens.UCSC.hg38, unify = TRUE)
  
  data(kmers)
  norms <- k3wg / k3we

  mutation_mm <- motifMatrix(mutation_motif, normalize = TRUE)
  if (nrow(mutation_mm) < 96) {
    temp <- matrix(0,nrow=96,ncol=ncol(mutation_mm))
    rownames(temp) <- rownames(signatures21)
    colnames(temp) <- colnames(mutation_mm)
    temp[rownames(mutation_mm),] <- mutation_mm
    temp <- temp + 1e-8
    mutation_mm <- temp
  }
  mutation_mm_n <- normalizeMotifs(mutation_mm,norms)
  mutation_mm_n <- mutation_mm_n/(rep(colSums(mutation_mm_n),each=nrow(mutation_mm_n)))

  ### issue
  if (is.matrix(mutsig)) {
    mutation_contri <- fcnnls(mutsig,mutation_mm_n)
    mutation_sig <- vector("list",1)
    temp <- identifySignatures(cbind(mutation_mm_n,runif(96)),2) 
    slot(temp,"signatures") <- mutsig
    slot(temp,"samples") <- t(mutation_contri$x)
    slot(temp,"observed") <- mutation_mm_n
    slot(temp,"fitted") <- mutation_contri$fitted
    slot(temp,"nSignatures") <- ncol(mutsig)
    mutation_sig[[1]] <- temp
  }

  if (is.vector(mutsig)) {
    mutation_sig <- vector("list",length(mutsig))
    for (i in 1:length(mutsig)) {
      mutation_sig[[i]] <- identifySignatures(mutation_mm_n,mutsig[i])
    }
  }
  return(mutation_sig)
}

##########
### Find modified SNV file and cosmic signature file
setwd("/Users/ahorn720/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP")
cosmic_sig = as.matrix(read.table(
  "/Users/ahorn720/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/COSMIC_30_MutSignatures.txt", 
                                  header=T, row.names=1,sep="\t"));

files <- list.files(path="/Users/ahorn720/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP",
                      pattern=".*_MutSigInput_groups_filtered.txt",
                    full.names=T, recursive=FALSE)

# files<- list.files(path="~/Downloads/",
#                    pattern = "FiveCases_MutSigInput_groups_filtered.txt",
#                    full.names = T, recursive = F)


# read four cases combined
mutation = read.table(files[1], header=T);
mutation<-mutate(mutation,Chr=as.factor(Chr))
mutation_sig = mutsignature(mutation, cosmic_sig);

#created by Katherine
#load(file = "~/Downloads/sigs_for_aaron.Rdata")
save(mutation_sig, file = "sigs_for_aaron.Rdata")

data_EP = mutation_sig[[1]]@samples[1:11,]
data_ES = mutation_sig[[1]]@samples[12:14,]
data_JC = mutation_sig[[1]]@samples[15:17,]
data_JP = mutation_sig[[1]]@samples[18:28,]
data_NBM = mutation_sig[[1]]@samples[29:30,] 




# read U55
# mutation.2 = read.table(files[3], header=T);
# mutation_sig.2 = mutsignature(mutation.2, cosmic_sig);
# data.2 = mutation_sig.2[[1]]@samples[1:3,]
# data.2 = data.2[,c(1,2,3,5,12,13,15,16,18,20,26)]
# data.2 = data.2[c(3,2,1),]


library(colorRamps)
pdf("/Users/ahorn720/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/SNV_mutSig/cosmic_sigs/SomaticSigSummary.pdf", width=7,height=6)

layout(rbind(1,2,3), heights=c(2,2,1))  # put legend on bottom 1/8th of the chart
col <- rep(matlab.like(ncol(data)),3)
par(mar=c(3, 1, 2, 1))
barplot(t(data_EP),col=col, main ="Patient EP Cosmic SNV Signatures", cex.names = .7)
barplot(t(data_JP),col=col, main ="Patient JP Cosmic SNV Signatures", cex.names = .7)
barplot(t(data_ES),col=col, main ="Patient ES Cosmic SNV Signatures", cex.names = .7)
barplot(t(data_JC),col=col, main ="Patient JC Cosmic SNV Signatures", cex.names = .7)
barplot(t(data_NBM),col=col, main ="Patient NBM Cosmic SNV Signatures", cex.names = .7)


# 
par(mar=c(3, 1, 2, 1))
# barplot(t(data_JC),col=col, main = "U55", cex.names = 1.3)

#par(mar=c(0, 0, 0, 0))

plot.new()
legend("center", legend = colnames(data), col = matlab.like(ncol(data)), cex=.7,
         bty="n", pch = 15, ncol=ncol(data), title = "COSMIC somatic signatures",
         x.intersp=0.45, pt.cex=2)

dev.off();
  
## Create Heatmap of Mutation Signatures for each Sample
setwd("/Users/ahorn720/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/SNV_mutSig/cosmic_sigs")
pdf("CosmicMutSig_bySample.pdf")
plotSampleMap(s = mutation_sig[[1]]) + 
ggtitle("Cosmic Somatic Signatures by Sample")
dev.off();
