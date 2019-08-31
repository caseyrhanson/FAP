##Combine multi-sample data
##One point for each tumor site/sample
##Include indel calls from Strelka
##Use more stringent criteria to include genes
##Require a SNV to be in all samples for one tumor site and have AF >= 0.1 in at least one such sample

##Shorter list of genes
##Add ploidy barplots underneath the CNA panel
rm(list=ls())
##Add LOH annotations for APC, TP53 and SMAD4
library(GenVisR);library(dplyr);library(tidyr);library(stringr);library(BiocGenerics);library(stringr)
# setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/WES_Nov2018/")
setwd("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/Bulk_A001_A002/")

mycolors <- c("#D55E00","#0072B2","#66c2a5","#E69F00","#e78ac3","#8da0cb",
	      "#999999","#FFFFFF","#f4cae4")
plot(1:length(mycolors),1:length(mycolors),col = mycolors,pch = 17)

# dir1 <- "dir1"
# samples <- read.table("alltumorlabels2.txt",as.is=TRUE, header=T)$Bam
# samplelabels <- read.table("alltumorlabels2.txt",as.is=TRUE, header=T)$Label

samples = read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.clonal0.6.EP.txt",header = T,sep = "\t")
samples = read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified",header = T,sep = "\t")



# SNVtypes <- c("synonymous","nonsynonymous","splicing","stopgain","stoploss")
# SNVtypes2 <- c("nonsynonymous","splicing","stopgain","stoploss")
# tissuetypes <- c("P","LN","LI","LU","BM")

SNVtypes = levels(unique(samples$functionalClass)) #creates chr string of types of SNVs.
#SNVtypes2 <- c("nonsynonymous","splicing","stopgain","stoploss")
tissuetypes <-  c("Polyp","Normal","AdenoCarcinoma") # is.normal = str_detect(samplelabels,"Asc|Trans")

# patients <- c("EP","JP")
patients <- c("A001","A002")

# idx1 <- which(sapply(strsplit(samples,"_",fixed=TRUE),"[",1) %in% patients)
# samples <- samples[idx1]
# samplelabels <- samplelabels[idx1]

all_snvs <- vector("list",length(patients))
names(all_snvs) <- patients
all_snvs2 <- all_snvs
i=1
for (i in 1:length(patients)) {
  patient1 <- patients[i]
  fn1 = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.clonal0.6.",patients[i],".txt")
  # fn1 <- paste(dir1,patient1,"_MuTectSNV_Indel_Coding_Filtered.txt",sep="")
  #patient1 <- sub("mCRCTB","TB",patient1)
  #patient1 <- sub("Case","",patient1)
  count1 <- read.delim(fn1,as.is=TRUE)
  samplelabels = colnames(count1)[grepl("ccf$",colnames(count1))]
  samplelabels = str_remove(samplelabels,"ccf")
  # count1 <- count1[count1$CtoT_filter == "Pass",]
  # freqs <- count1[,grep("freq",names(count1))[-1]]
  freqs <- count1[,grep("ccf$",names(count1))] #select all the columns that have ccf at the end of the name
  # count1 <- count1[apply(freqs,1,max)>=0.1,] #Variant allele frequency larger than 0.1
  count1 <- count1[apply(freqs,1,max)>=0.1,] #select rows with maximum VAF larger than 0.1
  rownames(count1) <- NULL
  # normal <- sub("_cover","",grep("cover",names(count1),value=TRUE)[1])
  # normal <- sub("d$","",grep("d$",names(count1),value=TRUE)[1])
  # sns <- sub("_mutect","",grep("mutect",names(count1),value=TRUE))
  # sns <- sub("c$","",grep("c$",names(count1),value=TRUE))
  # sns2 <- samplelabels[match(sns,samples)]


  # freqs <- count1[,grep("freq",names(count1))[-1]]
  freqs <- count1[,grep("ccf$",names(count1))]
  cover <- count1[,grep("altc$",names(count1))]
  # cover = grep(patients[i],grep("d$",names(count1), value = T),value = T) #all the depth/coverage
  reads <- round(freqs*cover)
  present <- (freqs >= 0.05 & reads >= 3) | freqs >= 0.1

  sns2 = samplelabels #add tissue type information
  sns = sns2

  #which samples are Normal
  idxTis_norm = grep("Asc|Trans|Desc",sns2)
  sns2[idxTis_norm] = paste0(sns2[idxTis_norm],"_Normal")
  #which samples are Cancer
  idxTis_adeno = grep("AdenoCarcinoma",sns2,)
  sns2[idxTis_adeno] = paste0(sns2[idxTis_adeno],"_adeno")
  #which samples are polyps
  idxTis_polyp = grep("Asc|Trans|Desc|Adeno",sns2,invert = T)
  sns2[idxTis_polyp] = paste0(sns2[idxTis_polyp],"_Polyp")

  tissues <- sapply(strsplit(sns2,"_",fixed=TRUE),"[",2)
  #tissues[! sns %in% samples] <- NA
  tissue_types <- tissuetypes[tissuetypes %in% tissues]
  count2 <- count1[,1:5] #make a little table with mutation specific info
  count2 <- data.frame(count2,Type = count1$functionalClass, Gene = count1$geneName) #add column Type (functionalClass) ######## Gene

  #tissue1 = tissue_types[1] #just for for/loop practice
  tissue1 = "Normal"
  for (tissue1 in tissue_types) {
    idx1 <- which(tissues == tissue1)
    freq1 <-freqs[,idx1,drop=FALSE]
    tissue_present <- rep("no",nrow(count1)) #set everything to no first
    tissue_present[apply(present[,idx1,drop=FALSE],1,any)] <- "yes" #if any mutation is present in the certain type of tissue.
    tissue_present[apply(present[,idx1,drop=FALSE],1,all) & apply(freq1,1,max) >= 0.1] <- "all"
    #print( paste("there are",sum(tissue_present=="all"),"all"))
    count2[[paste(patient1,"_",tissue1,"_present",sep="")]] <- tissue_present
  }
  all_snvs[[i]] <- count2[count2$Type != "synonymous SNV",] #non-synonymous variants
  all_snvs2[[i]] <- count2 #all variants
}

allcounts <- NULL
i=1
for (i in 1:length(patients)) {

  patient1 <- patients[i]
  count1 <- all_snvs2[[i]] #all variants from specific patient
  present <- count1[,grep("present",names(count1))] != "no" #creates bool index of present ("not absent")

  sns <- sub("_present","",grep("present",names(count1),value=TRUE))
  tissues <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

  snv_counts <- matrix(0,nrow=2,ncol=length(tissues))
  colnames(snv_counts) <- paste(patient1,tissues,sep="_")
  snv_counts[1,] <- sum(apply(present,1,all)) #number of mutations present in every type of tissue: Polyp and Normal
  snv_counts[2,] <- colSums(present) - sum(apply(present,1,all)) #number of present muts in each type of tissue minus the shared
  allcounts <- cbind(allcounts,snv_counts)
}



##mCRC genes

#Zheng's list of genes
# tmp <- c("APC","KRAS","TP53","SMAD4","PIK3CA","PTEN","ATM","TCF7L2",
# 	 "PIK3R1","ERBB3","AXIN2","FN1","DNAH8","APOB","PTPRT","SYNE1")
#
# tmp2 <- c("APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","NRAS","TCF7L2",
# 	  "FAM123B","SMAD2","CTNNB1","KIAA1804","SOX9","ACVR1B","GPC6","EDNRB")

#Aarons list of shared identical genes
tmpAaron_JP_non_exon = read.delim("JP_non-exonic_sharemutatedgenes_identical_ver4.csv",sep = ",")%>%filter(geneLoc!="intergenic")%>%select(geneName)
tmpAaron_JP_exon = read.delim("JP_exonic_sharemutatedgenes_identical_ver4.csv",sep = ",")%>%filter(geneLoc!="intergenic")%>%select(geneName)
tmpAaron_EP_non_exon = read.delim("EP_non-exonic_sharemutatedgenes_identical_ver4.csv",sep = ",")%>%filter(geneLoc!="intergenic")%>%select(geneName)
tmpAaron_EP_exon = read.delim("EP_exonic_sharemutatedgenes_identical_ver4.csv",sep = ",")%>%filter(geneLoc!="intergenic")%>%select(geneName)
tmp_identical = as.character(unique(rbind(tmpAaron_JP_non_exon,tmpAaron_JP_exon,tmpAaron_EP_non_exon,tmpAaron_EP_exon))[,1])

#get list of Amir's genes (most detrimental mutations). Find Mutations greater than 30 CADD score
tmpAmir = read.delim("SECURE_ EP_JP_WES_Nov2018_hg38_All - results-20190117-183218.csv.tsv",sep = "\t") #prefilterd for CADD score >= 30
tmpAmir = as.character(unique(tmpAmir%>%select(hg38_UCSC_refGene_name2))[-1,1])

#short_list <- tmp
short_list <- unique(c(tmp_identical,tmpAmir))#; rm(tmpAaron_EP_exon,tmpAaron_EP_non_exon,tmpAaron_JP_exon,tmpAaron_JP_non_exon)

#Which genes are known cancer drivers??
#just to remove the last one because its empty
crcgenes <- as.vector(unique(read.csv("PanCanDrivers_Cell2018.csv",header=T,skip = 3)$Gene)); crcgenes = crcgenes[-300]
#crcgenes <- scan("CRCGenes.txt","c")

#add the identical and detrimental genes together into shortlist
short_list <- unique(c(short_list,crcgenes))


#determines which genes were mutated in my samples
mutated_genes <- NULL

i=2
for (i in 1:length(patients))
{
  count1 <- data.frame(all_snvs[[i]]) #only non synonymous variants
  bool = apply(count1[,grep("present",names(count1))],MARGIN = 1,function(x) any(x=="yes"))
  count1 = count1[which(bool),]
    # count1 <- count1[apply(count1[,grep("present",names(count1))]=="yes",1,any),] #switched from "all" to "yes"
  #mutated_genes <- c(mutated_genes,unique(count1$Gene[count1$Gene %in% short_list])) #as.character???
  mutated_genes <- unique(c(mutated_genes,as.character(count1$Gene[count1$Gene %in% short_list]))) #as.character???
}

genes <- short_list[short_list %in% mutated_genes]

nmutation <- matrix(0,nrow=length(genes),ncol=length(patients))
rownames(nmutation) <- genes
colnames(nmutation) <- patients
nmutation2 <- nmutation

i=1  ################
for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  #count1 <- count1[apply(count1[,grep("present",names(count1))]=="yes",1,any),]
  bool = apply(count1[,grep("present",names(count1))],MARGIN = 1,function(x) any(x=="yes"))  #change "all" to "yes". mut present in 1 type of tissue
  count1 = count1[which(bool),]

  #determines number of unique mutations. Can have multiple different mutations for same gene
  table1 <- table(count1$Gene[count1$Gene %in% genes]) # table1[table1>=1]
  #nmutation[names(table1),i] <- pmax(nmutation[names(table1),i],table1) ### remade in the below for/loop
  #gene = "APC" # i dont thin this works right
  for (gene in genes) {
    print(gene)
    if (gene %in% names(table1)) {
      nmutation[gene,i] = table1[names(table1)==gene]
    } else {
      next()
    }
  }

  #determines number of Tissue Types the gene shows up in. One type of tissue gets one, two types of tissues gets 2...
  count2 <- data.frame(Gene=count1$Gene,
		       N_sample=apply(count1[,grep("present",names(count1))] != "no",1,sum)) #mutations which have any value but no for each tissue type
  table2 <- xtabs(count2$N_sample ~ count2$Gene) #table2[table2 > 0]
  table2 <- table2[names(table2) %in% genes]#table2[table2 > 0]
  nmutation2[names(table2),i] <- table2
}

#orders by presence, number of distinct mutations, and number of tissue repeats.
genes <- genes[order(rowSums(nmutation>0),rowSums(nmutation),
		     rowSums(nmutation2),decreasing=TRUE)]
genes <- genes[order(match(genes, c("APC","TP53","KRAS","SYNE1","SMAD4","PIK3CA")))]
nmutation <- nmutation[genes,]
#genes <- c(rownames(nmutation)[apply(nmutation>0,1,sum)>=2],"BRAF","PTEN","PIK3R1")

nmutation <- matrix(0,nrow=length(genes),ncol=length(patients))
rownames(nmutation) <- genes
colnames(nmutation) <- patients
nmutation2 <- nmutation

for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]

  table1 <- table(count1$Gene[count1$Gene %in% genes]) #table1[table1>0]
  # nmutation[names(table1),i] <- pmax(nmutation[names(table1),i],table1)
  # nmutation[,i] = table1[names(table1) %in% genes]
  # gene="APC"
  for (gene in genes) {
    print(gene)
    if (gene %in% names(table1)) {
      nmutation[gene,i] = table1[names(table1)==gene]
    } else {
      next()
    }
  }

  count2 <- data.frame(Gene=count1$Gene,
		       N_sample=apply(count1[,grep("present",names(count1))] != "no",1,sum))
  table2 <- xtabs(count2$N_sample ~ count2$Gene)
  table2 <- table2[names(table2) %in% genes]
  nmutation2[names(table2),i] <- table2
}

genes <- genes[order(rowSums(nmutation>0),rowSums(nmutation),
		     rowSums(nmutation2),decreasing=TRUE)]
genes <- genes[order(match(genes, c("APC","TP53","KRAS","SYNE1","SMAD4","PIK3CA")))]
nmutation <- nmutation[genes,]
max_mutation <- apply(nmutation,1,max)

present_code <- function(present) {
  code <- rep(0,nrow(present))
  ##present in all samples
  code[apply(present,1,all] <- 1
  return(code)
}

mutation_matrix <- matrix(0,nrow=sum(max_mutation),ncol=ncol(allcounts))
rownames(mutation_matrix) <- rep(genes,max_mutation)
colnames(mutation_matrix) <- colnames(allcounts)

#Creats a matrix based on how detrimental different types of mutations are. the values will determine shape/color
for (i in 1:length(patients)) {
  patient1 <- patients[i]
  count1 <- all_snvs[[i]]
  count1 <- count1[count1$Gene %in% genes,]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]

  if (nrow(count1) == 0) next
  sns <- sub("_present","",grep("present",names(count1),value=TRUE))

  present <- count1[,grep("present",names(count1))] != "no"

  j = 1
  for (j in 1:length(genes)) {
    gene1 <- genes[j]
    if (gene1 %in% count1$Gene) {
      idx <- which(count1$Gene == gene1)
      if (length(idx) > 1)
	idx <- idx[order(rowSums(present[idx,]),decreasing=TRUE)]
      for (k in 1:length(idx)) {
	idx1 <- idx[k]
	idx2 <- which(rownames(mutation_matrix) == gene1)[k]
	for (kk in 1:length(sns)) {
	  mutation_matrix[idx2,sns[kk]] <- as.integer(present[idx1,kk])
	}
	if (count1$Type[idx1] %in% c("nonframeshift_indel","frameshift_indel","stopgain_indel")) {
	  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 10
	}
	if (count1$Type[idx1] %in% c("stopgain")) {
	  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 20
	}
	if (count1$Type[idx1] %in% c("splicing")) {
	  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 30
	}
	##if (grepl("indel",count1$Type[idx1])) {
	##  mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] <- mutation_matrix[idx2,grep(patient1,colnames(mutation_matrix))] + 20
	##}
      }
    }
  }
}


gene_table <- rep(0,length(genes))
names(gene_table) <- genes

#determines number of patients affected by gene-specific mutation
for (i in 1:length(patients)) {
  count1 <- all_snvs[[i]]
  count1 <- count1[count1$Gene %in% genes,]
  ##count1 <- count1[apply(count1[,grep("present",names(count1))]=="all",1,any),]
  gene_table[as.character(unique(count1$Gene))] <- gene_table[as.character(unique(count1$Gene))] + 1 #added as.character
}

mutation_matrix2 <- mutation_matrix[nrow(mutation_matrix):1,]

sns2 <- colnames(mutation_matrix2)
pts2 <- sapply(strsplit(sns2,"_",fixed=TRUE),head,1)

#makes an index list of the Patients with an "NA" in between the titles
idx2 <- 1
for (i in 2:length(pts2)) {
  if (pts2[i] != pts2[i-1])
    idx2 <- c(idx2,NA)
  idx2 <- c(idx2,i)
}

mutation_matrix2 <- mutation_matrix2[,idx2]
sns2 <- colnames(mutation_matrix2)

#
#to see what the colors are: plot(1:length(mycolors),1:length(mycolors),col = mycolors,pch = 15)
mutation_point <- function(x,y,code) {
  if (code < 10) cols <- mycolors[3] ##3
  if (code > 10 & code < 20) cols <- mycolors[5] ##6
  if (code > 20 & code < 30) cols <- mycolors[6] ##3
  if (code > 30 & code < 40) cols <- mycolors[9] ##6
  ##if (code > 20) cols <- 5
  code <- code %% 10

  xys <- cbind(c(x-0.4,x-0.4,x+0.4,x+0.4,x-0.4),
	       c(y-0.4,y+0.4,y+0.4,y-0.4,y-0.4))

  if (code == 1) {
    polygon(xys[1:4,],border=NA,col=cols)
  }
}


library(GenomicRanges)
##library(SNPchip)

#########load CNA data from TitanCNA outputs###

alltitancnaresults = NULL
x = list.files(path = "titan",pattern = ".TitanCNA.RData",full.names = T)
for (i in 1:length(x)) {
  print(x[i])
  load(x[i])
  tmp = titancnaresults[1] #the titan results have 2 elements. So i've chosen to use the 1st one here
  alltitancnaresults = c(alltitancnaresults,tmp)
}

#
# load("mCRC_All_TitanCNA.RData")
# tmp1 <- alltitancnaresults
#
# load("mCRC_Leung_TitanCNA_Sequenza.RData")
# tmp2 <- alltitancnaresults
#
# load("CRC_Lim_TitanCNA_choice.RData")
# tmp3 <- alltitancnaresults
#
# load("CRC_Lee_TitanCNA_choice.RData")
# tmp4 <- alltitancnaresults
#
# load("mCRC_Kim2_TitanCNA.RData")
# tmp5 <- alltitancnaresults
#
# alltitancnaresults <- c(tmp1,tmp2,tmp3,tmp4,tmp5)


gc()

source("TitanCNA_04_TitanCNA2seg.R")
# chr_sizes <- read.table("hg38.chrom.sizes",as.is=TRUE)[1:24,]
# chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes[-25,2])))

# idx1 <- which(sapply(strsplit(samples,"_",fixed=TRUE),"[",1) %in% patients)
# idx1 <- which(sapply(strsplit(names(samples),"_",fixed=TRUE),"[",1) %in% patients)
# samples <- samples[idx1]
# samplelabels <- samplelabels[idx1]

samples = str_remove(str_remove(list.files(path = "titan",pattern = ".TitanCNA.RData",full.names = T),pattern = "titan/"),pattern = ".TitanCNA.RData")
#samplelabels = c(rep("Polyp",14),rep("Normal",2),rep("Polyp",11),rep("Normal",2))
samplelabels = samples

names(alltitancnaresults) = samples

##adjust segmean for overall ploidy and purity

all_cna <- vector("list",length=length(samples))
names(all_cna) <- samplelabels
ploidys <- rep(0,length(samples))
names(ploidys) <- samplelabels

i=1
for (i in 1:length(samples)) {
  print(samples[i])
  sn1 <- samples[i]
  cna <- titancna2seg(alltitancnaresults[[sn1]]$results,
		      alltitancnaresults[[sn1]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$allelicratio <- round(cna$allelicratio,3)
  cna$logcopynumberratio <- round(cna$logcopynumberratio,3)

  ploidy1 <- cna$ploidy[1]

  max_prev <- max(cna$cellularprevalence,na.rm=TRUE)
  ploidys[i] <- (ploidy1-2*(1-max_prev))/max_prev

  purity1 <- cna$normalproportion[1]
  ploidy2 <- ploidy1 * (1 - purity1) + 2 * purity1

  prev1 <- cna$cellularprevalence * (1 - purity1)
  purity2 <- max(prev1,na.rm=TRUE)
  prev1[is.na(prev1)] <- purity2
  seg_mean <- cna$seg.mean + log2(ploidy2/2)
  seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))
  ##seg_mean_adj[1+(2^seg_mean-1)/prev1 <= 0] <- -2 ##Homozygous deletions

  cna$seg.mean.adj <- round(seg_mean_adj,3)
  all_cna[[i]] <- cna
}
##save(all_cna,file="mCRC_All_TitanCNA_Adj.RData")

sns <- samplelabels
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
tissues <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

all_gr <- NULL
for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  print(sn1)
  cna <- all_cna[[i]]

  gr <- GRanges(seqnames=cna$chrom,
		ranges=IRanges(cna$loc.start,cna$loc.end),
		strand="*")

  if (is.null(all_gr)) {
    all_gr <- gr
  }
  else {
    all_gr <- disjoin(c(gr,all_gr))
  }
}

all_cna_seg <- data.frame(Chr=seqnames(all_gr),
			  Start=start(all_gr),
			  End=end(all_gr),
			  stringsAsFactors=FALSE)

for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  print(sn1)
  cna <- all_cna[[i]]

  cna$cellularprevalence[cna$LOHcall == "HET"] <-
    max(cna$cellularprevalence,na.rm=TRUE)

  seg.mean <- cna$seg.mean.adj

  cn1 <- rep(NA,length(all_gr))
  gr <- GRanges(seqnames=cna$chrom,
		ranges=IRanges(cna$loc.start,cna$loc.end),
		strand="*")
  idx <- findOverlaps(all_gr,gr)
  cn1[queryHits(idx)] <- seg.mean[subjectHits(idx)]

  all_cna_seg <- cbind(all_cna_seg,cn1)
  names(all_cna_seg)[length(all_cna_seg)] <- sn1
}

##save(all_cna,all_cna_seg,all_gr,file="mCRC_All_TitanCNA_Adj.RData")


seg_list <- read.table("CNA_Cytoband.txt") #specific regions of interest (zheng has 14)
seg_list_gr <- GRanges(seqnames=seg_list$V1,
		       ranges=IRanges(seg_list$V2,seg_list$V3),
		       strand="*")

weight_matrix <- matrix(0,nrow=nrow(seg_list),ncol=nrow(all_cna_seg))
disjoin_gr <- disjoin(c(seg_list_gr,all_gr))
for (i in 1:nrow(seg_list)) {
  idx1 <- which(overlapsAny(disjoin_gr,seg_list_gr[i]))
  idx2 <- findOverlaps(all_gr,disjoin_gr[idx1])
  weight_matrix[i,queryHits(idx2)] <- width(disjoin_gr[idx1[subjectHits(idx2)]])
}
weight_matrix <- weight_matrix/rowSums(weight_matrix)

seg_list_matrix <- data.frame(Chr=seqnames(seg_list_gr),
			      Start=start(seg_list_gr),
			      End=end(seg_list_gr),
			      Label=seg_list$V4,
			      stringsAsFactors=FALSE)

for (i in 1:length(samples)) {
  sn1 <- samplelabels[i]
  seg_cn <- all_cna_seg[[sn1]]
  seg_cn[is.na(seg_cn)] <- 0
  seg_cn <- weight_matrix %*% seg_cn
  seg_list_matrix <- cbind(seg_list_matrix,cn=seg_cn)
  names(seg_list_matrix)[ncol(seg_list_matrix)] <- sn1
}

##One value for each tumor site
seg_list_table <- NULL
patients <- unique(patients)
ploidys2 <- numeric(0)

for (i in 1:length(patients)) {

  patient1 <- patients[i]
  sns <- grep(patient1,samplelabels,value=TRUE)

  cna_list1 <- seg_list_matrix[sns]

  tissues <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)
  tissue_types <- unique(tissues[!is.na(tissues)])

  cna_table1 <- matrix(0,nrow=nrow(seg_list_matrix),ncol=length(sns))
  rownames(cna_table1) <- seg_list_matrix$Label
  colnames(cna_table1) <- sns
  cna_table1[,] <- as.matrix(seg_list_matrix[,sns])

  cna_table2 <-  matrix(0,nrow=nrow(seg_list_matrix),ncol=length(tissue_types))
  rownames(cna_table2) <- seg_list_matrix$Label
  colnames(cna_table2) <- paste(patient1,tissue_types,sep="_")

  for (j in 1:length(tissue_types)) {
    idx1 <- tissues %in% tissue_types[j]
    cna_table2[,j] <- rowMeans(cna_table1[,idx1,drop=FALSE])
    ploidys2 <- c(ploidys2,mean(ploidys[grep(paste(patient1,tissue_types[j],sep="_"),names(ploidys))]))
    names(ploidys2)[length(ploidys2)] <- paste(patient1,tissue_types[j],sep="_")
  }
  seg_list_table <- cbind(seg_list_table,cna_table2)
}

seg_list_table <- seg_list_table[nrow(seg_list_table):1,]

sns2 <- colnames(seg_list_table)
pts2 <- sapply(strsplit(sns2,"_",fixed=TRUE),head,1)

idx2 <- 1
for (i in 2:length(pts2)) {
  if (pts2[i] != pts2[i-1])
    idx2 <- c(idx2,NA)
  idx2 <- c(idx2,i)
}

seg_list_table <- seg_list_table[,idx2]
ploidys2 <- ploidys2[idx2]

sns2 <- colnames(seg_list_table)

#functions for making points
cna_point <- function(x,y,code) {
  if (code < 0) {
    cols <- colorRamp(c(mycolors[c(2,8)]))(pmax(1+code/2,0))
    cols <- rgb(red=cols[1,1],green=cols[1,2],blue=cols[1,3],
		maxColorValue=255)
      ##rgb(red=pmax(1+code/2,0),green=pmax(1+code/2,0),blue=1)
  }
  else if (code >= 0) {
    cols <- colorRamp(c(mycolors[c(1,8)]))(pmax(1-code/2,0))
    cols <- rgb(red=cols[1,1],green=cols[1,2],blue=cols[1,3],
		maxColorValue=255)
      ##rgb(green=pmax(1-code/2,0),blue=pmax(1-code/2,0),red=1)
  }

  xys <- cbind(c(x-0.4,x-0.4,x+0.4,x+0.4,x-0.4),
	       c(y-0.4,y+0.4,y+0.4,y-0.4,y-0.4))

  if (code != 0) {
    polygon(xys[1:4,],border=NA,col=cols)
  }
}

col_labels <- rbind(colorRamp(c(mycolors[c(2,8)]))(c(0.5)),
		    colorRamp(c(mycolors[c(1,8)]))(1-log2((3:8)/2)/2))
col_labels <- rgb(red=col_labels[,1],green=col_labels[,2],
		  blue=col_labels[,3],
		  maxColorValue=255)
names(col_labels) <- c(1,3:8)


load("mCRC_All_TitanCNA_CRCgenes.RData")
LOH_genes <- genes_cn[samplelabels,,rev(c("APC","TP53","SMAD4"))]
LOH_genes2 <- mutation_matrix2[1:3,]
rownames(LOH_genes2) <- dimnames(LOH_genes)[[3]]
for (i in which(!is.na(colnames(LOH_genes2)))) {
  sn1 <- colnames(LOH_genes2)[i]
  idx1 <- grep(sn1,rownames(LOH_genes))
  if (length(idx1) == 1) {
    LOH_genes2[,i] <- (LOH_genes[idx1,2,] == 0) + (LOH_genes[idx1,1,] == 1)
  }
  else {
    LOH_genes2[,i] <- (colSums(LOH_genes[idx1,2,] == 0) >= length(idx1)/2) +
      (colSums(LOH_genes[idx1,1,] == 1) >= length(idx1)/2)
  }
}

loh_point <- function(x,y,code) {
  xys <- cbind(c(x-0.4,x-0.4,x+0.4,x+0.4,x-0.4),
	       c(y-0.4,y+0.4,y+0.4,y-0.4,y-0.4))

  if (code == 1) {
    cols <- "lightgray"
    polygon(xys[1:4,],border=NA,col=cols)
  }
  if (code == 2) {
    cols <- "darkgray"
    polygon(xys[1:4,],border=NA,col=cols)
  }
}


pdf("mCRC_All_Filtered_SNV_Indel_CNA_LOH_Matrix.pdf",width=10,height=8)

layout(matrix(c(3,1,4,5,2,6),nrow=3),widths=c(8.5,1.5),heights=c(1.2,7.3,1.5))
layout.show(n=6)
## Figure A
par(mar=c(0.5,5,0.5,0.5))
plot(1,1,
     ylim=c(1,nrow(mutation_matrix2)+nrow(seg_list_table)+nrow(LOH_genes2)+2),
     xlim=c(1,ncol(mutation_matrix2)),
     xlab="",ylab="",type="n",axes=FALSE)
genes2 <- rownames(mutation_matrix)
genes2[-match(genes,genes2)] <- NA
genes2 <- rev(genes2)
#left side lables
axis(2,at=1:nrow(seg_list_table),labels=rownames(seg_list_table),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0,0,0))
axis(2,at=1:nrow(LOH_genes2)+nrow(seg_list_table)+1,labels=rownames(LOH_genes2),
     las=2,cex.axis=0.9,tick=FALSE,mgp=c(0,0,0))
axis(2,at=1:nrow(mutation_matrix2)+nrow(seg_list_table)+nrow(LOH_genes2)+2,
     labels=genes2,las=2,cex.axis=0.9,tick=FALSE,mgp=c(0,0,0))
##axis(1,at=1:ncol(mutation_matrix2),labels=sns2,
##     las=2,cex.axis=0.7,tick=FALSE)
abline(v=which(is.na(mutation_matrix2[1,])),lty=3)
abline(h=which(rownames(mutation_matrix2)[-1] != rownames(mutation_matrix2)[-nrow(mutation_matrix2)])+0.5+nrow(seg_list_table)+nrow(LOH_genes2)+2,lty=3)
abline(h=0:(nrow(LOH_genes2)+1)+nrow(seg_list_table)+0.5+1,lty=3) #makes dotted horizontal lines
abline(h=1:nrow(seg_list_table)+0.5,lty=3) #makes more dottend horizontal lines

#adds shapes to plot
for (i in 1:ncol(seg_list_table)) {
  if (is.na(seg_list_table[1,i])) next
  for (j in 1:nrow(seg_list_table)) {
    if (seg_list_table[j,i] != 0)
      cna_point(i,j,seg_list_table[j,i])
  }
}

for (i in 1:ncol(LOH_genes2)) {
  if (is.na(LOH_genes2[1,i])) next
  for (j in 1:nrow(LOH_genes2)) {
    loh_point(i,j+nrow(seg_list_table)+1,LOH_genes2[j,i])
  }
}

for (i in 1:ncol(mutation_matrix2)) {
  if (is.na(mutation_matrix2[1,i])) next
  for (j in 1:nrow(mutation_matrix2)) {
    if (mutation_matrix2[j,i] > 0)
      mutation_point(i,j+nrow(seg_list_table)+nrow(LOH_genes2)+2,
		     mutation_matrix2[j,i])
  }
}

par(mar=c(0.5,0,0.5,0.5))
par(mgp=c(-19,-21,-22))
idx3 <- rep(NA,nrow(mutation_matrix))
idx3[match(names(gene_table),rownames(mutation_matrix))] <- 1:length(gene_table)
gene_table2 <- gene_table[idx3]
gene_table2 <- gene_table2[length(gene_table2):1]
gene_table2 <- c(rep(NA,nrow(seg_list_table)+nrow(LOH_genes2)+2),gene_table2)

barplot(gene_table2,width=0.8,space=0.25,axisnames=FALSE,
	ylim=c(0.6,length(gene_table2)-0.4),cex.axis=1,cex.names=1,
	xlab="No of affected\npatients",horiz=TRUE)

#legend("bottomleft",c("Mutation","Missense","Nonsense","","Deletion/LOH","Het deletion","LOH","","Copy number",names(col_labels)),
legend("bottomleft",c("Mutation","Nonsynonymous","Indel","Stopgain","Splicing","","Deletion/LOH","Het deletion","LOH","","Copy number",names(col_labels)),
       fill=c("white",mycolors[c(3,5,6,9)],"white","white","darkgray","lightgray","white","white",col_labels),
       text.font=c(2,1,1,1,1,1,2,1,1,1,2,1,1,1,1,1,1,1),
       border=NA,bty="n",cex=1)


par(mar=c(0,5,0.5,0.5))
par(mgp=c(3,1,0))

allcounts2 <- allcounts[,idx2]

xlim <- c(0.6,ncol(allcounts2)-0.4)
ylim <- c(0,max(colSums(allcounts2)+100,na.rm=TRUE))

barplot(allcounts2,width=0.8,space=0.25,
	legend.text=c("Shared","Not Shared"), col=c("black","lightgrey"),
	names.arg=sns2,xlim=xlim,ylim=ylim,ylab="No of mutations",
	cex.axis=1,cex.names=1,las=2,
	args.legend=list(x="topleft",cex=1,bty="n"),
	axes=TRUE,axisnames=FALSE)


par(mar=c(5,5,0,0.5))
barplot(ploidys2,width=0.8,space=0.25,xlim=xlim,ylim=c(4,0),
	names.arg=sns2,ylab="Ploidy",
	cex.axis=1,cex.names=1,las=2,
	axes=TRUE,axisnames=TRUE)
lines(c(xlim[1]-0.5,xlim[2]+0.5),c(2,2),lty=3,lwd=2)


par(mar=c(0.5,0,0.5,0.5))

plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
##legend(0,1,c("Missense","Nonsense",names(col_labels)),
##       fill=c(mycolors[c(3,5)],col_labels),border=NA,bty="n",cex=1)

dev.off()

