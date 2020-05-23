cut_chr <- function(chr, chr_size, bin_size=1000000){
	num_bin=floor(chr_size/bin_size)
	starts <- c()
	ends <- c()
	for(i in 0:(num_bin-1)){
		starts <- c(starts,bin_size*i+1)
		ends <- c(ends,bin_size*(i+1))
	}
	ends[num_bin]=chr_size
	bin_table <- data.frame(Chr=rep(chr,num_bin),Start=starts, End=ends)

	return(bin_table)
}


polish_cna_seg <- function(cna_table) {
	#for (i in 2:nrow(cna_table)) {
	#	cna_table[i,2] = as.numeric(cna_table[i-1,3]+1)
	#}

	cum_rows <- 0
	for (i in 1:22){
		cna_chr <- cna_table[which(cna_table$Chr==i),]
		max_end <- max(cna_chr$End)
		max_index <- which(cna_table$Chr==i & cna_table$End==max_end)
		cna_table$End[max_index] = chr_sizes[i]

		min_start <- min(cna_chr$Start)
		if (min_start > 1){
			cur_seg <- cna_chr[which(cna_chr$Start==min_start),]
			cur_seg_start <- as.numeric(cur_seg$Start)[1]
			cur_seg_end <- as.numeric(cur_seg$End)[1]
			min_index <- which(cna_table$Start==cur_seg_start & cna_table$End==cur_seg_end)
			cna_table$Start[min_index] = 1
		}

		num_rows <- nrow(cna_chr)
		for (j in (cum_rows+2):(cum_rows+num_rows)){
			cna_table[j,2] = as.numeric(cna_table[j-1,3]+1)
		}
		
		cum_rows = cum_rows+num_rows

	}

	return(cna_table)
}


search_seg <- function(cna_data,chr,end_pos) {
	cna_data2 <- cna_data[which(cna_data$Chr==chr),]
	start_points <- cna_data2$Start
	end_points <- cna_data2$End
	
	for (i in 1:length(end_points)) {
		st <- start_points[i]
		ed <- end_points[i]
		if (end_pos<=as.numeric(ed)) {
			if(i > 1) {
				difference = end_pos-end_points[i-1]
				if (difference<5e5) {
					start_point = as.numeric(start_points[i-1])
					end_point = as.numeric(end_points[i-1])
				}

				else{
					start_point = as.numeric(st)
					end_point = as.numeric(ed)
				}
			}

			else{
				start_point = as.numeric(st)
				end_point = as.numeric(ed)
			}
			break
		}
	}
	return(c(start_point,end_point))
}


cna_seg2bin <- function(cna_table, bin_size=1000000) {
	#cna_table <- polish_cna_seg(cna_table)
	bin_table <- data.frame(Chr=as.numeric(),Start=as.numeric(),End=as.numeric())
	for (i in 1:length(chr_sizes)) {
		bin_table1 <- cut_chr(i,chr_sizes[i])
		bin_table <- merge(bin_table, bin_table1, all=TRUE, sort=FALSE)

	}

	names_cna <- names(cna_table)
	for(i in 4:length(names_cna)) {
		bin_table[,names_cna[i]] = rep(0, nrow(bin_table))
	}

	for(i in 1:nrow(bin_table)){
		cur_chr <- as.numeric(bin_table$Chr[i])
		cur_bin_end <- as.numeric(bin_table$End[i])
		cur_seg <- search_seg(cna_table,cur_chr,cur_bin_end)
		cur_seg_start <- cur_seg[1]
		cur_seg_end <- cur_seg[2]
		cur_cna <- cna_table[which(cna_table$Start==cur_seg_start & cna_table$End==cur_seg_end),]
		bin_table[i,4:ncol(bin_table)] = as.numeric(cur_cna[,4:ncol(cur_cna)])
	}
	
	return(bin_table)

}



chr_sizes <- read.table("broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE)[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

##Plot average CNA from MRS
##Highlight met unique CNAs
##Label cancer genes
##Remove small segments

seg_cutoff <- 1e6

setwd("/Users/huzheng/Desktop/Research/CurtisLab/BrainMetSeeding/DataSummary/mutationLandscape/PANCAN")

mycolors <- c("#D55E00","#0072B2","#009E73","#E69F00","#CC79A7","#F0E442",
              "#999999","#FFFFFF")

#tmp <- c("APC","KRAS","TP53","SMAD4","PIK3CA","PTEN","ATM","TCF7L2",
#	 "PIK3R1","ERBB3","AXIN2","FN1","DNAH8","APOB","PTPRT","SYNE1")

#tmp2 <- c("APC","TP53","KRAS","PIK3CA","FBXW7","SMAD4","NRAS","TCF7L2",
#	  "FAM123B","SMAD2","CTNNB1","KIAA1804","SOX9","ACVR1B","GPC6","EDNRB")

#short_list <- unique(c(tmp,tmp2))

#genelist1 <- as.character(read.table("intogen_drivers_data.txt",header=T)$GeneSymbol)
#genelist2 <- as.character(read.table("cancer_gene_census_grch38_v80.txt",header=T)$GeneSymbol)
#genelist1 <- scan("known_drivers_369genes.txt","c")
#genelist2 <- as.vector(unique(read.csv("PanCanDrivers_Cell2018.csv",header=T)$Gene))
short_list <- c()
genelist2 <- as.vector(read.table("Colon_cancer_driver_genes.txt",header=T)$Gene)

print(genelist2)

#short_list <- unique(c(short_list,genelist1,genelist2))
short_list <- unique(c(short_list,genelist2))
print(length(short_list))

library(GenomicRanges)

samples <- read.table("alltumorlabels_CRC.txt",as.is=TRUE,header=T)


source("TitanCNA_05_cna2genes.R")
source("TitanCNA_04_TitanCNA2seg.R")
chr_sizes <- read.table("broad.Homo_sapiens_assembly19.fa.sizes",as.is=TRUE)[1:22,2]
chr_cumsize <- cumsum(c(0,as.numeric(chr_sizes)))

load("mCRC_All_TitanCNA.RData")
tmp1 <- alltitancnaresults

load("mCRC_Leung_TitanCNA_Sequenza.RData")
tmp2 <- alltitancnaresults

load("mCRC_Lim_TitanCNA_choice.RData")
tmp3 <- alltitancnaresults

load("mCRC_Lee_TitanCNA_choice.RData")
tmp4 <- alltitancnaresults

load("mCRC_Kim2_TitanCNA.RData")
tmp5 <- alltitancnaresults

alltitancnaresults <- c(tmp1,tmp2,tmp3,tmp4,tmp5)


rm(tmp)
gc()
##adjust segmean for overall ploidy and purity

all_cna <- vector("list",length=nrow(samples))
names(all_cna) <- samples$Label

for (i in 1:nrow(samples)) {
  sn1 <- samples$Bam[i]
  cna <- titancna2seg(alltitancnaresults[[sn1]]$results,
		      alltitancnaresults[[sn1]]$convergeParams)
  cna$chrom <- as.integer(as.character(cna$chrom))
  cna$loc.end <- pmin(cna$loc.end,chr_sizes[cna$chrom])
  cna$allelicratio <- round(cna$allelicratio,3)
  cna$logcopynumberratio <- round(cna$logcopynumberratio,3)

  ploidy1 <- cna$ploidy[1]
  normal1 <- cna$normalproportion[1]
  if(normal1>0.6) {normal1=0.6}
  ploidy2 <- ploidy1 * (1 - normal1) + 2 * normal1

  prev1 <- cna$cellularprevalence * (1 - normal1)
  purity2 <- max(prev1,na.rm=TRUE)
  prev1[is.na(prev1)] <- purity2
  seg_mean <- cna$seg.mean + log2(ploidy2/2)
  seg_mean_adj <- log2(pmax(2^-2,1+(2^seg_mean-1)/purity2))

  cna$seg.mean.adj <- round(seg_mean_adj,3)
  all_cna[[i]] <- cna
}

sns <- samples$Label
patients <- sapply(strsplit(sns,"_",fixed=TRUE),"[",1)
sites <- sapply(strsplit(sns,"_",fixed=TRUE),"[",2)

allpatients <- unique(patients)

all_diff_genes <- NULL


########################################CNA frequency#######################################
cna_freq <- data.frame(Chr=as.numeric(),Start=as.numeric(),End=as.numeric())
for (i in 1:length(chr_sizes)) {
	cna_freq1 <- cut_chr(i,chr_sizes[i])
	cna_freq <- merge(cna_freq, cna_freq1, all=TRUE, sort=FALSE)

}
cna_freq <- data.frame(cna_freq, P_amp_freq=rep(0,nrow(cna_freq)),P_del_freq=rep(0,nrow(cna_freq)), M_amp_freq=rep(0,nrow(cna_freq)), M_del_freq=rep(0,nrow(cna_freq)),ACN_diff=rep(0,nrow(cna_freq)),RCN_diff=rep(0,nrow(cna_freq))) 

num_primary <- 0
num_met <- 0

cna_list <- vector(mode="list", length=nrow(cna_freq))
#for(i in 1:length(cna_list)){
#	cna_list[[i]] <- c()
#}

for (pt1 in allpatients) {
  print(pt1)
  all_gr <- NULL
  sns1 <- sns[patients == pt1]
  for (sn1 in sns1) {
    cna <- all_cna[[sn1]]
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
  all_gr <- all_gr[width(all_gr)>=seg_cutoff]
  
  cna_seg <- data.frame(Chr=seqnames(all_gr),Start=start(all_gr),
			End=end(all_gr))
  sites1 <- sites[patients == pt1]
  for (site1 in unique(sites1)) {
    n_site1 <- table(sites1)[site1]
    site_seg <- rep(0,nrow(cna_seg))
    for (sn1 in sns1[sites1 == site1]) {
      cna <- all_cna[[sn1]]
      gr <- GRanges(seqnames=cna$chrom,
		    ranges=IRanges(cna$loc.start,cna$loc.end),
		    strand="*")
      idx <- findOverlaps(all_gr,gr)
      site_seg[queryHits(idx)] <- site_seg[queryHits(idx)] + cna$seg.mean.adj[subjectHits(idx)]/n_site1
    }
    cna_seg <- cbind(cna_seg,site_seg,2^site_seg*2,
		     2^site_seg*2-median(2^site_seg*2))
    names(cna_seg)[ncol(cna_seg)-(2:0)] <- paste0(site1,c("_lr","_cn","_rcn"))
  }

  cna_seg <- polish_cna_seg(cna_seg)
  cna_bin <- cna_seg2bin(cna_seg)

  tissues_rcn <- grep("_rcn",names(cna_seg),value=TRUE)
  primary_rcn <- c("P_rcn")
  met_rcn <- setdiff(tissues_rcn,primary_rcn)
  met_rcn <- met_rcn[which(met_rcn != "LN_rcn")]

  tissues_cn <- grep("_cn",names(cna_seg),value=TRUE)
  primary_cn <- c("P_cn")
  met_cn <- setdiff(tissues_cn,primary_cn)
  met_cn <- met_cn[which(met_cn != "LN_cn")]

  if(length(met_rcn>0)) {

  	num_primary <- num_primary+1
  	num_met <- num_met+length(met_rcn)

  	for (k in 1:nrow(cna_bin)){
  		if (cna_bin[k,primary_rcn[1]] > 0.8) { cna_freq[k,"P_amp_freq"]=cna_freq[k,"P_amp_freq"]+1 }
  		#if (cna_bin[k,primary[1]] > 0.8) { cna_freq[k,"P_amp_freq"]=cna_freq[k,"P_amp_freq"]+cna_bin[k,primary[1]] }
  		if (cna_bin[k,primary_rcn[1]] < -0.8) { cna_freq[k,"P_del_freq"]=cna_freq[k,"P_del_freq"]+1 }
  		#if (cna_bin[k,primary[1]] < -0.8) { cna_freq[k,"P_del_freq"]=cna_freq[k,"P_del_freq"]+cna_bin[k,primary[1]] }
  		for (m in met_rcn){
  			if (cna_bin[k,m] > 0.8) { cna_freq[k,"M_amp_freq"]=cna_freq[k,"M_amp_freq"]+1 }
  			#if (cna_bin[k,m] > 0.8) { cna_freq[k,"M_amp_freq"]=cna_freq[k,"M_amp_freq"]+cna_bin[k,m] }
  			if (cna_bin[k,m] < -0.8) { cna_freq[k,"M_del_freq"]=cna_freq[k,"M_del_freq"]+1 }
  			#if (cna_bin[k,m] < -0.8) { cna_freq[k,"M_del_freq"]=cna_freq[k,"M_del_freq"]+cna_bin[k,m] }

  			cna_freq[k,"RCN_diff"]=cna_freq[k,"RCN_diff"]+cna_bin[k,m]-cna_bin[k,primary_rcn[1]]
  			cna_list[[k]] = c(cna_list[[k]],cna_bin[k,m]-cna_bin[k,primary_rcn[1]])
  		}
  	}

  	for (k in 1:nrow(cna_bin)){
  		for (m in met_cn){
  			cna_freq[k,"ACN_diff"]=cna_freq[k,"ACN_diff"]+cna_bin[k,m]-cna_bin[k,primary_cn[1]]			
  		}
  	}

  }
}

cna_freq$P_amp_freq <- cna_freq$P_amp_freq/num_primary
cna_freq$P_del_freq <- cna_freq$P_del_freq/num_primary
cna_freq$M_amp_freq <- cna_freq$M_amp_freq/num_met
cna_freq$M_del_freq <- cna_freq$M_del_freq/num_met

cna_freq$ACN_diff <- cna_freq$ACN_diff/num_met
cna_freq$RCN_diff <- cna_freq$RCN_diff/num_met

median_cn_diff <- c()
for(i in 1:length(cna_list)){
	mcn <- median(cna_list[[i]])
	median_cn_diff <- c(median_cn_diff,mcn)
}

cna_freq <- data.frame(cna_freq, Median_cn_diff=median_cn_diff)

#cna_freq$P_del_freq <- -cna_freq$P_del_freq
#cna_freq$M_del_freq <- -cna_freq$M_del_freq

cum_location <- c()
for(i in 1:nrow(cna_freq)){
	chr <- as.numeric(cna_freq$Chr[i])
	position <- as.numeric(cna_freq$Start[i])
	cum_site <- chr_cumsize[chr]+position
	cum_location <- c(cum_location,cum_site)

}

cna_freq2 <- data.frame(cna_freq, CumLoc=cum_location)
write.table(cna_freq2, "CNA_frequency_CRC.txt", sep="\t", row.names=FALSE, quote=FALSE)

#####################
pdf("CNA_frequency_CRC.pdf",width=10,height=4)
par(mfrow = c(2, 1))
par(mar=c(2,6,2,0))

plot(0,0,ylim=c(-1,1.2),xlim=c(0,sum(as.numeric(chr_sizes))),
	 type="n",xaxt="n",bty="n",ylab="SCNA frequency",xlab="",las=2, main="Primary tumor")
abline(v=chr_cumsize[1:22],lty=2,lwd=1,col="gray60")
text(x=chr_cumsize[1:22] + chr_sizes/2, y=rep(1.1,22),labels=1:22,cex=0.7)

for (i in 1:nrow(cna_freq2)){
	rect(cna_freq2$CumLoc[i],0,cna_freq2$CumLoc[i]+1000000,cna_freq2$P_amp_freq[i],col="#dd3497",border=NA)
	rect(cna_freq2$CumLoc[i],0,cna_freq2$CumLoc[i]+1000000,-cna_freq2$P_del_freq[i],col="#41ab5d",border=NA)
}

plot(0,0,ylim=c(-1,1.2),xlim=c(0,sum(as.numeric(chr_sizes))),
	 type="n",xaxt="n",bty="n",ylab="SCNA frequency",xlab="",las=2, main="Metastasis")
abline(v=chr_cumsize[1:22],lty=2,lwd=1,col="gray60")
#text(x=chr_cumsize[1:22] + chr_sizes/2, y=rep(2.0,22),labels=1:22,cex=0.8)

for (i in 1:nrow(cna_freq2)){
	rect(cna_freq2$CumLoc[i],0,cna_freq2$CumLoc[i]+1000000,cna_freq2$M_amp_freq[i],col="#dd3497",border=NA)
	rect(cna_freq2$CumLoc[i],0,cna_freq2$CumLoc[i]+1000000,-cna_freq2$M_del_freq[i],col="#41ab5d",border=NA)
}

dev.off()

