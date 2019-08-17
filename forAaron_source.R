getSampMutMulti <- function(samples, normal, d, cmedianTh=2.2, original=4.3) {
    numsamples = length(samples)
    rindexTF = vector()
    cindex = vector()
    for (si in 1:length(samples)) {
        if (length(rindexTF) > 0) {
            rindexTF = rindexTF | grepl(paste(samples[si],"\\[",sep=""), d$somatic)
        } else {
            rindexTF = grepl(paste(samples[si],"\\[",sep=""), d$somatic)
        }
        cindex = append(cindex, match(paste(samples[si], "maf", sep=""), colnames(d)))
    }
    rindex = which(rindexTF)
    cindex = as.vector(sapply(cindex, function(x){c(x-1,x,x+1)}))
    
    ncindex = match(paste(normal, "maf", sep=""), colnames(d))
    ncindex = as.vector(sapply(ncindex, function(x){c(x,x+1)}))  #normal samples maf and depth indexes
    
    dronindex = match("dron", colnames(d))
    gnindex = match("geneName", colnames(d))
    glindex = match("geneLoc", colnames(d))
    gfindex = match("functionalClass", colnames(d))
    caddindex = match("CADD_phred", colnames(d))
    gerpindex = match("GERP_RS", colnames(d))
    siftindex = match("SIFT_score", colnames(d))
    polyphenindex = match("Polyphen2_HVAR_pred", colnames(d))
    somindex = match("somatic", colnames(d))
    repindex = match("rep",colnames(d))

    if (! is.na(caddindex)) {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,caddindex,gerpindex,siftindex,polyphenindex,dronindex,repindex,somindex,ncindex)]
    } else {
        res = d[rindex,c(1,2,3,4,5,cindex,gnindex,glindex,gfindex,dronindex,repindex,somindex,ncindex)]
    }
    
    resColnames = colnames(res)
    resAdd = t(data.frame(apply(res, 1, function(x, cindex, original, resColnames) {     #x is every row
              maxMaf = 0
              maxTlod = 0
              ssb = 0                    #combined strand bias check
              ssbc = 0
              ssbp = 0
              refdav = 0
              altdav = 0
              for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample get maxMaf
                  ss = strsplit(as.character(x[j+1]), "\\|")
                  mafTmp = as.numeric(ss[[1]][1])
                  oriTlodTmp = 0
                  if ( grepl("\\|", as.character(x[j])) ) {
                      ori = strsplit(as.character(x[j]), "\\|")
                      oriLod = strsplit(ori[[1]][3], ",")            
                      oriTlodTmp = as.numeric(oriLod[[1]][1])
                  }
                  if (mafTmp > maxMaf) {
                      maxMaf = mafTmp    #define max maf
                  }
                  if (original > 0 & oriTlodTmp > maxTlod) {
                      maxTlod = oriTlodTmp
                  }
              }
              resVector = vector()
              for (j in seq(6,(3+length(cindex)), by=3)) {   #foreach sample
                  sampleName = resColnames[j]
                  sampleNameMaf = paste(sampleName, "mafc", sep="")
                  sampleNameMafa = paste(sampleName, "mafa", sep="")
                  sampleNameCcf = paste(sampleName, "ccf", sep="")
                  sampleNameCcfSD = paste(sampleName, "ccfSD", sep="")
                  sampleNameTime = paste(sampleName, "time", sep="")
                  sampleNameAlt = paste(sampleName, "altc", sep="")
                  sampleNameRef = paste(sampleName, "refc", sep="")
                  oriTlod = 0
                  if ( grepl("\\|", as.character(x[j])) ) {
                      ori = strsplit(as.character(x[j]), "\\|")
                      oriLod = strsplit(ori[[1]][3], ",")
                      oriTlod = as.numeric(oriLod[[1]][1])
                  }
                  ss = strsplit(as.character(x[j+1]), "\\|")
                  mafTmp = as.numeric(ss[[1]][1])
                  cmeme = strsplit(ss[[1]][3], ",")
                  endBias = as.numeric(ss[[1]][2])
                  strandBiases = strsplit(ss[[1]][4], ",")
                  strandBias = as.numeric(strandBiases[[1]][1])
                  strandBiasRef = 0
                  strandBiasFisherP = -1
                  if (length(strandBiases[[1]]) > 1){
                      strandBiasRef = as.numeric(strandBiases[[1]][2])
                      strandBiasFisherP = as.numeric(strandBiases[[1]][3])
                  }
                  mappingBias = as.numeric(ss[[1]][5])
                  cmedianav = as.numeric(cmeme[[1]][2])
                  cmemeSum = sum(as.numeric(cmeme[[1]]))
                  
                  #decide a b allele count
                  refnow = round(as.numeric(x[j+2])*(1-mafTmp))
                  altnow = round(as.numeric(x[j+2])*mafTmp)

                  if ( mafTmp > 0 ) {
                      ssb = ssb + strandBias*altnow
                      ssbp = ssbp + strandBiasFisherP*altnow
                      refdav = refdav + refnow
                      altdav = altdav + altnow
                      ssbc = ssbc + 1
                  }
                  
                  #decide mafNow
                  if ( mafTmp == 0 ) {
                      mafNow = mafTmp
                  } else if (endBias < 0.9 & ((strandBias != 0 & strandBias != 1) | (strandBiasFisherP > 0.7 & refnow >= 10 & altnow >= 5 & mafTmp >= 0.1)) & mappingBias < 0.8 & cmemeSum < 5.2 & cmedianav < cmedianTh) {  
                      mafNow = mafTmp
                      if (original > 0) {   #if oriLod
                          if (oriTlod < original) {
                              if (!(maxMaf > 0.1 & (original == 0 | (original > 0 & maxTlod >= original)))) {   #not called
                                  mafNow = 0
                              }
                          }
                      }
                  } else {
                      if (maxMaf > 0.2 & (original == 0 | (original > 0 & maxTlod >= original)) & mafTmp >= 0.01) {
                          if (cmedianav < 4) {
                              mafNow = mafTmp
                          } else {
                              mafNow = 0
                          }
                      } else if (maxMaf > 0.05 & (original == 0 | (original > 0 & maxTlod < original)) & mafTmp >= 0.04) {
                          if (cmedianav == 1 & mappingBias == 0) {
                              mafNow = mafTmp
                          } else {
                              mafNow = 0
                          }
                      } else {
                          mafNow = 0
                      }
                  }
                  resVector = c(resVector, c(mafNow, 0, 0, 0, 0, refnow, altnow))
                  names(resVector)[(length(resVector)-6):length(resVector)] = c(sampleNameMaf,sampleNameMafa,sampleNameCcf,sampleNameCcfSD,sampleNameTime,sampleNameRef,sampleNameAlt)   #time added
              } #for each sample
              
              ssb = ssb/altdav
              ssbp = ssbp/altdav
              #message(paste(x[1],x[2], sep="  "))
              #message(paste(ssb, ssbp, sep="  "))
              altdav = altdav/ssbc
              refdav = refdav/ssbc
              #if (((ssb >= 0.95 | ssb <= 0.05) & ssbp < 0.1) | ssbp < 0.001) {                         #multiple sample strand bias
              if (((ssb > 0.88 | ssb < 0.12) & numsamples > 1) | ((ssb > 0.9 | ssb < 0.1) & numsamples == 1)) {                         #multiple sample strand bias
              #if ((ssb >= 0.95 | ssb <= 0.05)| ssbp < 0.001) {                                        #multiple sample strand bias 
                  for (j in seq(6,(3+length(cindex)), by=3)) {
                      sampleName = resColnames[j]
                      sampleNameMaf = paste(sampleName, "mafc", sep="")
                      resVector[sampleNameMaf] = 0
                  }
              }

              totalMaf = sum(as.numeric(resVector[grepl("mafc", names(resVector))]))
              totalAlt = sum(as.numeric(resVector[grepl("altc", names(resVector))]))
              totalRef = sum(as.numeric(resVector[grepl("refc", names(resVector))]))
              mergeMAFC = round(totalAlt/(totalAlt+totalRef), 4)
              mergeMAFA = mergeMAFC
              mergedCCF = mergeMAFC
              mergedCCFsd = mergeMAFC
              resVector = c(resVector, totalMaf, totalAlt, totalRef, mergeMAFC, mergeMAFA, mergedCCF, mergedCCFsd)
              names(resVector)[(length(resVector)-6):length(resVector)] = c("totalMaf", "totalAlt", "totalRef", "mergeMAFC", "mergeMAFA", "mergeCCF", "mergeCCFsd")

              resVector                #return the result

          }, cindex=cindex, original=original, resColnames=resColnames)))
    

    res = cbind(res, resAdd)
    res = res[which(res$totalMaf != 0),]
    return(res)
}


#adjust CCF titan for multi samples

adjust.ccf.titan.multi <- function(sampAB, samples, t, titanPath="./titan/", correctColname=FALSE, overadj=1.6, sigTh=0.9) {
    
    numsamples = length(samples)
    purities = vector()
    
    for (i in 1:length(samples)) {
        sn = samples[i]
        message(sn)
        cnv.inputA = paste(titanPath,sn, "_nclones","1",".TitanCNA.segments.txt",sep="")
        if(file_test("-f",cnv.inputA)){
          "NA"
          } else {cnv.inputA = paste(titanPath,sn, "_nclones","2",".TitanCNA.segments.txt",sep = "")}
        #message(cnv.inputA)
        cnvA = read.delim(cnv.inputA)
        cnvA = cnvA[which(!is.na(cnvA$cellularprevalence)),]                #skip NA
        cnvA = cnvA[which(cnvA$num.mark > 9),]                              #skip two few marks
        cnvA$nt = cnvA$copynumber
        if ("minor_cn" %in% colnames(cnvA)) {
            cnvA$nb = cnvA$minor_cn
        } else {
            cnvA$nb = partialRound(cnvA$copynumber*(
                1-((cnvA$allelicratio - cnvA$normalproportion*0.5)/(1-cnvA$normalproportion) - (1-cnvA$cellularprevalence)*0.5)
                /cnvA$cellularprevalence))
        }
        #cnvA$cellularprevalence = sapply(cnvA$cellularprevalence, function(x){if (x == 1){0.99} else {x}})
        
        cnvSeqNames = cnvA$chrom
        if (!grepl("chr",cnvA$chrom[1])){
            cnvSeqNames = paste("chr",cnvA$chrom,sep="")
        }
        snvSeqNames = sampAB$chr
        if (!grepl("chr",sampAB$chr[1])){
            snvSeqNames = paste("chr",sampAB$chr,sep="")
        }
        cnvRangeA = GRanges(seqnames = cnvSeqNames, ranges = IRanges(cnvA$loc.start, end=cnvA$loc.end), strand=rep('+',dim(cnvA)[1]))
        snvRange = GRanges(seqnames = snvSeqNames, ranges = IRanges(sampAB$pos, end=sampAB$pos), strand=rep('+',dim(sampAB)[1]))
        foA = findOverlaps(snvRange, cnvRangeA)
        
        #pa,pu,nt,nb,seg
        message("info table building")
        queHits = queryHits(foA)
        subHits = subjectHits(foA)
        infos = t(sapply(1:dim(sampAB)[1], function(x, queHits, subHits) {                           
                             if (x %in% queHits){
                                 queMIndex = match(x, queHits)
                                 subMIndex = subHits[queMIndex]
                                 c(cnvA$cellularprevalence[subMIndex],     #pa
                                   cnvA$nt[subMIndex],                     #nt
                                   cnvA$nb[subMIndex],                     #nb
                                   subMIndex)                              #seg
                             } else {
                                 c(0,2,1,0)
                             }}, queHits = queHits, subHits = subHits))
        infos = data.frame(infos)
        message("info table built")

        pa1 = infos[,1]
        nt1 = infos[,2]
        nb1 = infos[,3]
        seg1 = infos[,4]
        con1 = as.numeric(cnvA$normalproportion[1])          #pu1 = 1-cnvA$normalproportion[1]
        pu1 = 1 - (con1 + (1-max(as.numeric(pa1)))*(1-con1))
        sAGP = pa1*(1-con1)
        message(pu1)
        purities = append(purities, pu1)
        names(purities)[length(purities)] = sn
        
        sampAB = data.frame(sampAB, pu=pu1, pa=pa1, sAGP=sAGP, nt=nt1, nb=nb1, seg=seg1)
        colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]] = paste(sn,colnames(sampAB)[(dim(sampAB)[2]-5):dim(sampAB)[2]], sep="")
        if ( correctColname == TRUE ) {
            colnames(sampAB) = gsub("\\.","-",colnames(sampAB))
        }
    }
    
    for( i in 1:dim(sampAB)[1]) {  # rescale the maf and calculate CCF

        if (i %% 1000 == 0) {
            message(i)
        }
        
        foundSites = 0             # count how many sites found
        depthTotal = 0
        mafaTotal = 0
        ccfTotal = 0
        ccfsdTotal = 0
        for (j in 1:length(samples)) {
            sn = samples[j]
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            if (maf1 > t) {
                foundSites = foundSites+1
            }
        }
        
        for (j in 1:length(samples)) {
            sn = samples[j]

            pa1 = as.numeric(sampAB[i, match(paste(sn, "pa", sep=""), colnames(sampAB))])
            nt1 = as.numeric(sampAB[i, match(paste(sn, "nt", sep=""), colnames(sampAB))])
            nb1 = as.numeric(sampAB[i, match(paste(sn, "nb", sep=""), colnames(sampAB))])
            maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
            refc1 = as.numeric(sampAB[i, match(paste(sn, "refc", sep=""), colnames(sampAB))])
            altc1 = as.numeric(sampAB[i, match(paste(sn, "altc", sep=""), colnames(sampAB))])
            pu1 = as.numeric(sampAB[i, match(paste(sn, "pu", sep=""), colnames(sampAB))])       #cell purity
            #if (nt1 > 0) pu1 = nt1*pu1/(nt1*pu1+2*(1-pu1))                                     #effective purity
            sAGP = as.numeric(sampAB[i, match(paste(sn, "sAGP", sep=""), colnames(sampAB))])    #segmental aneu- ploid genome proportion
            if (nt1 > 0 & nt1 != 2) pu1 = nt1*sAGP/(nt1*sAGP+2*(1-sAGP))                                   #effective purity
            
            if (maf1 > 0) {
                if ((maf1 > t & foundSites >= 2) | numsamples == 1) {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"unknown",overadj=overadj,sigTh=sigTh)
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                    sampAB[i, match(paste(sn, "time", sep=""), colnames(sampAB))] = CCF1[2]   #need to change
                } else {
                    CCF1 = computeCCF(maf1,refc1,altc1,pu1,pa1,sAGP,nt1,nb1,"late",overadj=overadj,sigTh=sigTh)
                    sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
                    sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
                    sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
                    sampAB[i, match(paste(sn, "time", sep=""), colnames(sampAB))] = CCF1[2]
                }
            }

            #for merged MAF
            dep1 = as.numeric(sampAB[i, match(paste(sn, "d", sep=""), colnames(sampAB))])
            mafa1 = as.numeric(sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))])
            ccf1 = as.numeric(sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))])
            ccfsd1 = as.numeric(sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))])
            depthTotal = depthTotal + dep1
            mafaTotal = mafaTotal + dep1*mafa1
            ccfTotal = ccfTotal + dep1*ccf1
            ccfsdTotal = ccfsdTotal + dep1*ccfsd1
            #for merged MAF
            
        }     #for each sample
        sampAB[i, match("mergeMAFA", colnames(sampAB))] = round(mafaTotal/depthTotal, 5)
        sampAB[i, match("mergeCCF", colnames(sampAB))] = round(ccfTotal/depthTotal, 5)
        sampAB[i, match("mergeCCFsd", colnames(sampAB))] = round(ccfsdTotal/depthTotal, 5)
    }
    return(sampAB)
}


computeCCF <- function(f, A, S, pu, pa, sAGP, nt, nb, prior="unknown", overadj=1.6, sigTh=0.90) {

    ccf = 0
    ccf2 = 0
    sd = 0
    cc <- seq(0.02, 1, by = 0.01)
    evoType = "A1/A2/B/C"
    N = A + S
    nc = nt * pa + 2 * (1 - pa)
    nc2 = nt * sAGP + 2 * (1 - sAGP)

    #message(nc)
    #message(paste(c(f,A,S,pu,pa,sAGP,nt,nb),collapse=" "))
    if (nb == 1 & nt == 2) {   #normal diploid
        ccf = 2*(f/pu)
        ff = pu*cc/2
        Ms = computeSD(N, S, ff)
        ccf2 <- Ms$M1
        sd <- Ms$SD
    }
    else if (nt == 1) {        #het deletion
        ccf = (f/pu)*nc
        ff.C <- pu*cc/nc                       #dbinom
        Ms.C <- computeSD(N, S, ff.C)          #dbinom
        ccf2 <- Ms.C$M1                        #dbinom
        sd <- Ms.C$SD
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cp.A <- pEarly.a/Ptot
        cp.CD <- 1 - cp.A
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        cp.AC <- 1 - cp.D
        cp.AD <- 1 - cp.C
        if (cp.A >= sigTh) {
            evoType <- "A1"
        } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh) {
            evoType <- "B/C"
        } else if (cp.C >= sigTh) {
            evoType <- "B"
        } else if (cp.D >= sigTh) {
            evoType <- "C"
        } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh) {
            evoType <- "A1/B"
        } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh) {
            evoType <- "A1/C"
        }
    } else if (nb == 0 & nt == 0) {        #homozygous deletion
        evoType <- "C"
        ccf = (f*nc2)/sAGP
        ff.C = cc[cc<=(1-sAGP)]/nc2
        Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=(1-sAGP)])
        if (sAGP > 0.98) {
            ff.C = cc[cc<=0.02]/nc2
            Ms.C = computeSD(N, S, ff.C, cc=cc[cc<=0.02])
        }
        ccf2 <- (Ms.C$M1)/pu                    #dbinom
        sd <- Ms.C$SD
        if (is.nan(sd)) {
            sd = 0.01
        }
        #message(paste(c(ccf,ccf2,sd,f),collapse=" "))
    } else if (nb == 0 | nt == 2 * nb) {   #NLOH or other balanced CNAs
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pLate + pEuploid
        cpEarly.a <- pEarly.a/Ptot
        cpLate.eup <- 1 - cpEarly.a
        cpLate <- pLate/Ptot
        cpEup <- pEuploid/Ptot
        if (Ptot > 0) {
            if (cpEarly.a >= sigTh){
                evoType <- "A1"
            } else if (cpLate.eup >= sigTh){
                evoType <- "B/C"
            } else if (cpLate >= sigTh){
                evoType <- "B"
            } else if (cpEup >= sigTh){
                evoType <- "C"
            }
        }
        allprobs = c(pEarly.a, pLate, pEuploid)
        names(allprobs) = c("pEarly.a", "pLate", "pEuploid")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
            ccf = (f/pu)*nc - (nt - nb - 1)*pa
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            sd <- Ms.A$SD
        } else {
            ccf = (f/pu)*nc
            ff.C <- pu*cc/nc                        #dbinom
            Ms.C <- computeSD(N, S, ff.C)           #dbinom
            ccf2 <- Ms.C$M1                         #dbinom
            sd <- Ms.C$SD
        }
    } else if (nb >= 1 & nt > 2) {
        fh.ea <- (sAGP * (nt - nb) + 1 - sAGP)/nc2
        fl.ea <- (sAGP * (nt - nb))/nc2
        fh.eb <- (nb * sAGP + 1 - sAGP)/nc2
        fl.eb <- nb * sAGP/nc2
        fh.t <- sAGP/nc2
        fh.e <- (1 - sAGP)/nc2
        pEarly.a <- pbeta(fh.ea, S+1, A+1) - pbeta(fl.ea, S+1, A+1)
        pEarly.b <- pbeta(fh.eb, S+1, A+1) - pbeta(fl.eb, S+1, A+1)
        pLate <- pbeta(fh.t, S+1, A+1)
        pEuploid <- pbeta(fh.e, S+1, A+1)
        Ptot <- pEarly.a + pEarly.b + pLate + pEuploid
        cp.A <- pEarly.a/Ptot
        cp.B <- pEarly.b/Ptot
        cp.C <- pLate/Ptot
        cp.D <- pEuploid/Ptot
        #message(paste(c(cp.A, cp.B, cp.C, cp.D), collapse=" "))
        cp.AB <- 1 - cp.C - cp.D
        cp.AC <- 1 - cp.B - cp.D
        cp.AD <- 1 - cp.B - cp.D
        cp.BC <- 1 - cp.A - cp.D
        cp.BD <- 1 - cp.A - cp.C
        cp.CD <- 1 - cp.A - cp.B
        cp.ABC <- 1 - cp.D
        cp.ABD <- 1 - cp.C
        cp.ACD <- 1 - cp.B
        cp.BCD <- 1 - cp.A
        if (Ptot > 0) {
            if (cp.A >= sigTh) {                   # earl A
                evoType = "A1"
            } else if (cp.B >= sigTh){
                evoType <- "A2"
            } else if (cp.C >= sigTh){
                evoType <- "B"
            } else if (cp.D >= sigTh){
                evoType <- "C"
            } else if (cp.CD >= sigTh & cp.C < sigTh & cp.D < sigTh){
                evoType <- "B/C"
            } else if (cp.AB >= sigTh & cp.A < sigTh & cp.B < sigTh){
                evoType <- "A1/A2"
            } else if (cp.AC >= sigTh & cp.A < sigTh & cp.C < sigTh){
                evoType <- "A1/B"
            } else if (cp.AD >= sigTh & cp.A < sigTh & cp.D < sigTh){
                evoType <- "A1/C"
            } else if (cp.BC >= sigTh & cp.B < sigTh & cp.C < sigTh){
                evoType <- "A2/B"
            } else if (cp.BD >= sigTh & cp.B < sigTh & cp.D < sigTh){
                evoType <- "A2/C"
            } else if (cp.BCD >= sigTh & cp.BC < sigTh & cp.BD < sigTh & cp.CD < sigTh & cp.B < sigTh & cp.C < sigTh & cp.D < sigTh){
                evoType <- "A2/B/C"
            } else if (cp.ABC >= sigTh & cp.BC < sigTh & cp.AB < sigTh & cp.AC < sigTh & cp.B < sigTh & cp.C < sigTh & cp.A < sigTh){
                evoType <- "A1/A2/B"
            } else if (cp.ABD >= sigTh & cp.AB < sigTh & cp.AD < sigTh & cp.BD < sigTh & cp.B < sigTh & cp.D < sigTh & cp.A < sigTh){
                evoType <- "A1/A2/C"
            } else if (cp.ACD >= sigTh & cp.AC < sigTh & cp.AD < sigTh & cp.CD < sigTh & cp.A < sigTh & cp.D < sigTh & cp.C < sigTh){
                evoType <- "A1/B/C"
            }
        }
        allprobs = c(pEarly.a, pEarly.b, pLate, pEuploid)
        names(allprobs) = c("pEarly.a", "pEarly.b", "pLate", "pEuploid")
        maxType = names(allprobs[match(max(allprobs),allprobs)])
        #if (maxType == "pEarly.a" & prior != "late") {
        if (evoType == "A1" & prior != "late") {
            ccf = (f/pu)*nc - (nt - nb - 1)*pa          #early A1
            ff.A <- pu*(cc - pa + (nt - nb) * pa)/nc    #dbinom
            #ff.A <- (cc - sAGP + (nt - nb) * sAGP)/nc2
            Ms.A <- computeSD(N, S, ff.A)               #dbinom
            ccf2 <- Ms.A$M1                             #dbinom
            #ccf2 = ccf2/pu
            sd <- Ms.A$SD
        #} else if (maxType == "pEarly.b" & prior != "late") {
        } else if (evoType == "A2" & prior != "late") {    
            ccf = (f/pu)*nc - (nb - 1)* pa         #early A2
            ff.B <- pu*(cc - pa + nb * pa)/nc      #dbinom
            Ms.B <- computeSD(N, S, ff.B)          #dbinom
            ccf2 <- Ms.B$M1                        #dbinom
            sd <- Ms.B$SD
        } else {
            ccf = (f/pu)*nc                        #other
            ff.C <- pu*cc/nc                       #dbinom
            Ms.C <- computeSD(N, S, ff.C)          #dbinom
            ccf2 <- Ms.C$M1                        #dbinom
            sd <- Ms.C$SD
        }
    }
    if ( f > 0.1 & ccf >= overadj ) {    #correct for over-adjustment
        if (evoType != "A1") {
            if ( (nt-nb) >= 3 ) {
                ccf = (f/pu)*2
            } else if ( nt >= 2 & (nt-nb) < 3 ) {
                ccf = (f/pu)*nc - (nt - nb - 1)*pa
            } else {
                ccf = (f/pu)*nc
            }
        } else {
            ccf = f*2
        }
    }
    if ( f > 0.7 & ccf2 < 0.05 & is.nan(sd) ) {  #probably un-identified LOH
        ccf = f*2
        ccf2 = 1
        sd = 0.01
    }
    return(c(ccf, evoType, ccf2, sd))
}

computeSD <- function(N, S, f, cc=seq(0.02, 1, by = 0.01)) {
    M1list <- c()
    M2list <- c()
    MLElist <- c()
    for (ii in 1:length(N)) {
        PF <- sum(dbinom(S[ii], N[ii], f), na.rm = TRUE)
        M1 <- sum(dbinom(S[ii], N[ii], f) * cc, na.rm = TRUE)/PF
        M2 <- sum(dbinom(S[ii], N[ii], f) * cc^2, na.rm = TRUE)/PF
        M1list <- c(M1list, M1)
        M2list <- c(M2list, M2)
        MLElist <- c(MLElist, cc[which.max(dbinom(S[ii], N[ii], f))])
    }
    return(list(M1 = MLElist, SD = sqrt(M2list - M1list^2)))
}


#Aarons shorter version for diploid samples only
adjust.ccf.titan.multi.fordiploidsamples <- function(sampAB, samples,...) {
  for( i in 1:dim(sampAB)[1]) {
    for (j in 1:length(samples)) {
      sn = samples[j]
      print(sn)
      pa1 = as.numeric(sampAB[i, match(paste(sn, "pa", sep=""), colnames(sampAB))])
      nt1 = as.numeric(sampAB[i, match(paste(sn, "nt", sep=""), colnames(sampAB))])
      nb1 = as.numeric(sampAB[i, match(paste(sn, "nb", sep=""), colnames(sampAB))])
      maf1 = as.numeric(sampAB[i, match(paste(sn, "mafc", sep=""), colnames(sampAB))])
      refc1 = as.numeric(sampAB[i, match(paste(sn, "refc", sep=""), colnames(sampAB))])
      altc1 = as.numeric(sampAB[i, match(paste(sn, "altc", sep=""), colnames(sampAB))])
      pu1 = as.numeric(sampAB[i, match(paste(sn, "pu", sep=""), colnames(sampAB))])       #cell purity
      sAGP = as.numeric(sampAB[i, match(paste(sn, "sAGP", sep=""), colnames(sampAB))])    #segmental aneu- ploid genome proportion
      
      print(paste0("assigned mafc, refc and altc from sample ",sn))
      
      CCF1 = computeCCF(maf1,refc1,altc1,pu = 1,pa = 1,sAGP = 1,nt = 2,nb = 1)
      print("finished computing CCF values")
      sampAB[i, match(paste(sn, "ccf", sep=""), colnames(sampAB))] = as.numeric(CCF1[3])
      sampAB[i, match(paste(sn, "ccfSD", sep=""), colnames(sampAB))] = as.numeric(CCF1[4])
      sampAB[i, match(paste(sn, "mafa", sep=""), colnames(sampAB))] = as.numeric(CCF1[1])/2
      sampAB[i, match(paste(sn, "time", sep=""), colnames(sampAB))] = CCF1[2]   #need to change
      print(paste0(sn," ","ccfs added to table"))
    }
  }
  return(sampAB)
}
