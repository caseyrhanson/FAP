rm(list=ls())
library(dplyr);library(tidyr);library(yarrr);library(stringr)
layout(1)

setwd("~/Bulk_A001_A002/")

source("~/aarons_FAP_github_repository/recent_Annotation.R")
recent = recent_Annotation()
recent = recent[recent$DNA_WGS!="",]
recent$SampleName = str_remove_all(recent$SampleName,"-")
# recent$Stage..Polyp..Normal..AdCa.
# recent$Location

drivers = read.csv("../../../CancerDriverGenes/PanCanDrivers_Cell2018.csv",
                   sep = ",",skip = 3,header = T)
driver.colors = data.frame("geneName" = sort(unique(drivers$Gene)),
           "Colors" = piratepal(palette = "basel",length.out = length(unique(drivers$Gene))))
for (pt in c("A001","A002")) {
  ccf.table = read.table(file = paste0("mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt"),
                         header = T,sep = "\t",stringsAsFactors = F)
  ccf.table.cols = select(ccf.table,1:5,geneName,geneLoc,functionalClass,CADD_phred,ends_with("ccf"),-contains("bpccf"),-mergeCCF)
  ccf.cols = which(grepl("ccf",colnames(ccf.table.cols)))
  ccf.cols.name = grep("ccf",colnames(ccf.table.cols),value = T)
  
  #prepare to make at most 12 plots
  
  #CCF Table of just Driver genes
  ccf.table.cols.drivers = ccf.table.cols[(ccf.table.cols$geneName %in% drivers$Gene),]
  ccf.table.cols.drivers = filter(ccf.table.cols.drivers,geneLoc=="exonic")%>%
    filter(functionalClass!="synonymous SNV")
  
  #Save Image
  svg(filename = paste0(pt,"_ccf_Histogram.svg"),width = 12)
  
  par(mfrow=c(3,4),mar = c(4, 4, 3, 1))
  layout.show(n = 12)
  
  for (i in 1:length(ccf.cols)) {
    #Select Sample CCF
    plot.me = ccf.table.cols[,ccf.cols[i]]
    
    names.me = ccf.table.cols.drivers[,c(6,9,ccf.cols[i])]
    nonzero.index = names.me[,ccf.cols.name[i]]>0
    names.me.nonzero = names.me[nonzero.index,]
    names.me.nonzero = names.me.nonzero[names.me.nonzero$CADD_phred>=20,]
    print(length(names.me.nonzero$geneName))
    
    #Make Histogram
    h = hist(plot.me[plot.me>=0.1],
             xlim = c(0,1),breaks = 30,
             ylim = c(0,2000), 
             main = paste0(ccf.cols.name[i]," CCF distribution"),
             xlab = "Cancer Cell Fraction",
             ylab = "Frequency",bty = "l")
    text.line = sort(h$counts,decreasing = T)[2]
    names.me.nonzero$y = text.line + 20*names.me.nonzero$CADD_phred
    
    genes.colors = driver.colors[driver.colors$geneName %in% names.me.nonzero$geneName,]
    names.me.nonzero = full_join(names.me.nonzero,genes.colors,"geneName")
    
    if (nrow(names.me.nonzero)!=0) {
      text(x = names.me.nonzero[,ccf.cols.name[i]],
           y = names.me.nonzero$y,
           labels = names.me.nonzero$geneName,
           offset = .5,
           srt = "45",
           col = as.vector(names.me.nonzero[,"Colors"]))
    }
    next()
    
  }
  dev.off()
}



layout(1)
plot(1:7, abs(stats::rnorm(7)), type = "h", axes = FALSE)
axis(side = 1, at = 1:7, labels = letters[1:7])
box(lty = '1373', col = 'red')
box(which = "plot", lty = "solid")
# lines(density(x = plot.me[plot.me>=0.1],
#               weights = rep(.05,length(plot.me[plot.me>=0.1]))),
#       col = "red",lty = 3)


# CCF Histogram distribution #######
gghistogram(plot.me, x = ccf.cols.name[1], y = "..count..",
            xlab = "Number of citation",
            ylab = "Number of genes",
            binwidth = 5, 
            fill = "lightgray", color = "black",
            label = "gene", label.select = key.gns, repel = TRUE,
            font.label = list(color= "citation_index"),
            xticks.by = 20, # Break x ticks by 20
            gradient.cols = c("blue", "red"),
            legend = c(0.7, 0.6),                                 
            legend.title = ""       # Hide legend title
)
##################



x = data.frame("hey" = as.numeric(rnorm(n = 100,mean = 0.25,sd = .1)),"names" = c("dude","","meh","hey",rep("",96)))
x
h = hist(x$hey,ylim = c(0,50),breaks = 30)
text(x$hey, y=30, labels = x$names,srt = "45")

# lines(density(x$hey,na.rm = T,weights = rep(0.04,length(x$hey))),col = "green",lty = 3)
# density.line = density(x$hey,na.rm = T,weights = rep(0.04,length(x$hey)))
# density.line$
# round(as.numeric(x.dense$x),digits = 2)

################
attach(airquality)
par(yaxs="i",las=1)
hist(Ozone)
layout(c(1,2,3))
hist(Ozone)
plot(density(x = airquality$Ozone,na.rm = T))
hist(Ozone)
box(bty = "l")
lines(density(x = airquality$Ozone,na.rm = T),lwd = 4,col = "red")
layout(1)

hist(air$Respirable.Particles,
     prob=TRUE,col="black",border="white",
     xlab="Respirable Particle Concentrations",
     main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(air$Respirable.Particles,na.rm=T),col="red",lwd=4)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
##############

