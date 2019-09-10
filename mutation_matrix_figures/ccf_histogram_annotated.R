rm(list=ls())
library(dplyr);library(tidyr)
setwd("~/Bulk_A001_A002/")

drivers = read.csv("../../../CancerDriverGenes/PanCanDrivers_Cell2018.csv",
                   sep = ",",skip = 3,header = T)

ccf.table = read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt",
                       header = T,sep = "\t",stringsAsFactors = F)
ccf.table.cols = select(ccf.table,1:5,geneName,geneLoc,functionalClass,ends_with("ccf"),-bpccf,-mergeCCF)
ccf.cols = grep("ccf",colnames(ccf.table.cols))
ccf.cols.name = grep("ccf",colnames(ccf.table.cols),value = T)

#CCF Table of just Driver genes
ccf.table.cols.drivers = ccf.table.cols[(ccf.table.cols$geneName %in% drivers$Gene),]
ccf.table.cols.drivers = filter(ccf.table.cols.drivers,geneLoc=="exonic")%>%
  filter(functionalClass!="synonymous SNV")

#Select Sample CCF
plot.me = ccf.table.cols[,ccf.cols[1]]

names.me = ccf.table.cols.drivers[,c(6,ccf.cols[1])]
nonzero.index = names.me[,2]>0
names.me.nonzero = names.me[nonzero.index,]


#Make Histogram
h = hist(plot.me[plot.me>=0.1],xlim = c(0,1),breaks = 30,ylim = c(0,3000))
h
text.line = mean(h$counts)
driver.genes = 
text(x$hey, y=text.line, labels = x$names,srt = "45")

lines(density(x = plot.me[plot.me>=0.1],
              weights = rep(.05,length(plot.me[plot.me>=0.1]))),
      col = "red",lty = 3)


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
par(yaxs="i",las=1)
hist(air$Respirable.Particles,
     prob=TRUE,col="black",border="white",
     xlab="Respirable Particle Concentrations",
     main="Distribution of Respirable Particle Concentrations")
box(bty="l")
lines(density(air$Respirable.Particles,na.rm=T),col="red",lwd=4)
grid(nx=NA,ny=NULL,lty=1,lwd=1,col="gray")
##############

