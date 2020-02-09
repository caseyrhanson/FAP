rm(list=ls())


# Libraries
library(tidyverse)
library(stringr)
library(ggplot2);library(RColorBrewer)

setwd("~/Bulk_A001_A002/schematic_colons")


############ A001 #########
pt = "A001"
# path = "https://docs.google.com/spreadsheets/d/e/2PACX-1vSpRgOtkF7iQsDyG1_WUmzvXvE_4DsrUrISQsIA64403vI6n-YpSkWeYTs24f8LEWi0dwlOmHdWYMXT/pub?output=csv"
# mergedlocations = read.table(file=path,header = T,sep = ",")
# write.table(mergedlocations,file = paste0(pt,"_ColonSchematic_Polypslocations.txt"),append = F,sep = "\t",row.names = F)
mergedlocations = read.table(file = paste0(pt,"_ColonSchematic_Polypslocations.txt"),header = T,sep = "\t")

#Polyp values
polyps.locs = mergedlocations[!is.na(mergedlocations$Polyp),]
#Border values
border.locs = mergedlocations[is.na(mergedlocations$Polyp),]


#Create Vector for colors
colors = brewer.pal(n = length(unique(polyps.locs$Location)),name = "Set1")
locs = unique(polyps.locs$Location)
color.vector = vector()
for (i in 1:length(locs)) {
  x = rep(colors[i],sum(polyps.locs$Location == locs[i]))
  color.vector = append(color.vector,x)
}
color.vector

pdf(file = "Schematic_A001.pdf",width = 12,height = 7)
plot(polyps.locs$X_adjustedcoordinate,
     polyps.locs$Y_adjustedcoordinate,
     pch = 20,ylab = "", col = color.vector,bty = "n",
     xlab = "Descending            Transverse            Ascending",
     main = "Schematic of Patient A001's Polyp-Filled Colon", ylim = c(-8, 6))

text(polyps.locs$X_adjustedcoordinate,
     polyps.locs$Y_adjustedcoordinate,
     polyps.locs$Polyp,cex = 0.8,pos = 4,srt = 315)

#Add Borders
top = border.locs[border.locs$Location=="top",]
bottom = border.locs[border.locs$Location=="bottom",]
lines(top$X_adjustedcoordinate,top$Y_adjustedcoordinate,xpd = NA)
lines(sort(bottom$X_adjustedcoordinate),bottom$Y_adjustedcoordinate, xpd = NA)
dev.off()

# collect X and Y cartesian coordinates
distances.x.y = as.matrix(cbind(polyps.locs$X_adjustedcoordinate,polyps.locs$Y_adjustedcoordinate))
rownames(distances.x.y) = polyps.locs$Polyp

#calculate Euclidean distances and create matrix format then dataframe format
physical.distance = dist(distances.x.y,method = "euc")
physical.distance = as.matrix(physical.distance)
physical.distance = data.frame(Sample = rownames(physical.distance),physical.distance)
colnames(physical.distance) = str_remove_all(string = colnames(physical.distance),pattern = "X")

#melt the table down so that the Euc dist is shown next to its paired samples
physical.distance = physical.distance%>%gather(key = SampleID,value = Euc_dist,-Sample)

write.table(physical.distance,
            file = paste0(pt,"_ColonSchematic_Polypslocations_distances.txt"),append = F,sep = "\t")




############ A002
rm(list=ls())
pt = "A002"
# path = "https://docs.google.com/spreadsheets/d/e/2PACX-1vSTnvwx9bM1nNlADvvUgKdvIRGWJhjw5iXP875HflluhhDnlXziZnXvhVUAF5HK1BWtaX-r5WMktdwe/pub?gid=661794192&single=true&output=csv"
# mergedlocations = read.table(file=path,header = T,sep = ",")
# write.table(x = mergedlocations,file = paste0(pt,"_ColonSchematic_Polypslocations.txt"),append = F,sep = "\t",row.names = F)
mergedlocations = read.table(paste0(pt,"_ColonSchematic_Polypslocations.txt"),header = T,sep = "\t")

border.bool = grepl(pattern = "B",mergedlocations$Polyp)

#Polyp values
polyps.locs = mergedlocations[!border.bool,]
#Border values
border.locs = mergedlocations[border.bool,]
border.locs = border.locs[order(border.locs$X_Coordinate),]

#Create Vector for colors
colors = brewer.pal(n = length(unique(polyps.locs$Location)),name = "Set1")
locs = unique(polyps.locs$Location)
color.vector = vector()
for (i in 1:length(locs)) {
  x = rep(colors[i],sum(polyps.locs$Location == locs[i]))
  color.vector = append(color.vector,x)
}
color.vector

pdf(file = "Schematic_A002.pdf",width = 12,height = 7)
plot(polyps.locs$X_Coordinate,
     polyps.locs$Y_Coordinate,
     pch = 20,ylab = "", col = color.vector,bty = "n",
     xlab = "Descending            Transverse            Ascending",
     main = "Schematic of Patient A001's Polyp-Filled Colon", ylim = c(0, 6),xlim = c(0,50))

text(polyps.locs$X_Coordinate,
     polyps.locs$Y_Coordinate,
     polyps.locs$Polyp,cex = 0.8,pos = 4,srt = 315)

#Add Borders
top = border.locs[border.locs$Location=="top",]
top$Y_Coordinate[top$Y_Coordinate >= 6] = 6
bottom = border.locs[border.locs$Location=="bottom",]
bottom$Y_Coordinate[bottom$Y_Coordinate <= 0] = 0
lines(top$X_Coordinate,top$Y_Coordinate,xpd = NA)
lines(sort(bottom$X_Coordinate),bottom$Y_Coordinate, xpd = NA)


dev.off()



####### Calculate locations between polyps
distances.x.y = as.matrix(cbind(polyps.locs$X_Coordinate,
                                polyps.locs$Y_Coordinate))
rownames(distances.x.y) = polyps.locs$Polyp

#calculate Euclidean distances and create matrix format then dataframe format
physical.distance = dist(distances.x.y,method = "euc")
physical.distance = as.matrix(physical.distance)
physical.distance = data.frame(Sample = rownames(physical.distance),physical.distance)
colnames(physical.distance) = str_remove_all(string = colnames(physical.distance),pattern = "X")

#melt the table down so that the Euc dist is shown next to its paired samples
physical.distance = physical.distance%>%gather(key = SampleID,value = Euc_dist,-Sample)
physical.distance
write.table(physical.distance,
            file = paste0(pt,"_ColonSchematic_Polypslocations_distances.txt"),append = F,sep = "\t")
