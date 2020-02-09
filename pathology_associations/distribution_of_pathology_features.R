rm(list=ls())
library(dplyr);library(tidyr);library(ggplot2);library(stringr)

source("~/aarons_FAP_github_repository/recent_Annotation.R")

path = recent_Annotation(go_online_clinical = "no",go_online_pathology = "no")

path.red = path%>%select(SampleName,Patient,Stage = Stage..Polyp..Normal..AdCa.,Location,Size = Size..mm.,
              VAP_Names,NeoCellsInTumor,nonNeoStromaInTotal,TumorInTotal,
              PercentNormalOverall,PercentStromaOverall,PercentCancerOverall,
              PercentAdenomaOverall,PercentStromaInCancer,PercentStromaInAdenoma,
              PercentNecrosisInTumor)
cols.for.gather = c("NeoCellsInTumor","nonNeoStromaInTotal","TumorInTotal",
                    "PercentNormalOverall","PercentStromaOverall","PercentCancerOverall",
                    "PercentAdenomaOverall","PercentStromaInCancer","PercentStromaInAdenoma",
                    "PercentNecrosisInTumor")
path.red.melt = path.red%>%gather(key = Pathology.Features,value = Pathology.Values,cols.for.gather)
#if size is blank, it is missing
path.red.melt$Size[path.red.melt$Size== ""] = NA
#if size has mutliple dimesions, choose largest
path.red.melt$Size = as.integer(sapply(str_split(path.red.melt$Size,"x",simplify = F),max,1))
head(path.red.melt)

#recategorize sizes to be discrete
path.red.melt$Size.cat[path.red.melt$Size<5 & !is.na(path.red.melt$Size)] = "Small"
path.red.melt$Size.cat[path.red.melt$Size>=5 & path.red.melt$Size<10 & !is.na(path.red.melt$Size)] = "Medium"
path.red.melt$Size.cat[path.red.melt$Size>=10 & !is.na(path.red.melt$Size)] = "Large"

#recategorize Location
table(path.red.melt$Location)
grep(pattern = "ectum",)
path.red.melt$Location

head(path.red.melt,40)
str(path.red.melt)
str(mpg)
p = ggplot(path.red.melt,aes(x = Pathology.Features,y = Pathology.Values))
p+geom_jitter(aes(color = Location))
              
p + geom_dotplot(path.red.melt,
                 mapping = aes(color = Location,size = Size,shape = Stage))



f <- ggplot(mpg, aes(class, hwy))
f + geom_dotplot(binaxis = "y", stackdir = "center"), x, y, alpha, color, fill, group
