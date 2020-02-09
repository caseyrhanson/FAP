rm(list=ls())

pt = "A001"
pts = c("A001", "A002", "EP", "JP")

if (pt %in% c("EP","JP")) {
  path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/titan/"
} else if (pt %in% c("A001", "A002")) {
  path = "~/Bulk_A001_A002/titan/segments/"
}
path

seg.files = list.files(path = path,pattern = ".TitanCNA.segments.txt")
seg.files.samples = str_remove_all(string = seg.files, pattern = "_nclones1.TitanCNA.segments.txt")
seg.files.samples = str_remove_all(string = seg.files.samples,pattern = "_nclones2.TitanCNA.segments.txt")

seg.tables = list()

index.seg.files.pt = grep(pattern = pt,x = seg.files)

#Read in each segmentation file and rbind them all together per patient
for (i in 1:length(index.seg.files.pt)) {
  
  temp = read.table(file = paste0(path,seg.files[index.seg.files.pt[i]]),
                    header = T,sep = "\t",stringsAsFactors = F)
  
  temp = temp%>%
    select(chr = chrom, start = loc.start, end = loc.end, copy_number = copynumber, 
           seg.mean = seg.mean)
  
  temp$single_cell_id = seg.files.samples[index.seg.files.pt[i]]

  seg.tables[[i]] = temp
}

big.seg.table.pt = do.call(rbind,seg.tables)
head(big.seg.table.pt)

#remove the EP-AdenoCa and JPAdenoCa samples because the current copy number calls aren't right
big.seg.table.pt = big.seg.table.pt[big.seg.table.pt$single_cell_id!="EP-AdenoCa" & big.seg.table.pt$single_cell_id!="JPAdenoCa",]

x = big.seg.table.pt%>%
  select(chromosome = "chr",start,end,
         segmean =  "copy_number",
         sample = "single_cell_id")
hey = cnFreq(x,
       genome = "hg38",
       plot_title = paste0(pt,": Frequency of CN Gains or Losses for ",length(unique(x$sample)), " samples"),
       plotType = "frequency",x_title_size = 0.8)
ggsave()
#Plot the Copy Number Changes
cnSpec(x,
       genome = "hg38",
       CN_Loss_colour = "blue",CN_Gain_colour = "red",
       plot_title = paste0())

cnView()

library(GenVisR)
# cnSpec()
# cnSpec(LucCNseg, genome="hg19",out = "data")
# cnSpec(LucCNseg, genome="hg38")

# > LucCNseg
# chromosome     start       end probes segmean sample
# 1             1    232500    267500     15    3.31   Luc1
# 2             1    837500   2582500    699    1.06   Luc1
# 3             1   2587500   2630000     18    5.33   Luc1


# 
# # create copy number data
# chromosome <- "chr14"
# coordinate <- sort(sample(0:106455000, size = 2000, replace = FALSE))
# cn <- c(rnorm(300, mean = 3, sd = 0.2), rnorm(700, mean = 2, sd = 0.2), rnorm(1000, 
#                                                                               mean = 3, sd = 0.2))
# data <- as.data.frame(cbind(chromosome, coordinate, cn))
# 
# # create segment data
# dataSeg <- data.frame(chromosome = c(14, 14, 14), start = coordinate[c(1, 301, 
#                                                                        1001)], end = coordinate[c(300, 1000, 2000)], segmean = c(3, 2, 3))
# # call cnView with included segment data
# cnView(data, z = dataSeg, chr = "chr14", genome = "hg19", ideogram_txtSize = 4)
# load("~/Bulk_A001_A002/titan/RData/A001C007.TitanCNA.RData")
# x = titancnaresults[[1]]$results
# x = x%>%select(chromosome = "Chr",coordinate ="Position",cn = "CopyNumber")%>%
#   mutate(chromosome = paste0("chr",chromosome))
# head(x)
# 
# z = read.delim("~/Bulk_A001_A002/titan/segments/A001C007_nclones1.TitanCNA.segments.txt",header = T,sep = "\t")
# z = z%>%
#   select("chromosome" = chrom, "start" = loc.start, "end" =  loc.end, segmean = "copynumber")%>%
#   mutate(chromosome = paste0("chr",chromosome))
# head(z)
# 
# cnView(x,z = z,genome = "hg38", chr = c("chr1","chr2"),CNscale = "absolute")
# GenVisR::hg
# data(hg38chr)





# 
# 
# library(cellscape)
# browseVignettes("cellscape")
# example("cellscape")
# # run cellscape
# # EXAMPLE 2 - COPY NUMBER DATA
# 
# # single cell tree edges
# tree_edges <- read.csv(system.file("extdata", "cnv_tree_edges.csv", 
#                                    cllscp+     package = "cellscape"))
# # tree_edges
# # source target
# # 1      51     50
# # 2      51    c16
# # 3      50    c13
# 
# # cnv segments data
# cnv_data <- read.csv(system.file("extdata", "cnv_data.csv", package = 
#                                    cllscp+     "cellscape"))
# # cnv_data
# # chr     start       end copy_number single_cell_id
# # 1      1         0    977835          NA             c0
# # 2      1    977836   3562166   0.8669376             c0
# # 3      1   3562167   4151406          NA             c0
# # annotations
# sc_annot <- read.csv(system.file("extdata", "cnv_annots.tsv", package = 
#                                    cllscp+     "cellscape"), sep="\t")
# 
# # custom clone colours
# clone_colours <- data.frame( clone_id = c("1","2","3"), 
#                              cllscp+     colour = c("7fc97f", "beaed4", "fdc086"))
# 
# # run cellscape
# cellscape(cnv_data=cnv_data, tree_edges=tree_edges, sc_annot=sc_annot, 
#           cllscp+     width=800, height=475, show_warnings=FALSE, 
#           cllscp+     clone_colours = clone_colours)
# 
# cnv_data[,1:4]
# cellscape(cnv_data = cnv_data,tree_edges = tree_edges)
# cellscape(cnv_data=cnv_data, tree_edges=tree_edges, sc_annot=sc_annot, 
#           width=800, height=475, show_warnings=FALSE,
#           clone_colours = clone_colours)
# 
# # cnv_data
# #chr     start       end copy_number single_cell_id
# #1         0    977835          NA             c0
# #1    977836   3562166   0.8669376             c0
# #1   3562167   4151406          NA             c0


# cnFreq(LucCNseg)
# cnFreq(x)

# data <- cytoGeno[cytoGeno$genome == "hg38", ]
