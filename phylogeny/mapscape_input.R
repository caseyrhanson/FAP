

library("mapscape");library(dplyr);library(data.tree);library(ape);library(tidyr)

# EXAMPLE 1 - Patient A21, Gundem et al., 2015

# clonal prevalences
clonal_prev <- read.csv(system.file("extdata", "A21_clonal_prev.csv", 
                                    package = "mapscape"))
head(clonal_prev)
# sample_id clone_id  clonal_prev
# 1         A        1  0.001043161
# 2         A        2 -0.013192070
# 3         A        3  0.000371804
# mutations
mutations <- read.csv(system.file("extdata", "A21_mutations.csv", 
                                  package = "mapscape"))
mutations%>%arrange(chrom)%>%arrange(coord)
head(mutations%>%arrange(chrom)%>%arrange(coord))
# sample_id chrom  coord clone_id       VAF
# 1         G     5 434261        5 0.5000000
# 2         E     5 434261        5 0.2857143
# 3         C     5 434261        5 0.3750000
# 4         H     5 434261        5 0.6111111
# 5         D     5 434261        5 0.0000000
# 6         J     5 434261        5 0.0000000
# locations of each tumour sample on user-provided image
sample_locations <- read.csv(system.file("extdata", 
                                         "A21_sample_locations.csv", package = "mapscape"))
head(sample_locations)
# sample_id location_id   x   y
# 1         A           A  79 297
# 2         C           C 177 306
# 3         D           D 101 317
# 4         E           C 177 306
# 5         F           F 245 307
# 6         G           C 177 306

# genotype tree edges
tree_edges <- read.csv(system.file("extdata", "A21_tree.csv", package 
                                   = "mapscape"))
head(tree_edges)
# source target
# 1 ancestor      2
# 2        2      4
# 3        2      8
# 4        8      9
# 5        8     18
# 6       18      1

# image reference
img_ref <- system.file("extdata", "A21_anatomical_image.png", package 
                       = "mapscape")

# radial order of samples
sample_ids <- c("H","F","J","D","A","I","C","E","G")

# run mapscape
mapscape(clonal_prev = clonal_prev, tree_edges = tree_edges, 
         sample_locations = sample_locations, mutations = mutations, 
         img_ref = img_ref, sample_ids = sample_ids)

#########################

# tree_edges
tree = ape::read.tree(file = "~/Bulk_A001_A002/phylip/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.A001.txt.fasta.phylip")
x = data.tree::as.Node(tree)
tree_edges = ToDataFrameNetwork(x)
colnames(tree_edges)
tree_edges = select(tree_edges,source = "from", target = "to")

# clonal_prev = clonal_prev
# sample_id clone_id  clonal_prev
# 1         A        1  0.001043161
# 2         A        2 -0.013192070
# 3         A        3  0.000371804
clusters = read.delim(file = "~/Bulk_A001_A002/pyclone/A001/tables/cluster.tsv")
colnames(clusters)
head(clusters)
clonal_prev = select(clusters,sample_id,clone_id="cluster_id",clonal_prev="mean")

# sample_locations = sample_locations, 
# sample_id location_id   x   y
# 1         A           A  79 297
# 2         C           C 177 306
# 3         D           D 101 317
# 4         E           C 177 306
# 5         F           F 245 307
# 6         G           C 177 306
sample_id = sort(unique(clonal_prev$sample_id))
location_id = LETTERS[1:length(sample_id)]
x = 1:length(sample_id)
y = sample(x = 1:20,length(sample_id),replace = F)
sample_locations = data.frame(sample_id,location_id,x,y)

# mutations = mutations, 
# sample_id chrom  coord clone_id       VAF
# 1         G     5 434261        5 0.5000000
# 2         E     5 434261        5 0.2857143
# 3         C     5 434261        5 0.3750000
# 4         H     5 434261        5 0.6111111
# 5         D     5 434261        5 0.0000000
# 6         J     5 434261        5 0.0000000
x = read.delim(file = "~/Bulk_A001_A002/pyclone/A001/tables/loci.tsv")
head(x)
mutations = separate(x,col = mutation_id,into = c("chrom","coord",NA,NA,NA),sep = ":")%>%
  select(sample_id,chrom,coord,clone_id="cluster_id",VAF = "variant_allele_frequency")


# img_ref = img_ref
img_ref <- system.file("extdata", "A21_anatomical_image.png", package 
                       = "mapscape")
img_ref = "~/Google Drive/Stanford_postdoc/Research/HTAN/Aaron_Horning_CurtisLab.png"

# sample_ids = sample_ids
sample_ids = sort(unique(clonal_prev$sample_id))


# run mapscape
mapscape(clonal_prev = clonal_prev, tree_edges = tree_edges, 
         sample_locations = sample_locations, mutations = mutations, 
         img_ref = img_ref, sample_ids = sample_ids)

