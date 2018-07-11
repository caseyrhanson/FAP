


# Read in libraries and set seed

rm(list=ls())
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenVisR") #already installed on Aaron's macbookpro
# Load the GenVisR package
library("GenVisR");library(dplyr);library(tidyr);library(ggplot2);library(stringr);library(reshape2)
set.seed(426)


setwd("~/Google Drive/Stanford_postdoc/Research/FAP Samples_EdEsplin/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP")




### Read in and organize main datasets FAP samples: clinical and mutation

        #Read in clinical information for adding additional information to waterfall plots
        sampleData <- read.table("https://docs.google.com/spreadsheets/d/e/2PACX-1vRd0tND0K2rcfq_C6pte4_W42vkuwsemkv_A_jnrm6mm5N72nxxw5Bezumy898XmtzEQlcxfrYW9j2E/pub?output=csv",sep = ",",header = T)
        
        ### Read in and reorganize the Somatic Mutation table. Prepare for Waterfall plot creation.
        file<-read.table(file = "mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified", header = T, sep = "\t",quote = "")



### For loops to analyze different Exonic regions and Patient data.
patients <- c("JP","EP")
regions <- c("exonic","non-exonic")
for (pt in patients) {
        for (region in regions) {
                #choose only necessary column data. Replace "-" with "."
                #Chose patient EP to start with
                clinical<-select(sampleData, "Sample_Name", "Stage", "Location_of_Sampled_Specimen", "Size_.mm.", "Size_SML.S..5.5.M..10.L.10.", "Grade", "PatientSample_ez",  "mean_coverageMetrics")%>%
                        filter(grepl(paste(pt),PatientSample_ez))%>%
                        mutate(sample=str_replace_all(PatientSample_ez, pattern = "-",replacement = "."))%>%
                        select(-Sample_Name, Grade = Grade, Size = "Size_SML.S..5.5.M..10.L.10." , Location = Location_of_Sampled_Specimen , Stage = Stage)
                #Create the long melted version of the clinical table
                clinical.melt<-melt(clinical,id.vars = "sample")
                
                target <- c( "Size", "Stage")
                clinical.melt<-filter(clinical.melt, variable %in% target)
                
                
                
                #Select only statistically interesting columns. Rename the columns so that are compatible with sample information
                #also, create new columns that represent the number of reads that a variant is represented by: (maf x depth = varreads).
                file_2<-select(file, chr, pos, id, ref, alt, geneName, geneLoc, functionalClass, somatic, germline,ends_with(match = "maf"), ends_with(match = "d"),-CADD_phred, -Polyphen2_HVAR_pred, -mutect.snv.res.filtered.classified.founds.flanking.bam.out.bad)
                
                # maf <- select(file, contains('maf')) #select only "maf" minor allele frequency columns
                # d <- select(file, ends_with('d'),-id,) #select only "d" depth columns
                # num_var_reads <- maf*d #multiplying the two determines the number of reads which contained the variant.
                # colnames(num_var_reads)<-str_replace_all(colnames(num_var_reads),pattern = "maf",replacement = "") #relabel the columns for the num_var_reads
                # 
                # file<-cbind(file,num_var_reads) #adds num_var_reads columns to file table
                # file_varreads<-file[,c(1:10,71:100)] #selects variant descriptor columns and num_var_reads
                
                
                file_varreads<-file_2[,c(1:40)] #selects variant descriptor columns and maf columns
                colnames(file_varreads)<-str_replace_all(colnames(file_varreads),pattern = "maf",replacement = "")
                
                #select only variants per polyp which are known to be somatic
                file_varreads_melt<-gather(file_varreads,"SampleName", "MAF",11:40)%>%
                        filter(grepl(paste(pt),SampleName))%>%
                        filter(str_detect(somatic, paste0(SampleName,"\\","[")))%>% #filters for only somatic mutations per sample
                        filter(MAF>=0.05)
                
                # #Create new column for variant class - if in geneLoc exon, add the functional class.
                if (region=="exonic") {
                        #Create table with only exon functional variations (remove the non-exonic variants)
                        file_varreads_melt$variant_class <- file_varreads_melt$functionalClass
                        file_varreads_melt_2<-filter(file_varreads_melt, variant_class!="NA") 
                } else if (region=="non-exonic") {
                        #Create table without exon functional variations (remove the exonic variants)
                        file_varreads_melt$variant_class <- file_varreads_melt$functionalClass
                        file_varreads_melt_2<-file_varreads_melt[is.na(file_varreads_melt$variant_class),]
                }
                
                
                #Making custom MAF plots per patient. make each varialbe a factor
                customMAF <- select(file_varreads_melt_2, "sample"="SampleName","gene"="geneName", variant_class )%>% #sample,gene,variant_class
                        mutate_all(funs(factor))
                
                #
                if (pt=="EP") {
                        #EP samples to plot - in order
                        samples_to_plot <- c("EP_Dec_NL", "EP.11",  "EP.31", "EP.57",  "EP.84","EP.26","EP.6", "EP.74B", "EP.88B", "EP.AdenoCa")        
                } else if (pt=="JP") {
                        #JP samples to plot - in order
                        samples_to_plot <- c("JP_Dec_NL.1", "JP69", "JP9",  "JP38","JP3", "JP31", "JP63", "JP34B", "JP61B", "JP6B", "JPAdenoCa")        
                }
                
                
                
                #for exonic regions only
                main_layer <- theme_grey()
                clinical_layer <- theme_classic()+
                        theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10))
                
                #pdf(file="EP_exonicvariants_waterfall.pdf", height=20, width=15)
                pdf(file=paste0(pt,"_",region,"_variants_waterfall.pdf"), height=20, width=15)
                waterfall(customMAF, fileType = "Custom",sampOrder = samples_to_plot, variant_class_order = c("nonsynonymous SNV", "stopgain", "synonymous SNV","unknown"), mainDropMut = T,mainXlabel = T,mainRecurCutoff = .03,plotMutBurden = T,mainPalette = c("purple","red", "green", "grey" ), mainLayer = main_layer, clinData = clinical.melt, clinLayer = clinical_layer, clinLegCol = 2)
                dev.off()
                
                
                #determine which genes are mutated more than once
                length(customMAF$gene) #Total number of mutations found in patient samples
                length(unique(customMAF$gene)) #number of distinct mutated genes in patient samples
                repeat_gene<-group_by(file_varreads_melt_2,geneName)%>%
                        arrange(geneName)%>%
                        count(geneName)%>%
                        arrange(geneName)%>%
                        filter(n>1)%>%
                        arrange(geneName)
                        
                #creates table with shared mutated genes.
                table <- inner_join(x = file_varreads_melt_2, y = repeat_gene,"geneName")%>%
                        group_by(geneName)%>%
                        arrange(pos)%>%
                        arrange(geneName)
                write.csv(table,file = paste0(pt,"_",region,"_sharedmuts.csv"),append = F,col.names = T,sep=",",row.names = F)
                
                #creates table with shared point mutations
                table_2 <- add_count(table,pos, id, sort = T)%>%
                        filter(nn>1)
                write.csv(table_2,file = paste0(pt,"_",region,"_sharedmuts_exact.csv"),append = F,col.names = T,sep=",",row.names = F)
                
                # ######### Germline mutations #############
                # 
                # APC<-filter(file, chr==5,grepl(pattern = "APC",x = geneName))
                # col_APC<-colnames(APC)
                #         select(APC,col_APC,ends_with("maf"))
                # 
                # file_varreads_melt<-gather(file_varreads,"SampleName", "MAF",11:40)%>%
                #         filter(grepl("EP",SampleName))%>%
                #         filter(str_detect(somatic, paste0(SampleName,"\\","[")))%>% #filters for only somatic mutations per sample
                #         filter(MAF>=0.05)



        }
}
        