rm(list=ls())

#http://www.sequenceontology.org/browser/current_svn/term/SO:0001587
#https://www.ncbi.nlm.nih.gov/variation/docs/glossary/


library(g3viz);library(GenomicRanges);library(yarrr)
####
pt = "JP"
gene = "APC"

#collect protein domain information
protein.pfam = g3viz::hgnc2pfam(gene,guess = T,output.format = "list")
protein.domains = protein.pfam$pfam

# Collect 10 different colors for each of the domains
mutant.colors = c(yarrr::piratepal(palette = "basel"),yarrr::piratepal("pony")[6:9])
names(mutant.colors) = unique(protein.domains$hmm.name) #give each unique domain a color

#place the colors into the table according to their matched domain
index.prot.domains.for.colors = match(protein.domains$hmm.name,names(mutant.colors))
protein.domains$color = mutant.colors[index.prot.domains.for.colors]

#add sizes (heights) determined empirically by looking at the graph output
protein.domains$size = c(0.03, 0.02, 0.05, 0.05, 0.04, 0.08, 0.06, 0.08, 0.08, 0.08, 0.08, 0.03, 0.08, 0.08, 0.04, 0.08, 0.02, 0.04, 0.02, 0.08, 0.02, 0.08, 0.08, 0.04, 0.03, 0.03)

# protein.NM_000038 = P25054 = "MAAASYDQLLKQVEALKMENSNLRQELEDNSNHLTKLETEASNMKEVLKQLQGSIEDEAMASSGQIDLLERLKELNLDSSNFPGVKLRSKMSLRSYGSREGSVSSRSGECSPVPMGSFPRRGFVNGSRESTGYLEELEKERSLLLADLDKEEKEKDWYYAQLQNLTKRIDSLPLTENFSLQTDMTRRQLEYEARQIRVAMEEQLGTCQDMEKRAQRRIARIQQIEKDILRIRQLLQSQATEAERSSQNKHETGSHDAERQNEGQGVGEINMATSGNGQGSTTRMDHETASVLSSSSTHSAPRRLTSHLGTKVEMVYSLLSMLGTHDKDDMSRTLLAMSSSQDSCISMRQSGCLPLLIQLLHGNDKDSVLLGNSRGSKEARARASAALHNIIHSQPDDKRGRREIRVLHLLEQIRAYCETCWEWQEAHEPGMDQDKNPMPAPVEHQICPAVCVLMKLSFDEEHRHAMNELGGLQAIAELLQVDCEMYGLTNDHYSITLRRYAGMALTNLTFGDVANKATLCSMKGCMRALVAQLKSESEDLQQVIASVLRNLSWRADVNSKKTLREVGSVKALMECALEVKKESTLKSVLSALWNLSAHCTENKADICAVDGALAFLVGTLTYRSQTNTLAIIESGGGILRNVSSLIATNEDHRQILRENNCLQTLLQHLKSHSLTIVSNACGTLWNLSARNPKDQEALWDMGAVSMLKNLIHSKHKMIAMGSAAALRNLMANRPAKYKDANIMSPGSSLPSLHVRKQKALEAELDAQHLSETFDNIDNLSPKASHRSKQRHKQSLYGDYVFDTNRHDDNRSDNFNTGNMTVLSPYLNTTVLPSSSSSRGSLDSSRSEKDRSLERERGIGLGNYHPATENPGTSSKRGLQISTTAAQIAKVMEEVSAIHTSQEDRSSGSTTELHCVTDERNALRRSSAAHTHSNTYNFTKSENSNRTCSMPYAKLEYKRSSNDSLNSVSSSDGYGKRGQMKPSIESYSEDDESKFCSYGQYPADLAHKIHSANHMDDNDGELDTPINYSLKYSDEQLNSGRQSPSQNERWARPKHIIEDEIKQSEQRQSRNQSTTYPVYTESTDDKHLKFQPHFGQQECVSPYRSRGANGSETNRVGSNHGINQNVSQSLCQEDDYEDDKPTNYSERYSEEEQHEEEERPTNYSIKYNEEKRHVDQPIDYSLKYATDIPSSQKQSFSFSKSSSGQSSKTEHMSSSSENTSTPSSNAKRQNQLHPSSAQSRSGQPQKAATCKVSSINQETIQTYCVEDTPICFSRCSSLSSLSSAEDEIGCNQTTQEADSANTLQIAEIKEKIGTRSAEDPVSEVPAVSQHPRTKSSRLQGSSLSSESARHKAVEFSSGAKSPSKSGAQTPKSPPEHYVQETPLMFSRCTSVSSLDSFESRSIASSVQSEPCSGMVSGIISPSDLPDSPGQTMPPSRSKTPPPPPQTAQTKREVPKNKAPTAEKRESGPKQAAVNAAVQRVQVLPDADTLLHFATESTPDGFSCSSSLSALSLDEPFIQKDVELRIMPPVQENDNGNETESEQPKESNENQEKEAEKTIDSEKDLLDDSDDDDIEILEECIISAMPTKSSRKAKKPAQTASKLPPPVARKPSQLPVYKLLPSQNRLQPQKHVSFTPGDDMPRVYCVEGTPINFSTATSLSDLTIESPPNELAAGEGVRGGAQSGEFEKRDTIPTEGRSTDEAQGGKTSSVTIPELDDNKAEEGDILAECINSAMPKGKSHKPFRVKKIMDQVQQASASSSAPNKNQLDGKKKKPTSPVKPIPQNTEYRTRVRKNADSKNNLNAERVFSDNKDSKKQNLKNNSKVFNDKLPNNEDRVRGSFAFDSPHHYTPIEGTPYCFSRNDSLSSLDFDDDDVDLSREKAELRKAKENKESEAKVTSHTELTSNQQSANKTQAIAKQPINRGQPKPILQKQSTFPQSSKDIPDRGAATDEKLQNFAIENTPVCFSHNSSLSSLSDIDQENNNKENEPIKETEPPDSQGEPSKPQASGYAPKSFHVEDTPVCFSRNSSLSSLSIDSEDDLLQECISSAMPKKKKPSRLKGDNEKHSPRNMGGILGEDLTLDLKDIQRPDSEHGLSPDSENFDWKAIQEGANSIVSSLHQAAAAACLSRQASSDSDSILSLKSGISLGSPFHLTPDQEEKPFTSNKGPRILKPGEKSTLETKKIESESKGIKGGKKVYKSLITGKVRSNSEISGQMKQPLQANMPSISRGRTMIHIPGVRNSSSSTSPVSKKGPPLKTPASKSPSEGQTATTSPRGAKPSVKSELSPVARQTSQIGGSSKAPSRSGSRDSTPSRPAQQPLSRPIQSPGRNSISPGRNGISPPNKLSQLPRTSSPSTASTKSSGSGKMSYTSPGRQMSQQNLTKQTGLSKNASSIPRSESASKGLNQMNNGNGANKKVELSRMSSTKSSGSESDRSERPVLVRQSTFIKEAPSPTLRRKLEESASFESLSPSSRPASPTRSQAQTPVLSPSLPDMSLSTHSSVQAGGWRKLPPNLSPTIEYNDGRPAKRHDIARSHSESPSRLPINRSGTWKREHSKHSSSLPRVSTWRRTGSSSSILSASSESSEKAKSEDEKHVNSISGTKQSKENQVSAKGTWRKIKENEFSPTNSTSQTVSSGATNGAESKTLIYQMAPAVSKTEDVWVRIEDCPINNPRSGRSPTGNTPPVIDSVSEKANPNIKDSKDNQAKQNVGNGSVPMRTVGLENRLNSFIQVDAPDQKGTEIKPGQNNPVPVSETNESSIVERTPFSSSSSSKHSSPSGTVAARVTPFNYNPSPRKSSADSTSARPSQIPTPVNNNTKKRDSKTDSTESSGTQSPKRHSGSYLVTSV"
substr(x = protein.NM_000038,start = 1310,stop = 1310) #checking for A001 germline mutation at 1822, should be V.
# :p.K1310N

#VAP data has a lot of annotations 
if(pt=="A001|A002"){
  vap.path = "~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt"
} else if(any(grepl(pattern = pt,x = c("EP","JP")))){
  vap.path = "~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.simplified.txt"
}
vap = read.delim(file = vap.path,header = T,sep = "\t",stringsAsFactors = F)
#add chr to chromosome if it is not already added
if(!grepl(pattern = "chr",x = vap$chr[1])){
  vap$chr = paste0("chr",vap$chr)
}


# nrow(vap)
# x = vap$geneName
# gene.index = grep(pattern = "APC",x = x)
# vap[gene.index,]
# vap = tbl_df(vap)

vap$geneName = as.character(vap$geneName)


vap.gene = vap%>%filter(geneName=="APC")%>%
  select(1:2,4:6,geneName,geneLoc,functionalClass, AAChange,germline,somatic)%>%
  unite(mut.ID,"chr", "pos", "id", "ref", "alt",sep = ":")

vap.gene$AAChange
vap.gene
vap.gene$mut.ID
#Specific information for the Patient Samples
if(pt=="A001|A002"){
  ccf.path = paste0("~/Bulk_A001_A002/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.ccf.",pt,".txt")
} else if(any(grepl(pattern = pt,x = c("EP","JP")))){
  ccf.path = paste0("~/DNAseq_WGS/scripts/RupingPipelineLocalCopies/post-VAP/mutect.snv.res.filtered.classified.founds.nopara.somatic.table.adjustedCCF.q100.",pt,".txt")
}

ccf.data = read.table(file = ccf.path,header = T,sep = "\t")
ccf.gene = ccf.data%>%filter(geneName==gene)%>%
  select(1:5,functionalClass,somatic)%>%
  unite(mut.ID,"chr", "pos", "id", "ref", "alt",sep = ":")

ccf.gene
#"stopgain" same as "Nonsense"
ccf.gene$functionalClass = str_replace(string = ccf.gene$functionalClass,pattern = "stopgain",replacement = "Nonsense")
ccf.gene
ccf.gene$mut.ID

ccf.gene$mut.ID %in% vap.gene$mut.ID

patient.gene.mutations = left_join(ccf.gene,vap.gene,"mut.ID")

patient.gene.mutations

write.table(x = patient.gene.mutations,file = paste0("~/Bulk_A001_A002/",gene,"_MutationsIn_",pt,".txt"),
            append = F,sep = "\t",row.names = F)



# Get the amino acid change that makes sense for the protein. In my case it was the 2nd
# record int he AAChange column
patient.gene.mutations.AA2 = separate(patient.gene.mutations,col = "AAChange",into = c(NA,"AAChange2"),sep = ",")%>%
  separate("AAChange2",c(NA,NA,NA,NA,"AAChange2"),sep = ":")

#get just the number from the amino acid change
x = str_sub(string = patient.gene.mutations.AA2$AAChange2,start = 4,end = -2)
patient.gene.mutations.AA2$AAChange2.loci = as.numeric(as.character(x))
patient.gene.mutations.AA2$germline = as.character(patient.gene.mutations.AA2$germline)


str(patient.gene.mutations.AA2)

#Manually added Germline Mutation information
if(pt=="A001"){ #A001 Germline mutation
  germline = c("5:112841059:rs459552:T:A","Missense", NA,"APC","exonic","Missense", "p.V1822D", "A001_blood", NA,1822)
  patient.gene.mutations.AA2 = rbind(patient.gene.mutations.AA2,germline)
} else if(pt=="A002"){
  germline = c("5:112838777:rs587779352:ACAAA:","Frameshift",NA,"APC", "exonic","Frameshift", "p.Lys1061_Gln1062insTer","A002_blood",NA,1061)
  patient.gene.mutations.AA2 = rbind(patient.gene.mutations.AA2,germline)
} else if(pt=="EP"){
  germline = c("5:112828001:rs137854572:C:T","Nonsense",NA,"APC","exonic","Nonsense", "p.Gln541X","F_blood",NA,541)
  patient.gene.mutations.AA2 = rbind(patient.gene.mutations.AA2,germline)
  patient.gene.mutations.AA2$somatic.x = str_replace_all(patient.gene.mutations.AA2$somatic.x,pattern = "EP",replacement = "F")
} else if(pt=="JP"){
  germline = c("5:112828001:rs137854572:C:T","Nonsense",NA,"APC","exonic","Nonsense", "p.Gln541X","G_blood",NA,541)
  patient.gene.mutations.AA2 = rbind(patient.gene.mutations.AA2,germline)
  patient.gene.mutations.AA2$somatic.x = str_replace_all(patient.gene.mutations.AA2$somatic.x,pattern = "JP",replacement = "G")
}



######### lolliplot #######
library(Gviz)
library(rtracklayer)
library(trackViewer)


svg(filename = paste0("~/Bulk_A001_A002/proteinFigures/lolliplot_",pt,".svg"),height = 10)

# Add colors to SNP pies by patient.
# https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#plot_multiple_samples

#Make new column with all the sample names only from Somatic and Germline
patient.gene.mutations.AA2$samplesWmut = str_remove_all(patient.gene.mutations.AA2$somatic.x,pattern = "\\[good]\\,|\\[goodSub]\\,")
index.na.somatic = which(is.na(patient.gene.mutations.AA2$samplesWmut))
#replace the NA in the somatic column with the corresponding value in the germline column
patient.gene.mutations.AA2$samplesWmut[index.na.somatic] = patient.gene.mutations.AA2$germline[index.na.somatic]

#only keep samples with AAmutations for the plots
patient.gene.mutations.AA2 = patient.gene.mutations.AA2%>%filter(!is.na(AAChange2))


snp.labels = paste(patient.gene.mutations.AA2$AAChange2,
                   patient.gene.mutations.AA2$functionalClass.x,
                   patient.gene.mutations.AA2$samplesWmut)
SNP = patient.gene.mutations.AA2$AAChange2.loci
SNP = SNP[!is.na(SNP)]
SNP = as.numeric(SNP)

snp.labels = snp.labels[!grepl(pattern = "NA",x = snp.labels)]

sample.gr = GRanges("chr5",
                    IRanges(SNP,width=1,names = snp.labels),
                    GermOrSom = ifelse(!is.na(patient.gene.mutations.AA2$somatic.x),"somatic","germline"))

#add border color
sample.gr$color <- ifelse(sample.gr$GermOrSom=="somatic","lightgrey","black")

features = GRanges("chr5",
                   IRanges(start = protein.domains$start,
                                  end = protein.domains$end,
                                  names = protein.domains$hmm.name),
                   color = "black", fill = protein.domains$color, height = protein.domains$size)

lolliplot(sample.gr, features,type = "circle")

dev.off()

############################