#normalize and graph HML-2 expression from GTEx: HML-2 per tissue

library(data.table)
library(Biobase)
library(edgeR)
library(org.Hs.eg.db)
library(DESeq2)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(refGenome)
library(plyr)
library(stringr)
library(ggplot2)
library(matrixStats)
library(dplyr)
library(pheatmap)
library(biomaRt)
library("affy")
library(affycoretools)
library(EnvStats)
library(tidyverse)
library(qusage)
library("Hmisc")
library(gplots)
library("ggpubr")
library(gridExtra)
library(grid)
library(lattice)
library(stringr)

#variables
dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
raw_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
TPM_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures"
metadata = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
setwd(dir)
getwd()

#TPM normalize expression data and extract HML-2, GAPDH, and ACTB counts
#load in individual expression files for each tissue and assign to a list
filelist = list.files(path = raw_counts, pattern = "*V15_Telescope_output.csv") #save data as a list
filelist
datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "_V15_Telescope_output.csv"))
filelist_edited
names(datalist) = filelist_edited
head(datalist)
lapply(datalist,dim)

#convert raw counts to TPM for every data frame in the list
for (i in (1:length(datalist))) {
  Tissue = stringr::str_replace((names(datalist)[[i]]), '_Telescope_output.csv', '')
  Tissue
  df <- datalist[[i]]
  head(df)
  
  #convert to TPM 
  #read in hg38.gtf for transcript coordinates
  ens = ensemblGenome()
  read.gtf(ens, "hg38.gtf")
  class(ens)
  Identifier = "exon"
  hg38_Annotations = extractFeature(ens, Identifier)
  hg38_Annotations
  hg38_Annotation_df = data.frame(start=getGtf(hg38_Annotations)$start,end=getGtf(hg38_Annotations)$end,gene_id=getGtf(hg38_Annotations)$gene_id, transcript_id=getGtf(hg38_Annotations)$transcript_id)
  hg38_Annotation_df[1:5,]
  
  #for each row, get the difference between start and end
  hg38_Annotation_df$difference = hg38_Annotation_df$end-hg38_Annotation_df$start
  hg38_Annotation_df[1:5,]
  
  #for each gene_id, sum the difference to get transcript length. This also converts length to kb.
  #hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(transcript_id), summarise, difference=(sum(difference)/1000)) #FENRIR
  hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(gene_id), summarise, difference=(sum(difference)/1000)) #with Telescope
  hg38_Annotation_df_lengthSum[1:5,]
  
  # Divide the read counts by the length of each gene (transcript?) in kilobases. This gives you reads per kilobase (RPK).
  #stuff = merge(counts, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "transcript_id") #FENRIR
  stuff = merge(df, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "gene_id") #with Telescoppe
  head(stuff)
  tail(stuff)
  stuff_mod = stuff[,-1]
  rownames(stuff_mod) = stuff[,1]
  head(stuff_mod)
  tail(stuff_mod)
  RPK = stuff_mod/(stuff_mod$difference)
  head(RPK)
  
  # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
  Per_Million=(colSums(RPK))/1e6
  Per_Million
  
  # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
  Counts_HML2_TPM <- t(t(RPK)/Per_Million)
  head(Counts_HML2_TPM)[1:5,1:5]
  
  #filter counts file to extract HML-2 and some house keeping genes (GAPDH, ACTB)
  HML2_counts = Counts_HML2_TPM[grepl("^HML-2", rownames(Counts_HML2_TPM)),]
  head(HML2_counts)
  length(rownames(HML2_counts))
  rownames(HML2_counts) <- stringr::str_replace(rownames(HML2_counts), 'HML-2_', '')
  rownames(HML2_counts) <- stringr::str_replace(rownames(HML2_counts), '_new', '')
  head(HML2_counts)
  GAPDH = Counts_HML2_TPM[rownames(Counts_HML2_TPM)=="GAPDH",]
  head(GAPDH)
  ACTB = Counts_HML2_TPM[rownames(Counts_HML2_TPM)=="ACTB",]
  head(ACTB)
  HML2_HKG_counts = rbind(HML2_counts,GAPDH,ACTB)
  head(HML2_HKG_counts)
  rownames(HML2_HKG_counts)
  
  write.csv(HML2_HKG_counts, file=paste(TPM_counts, paste(Tissue, "TPM_HML2.csv", sep="_"), sep = "/"))
  
  #log2(counts+1) transform 
  Counts_HML2_TPM_Filtered_Log2Trans = log2(HML2_HKG_counts + 1)
  head(Counts_HML2_TPM_Filtered_Log2Trans)
  rownames(Counts_HML2_TPM_Filtered_Log2Trans)
  write.csv(Counts_HML2_TPM_Filtered_Log2Trans, file=paste(TPM_counts, paste(Tissue, "TPM_Log2_HML2.csv", sep="_"), sep ="/"))
}

#graph HML-2 expression per tissue in a box plot
HML2_filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
HML2_filelist
HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
HML2_datalist
HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
HML2_filelist_edited
names(HML2_datalist) = HML2_filelist_edited
head(HML2_datalist)
lapply(HML2_datalist,dim)

HML2_datalist_T = list()

for (i in (1:length(HML2_datalist))) {
  HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- HML2_datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
  HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
  HML2_DF$X
  HML2_DF$difference = NULL
  head(HML2_DF)
  colnames(HML2_DF)
  HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
  HML2_DF_SampleSum$X
  HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
  #grab sum row and colnames
  HML2_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
  rownames(HML2_SampleSum) = HML2_SampleSum$X
  HML2_SampleSum$X = NULL
  head(HML2_SampleSum)
  #transpose and add new column with tissue id
  ERV_Counts_modT = as.data.frame(t(HML2_SampleSum))
  ERV_Counts_modT$tissue = HML2_Tissue
  head(ERV_Counts_modT)
  HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT # add it to your list
}
HML2_datalist_T
lapply(HML2_datalist_T,dim)

big_data = do.call(rbind, HML2_datalist_T)
head(big_data)
big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                      header=FALSE, col.names = paste("col", 1:2), stringsAsFactors=FALSE))
big_data$col1 = NULL
rownames(big_data)  = c()
colnames(big_data) = c("HML2_TPM", "tissue")
head(big_data)
test = big_data[c("HML2_TPM", "tissue")]
head(test)
table(test$tissue)
write.csv(table(test$tissue), paste(dir, "Tissue_Distribution.csv", sep ="/"))

#Check values on groups
Cbelum <- test %>% filter(tissue == "Brain_Cerebellum_V15")
Tibnerve <- test %>% filter(tissue == 'Nerve_Tibial_V15')
Prostate <- test %>% filter(tissue == "Prostate_V15")
Thyroid <- test %>% filter(tissue == "Thyroid_V15")
Testis <- test %>% filter(tissue == "Testis_V15")
Corart <- test %>% filter(tissue == "Artery_Coronary_V15")
Esogastro <- test %>% filter(tissue == "Esophagus_Gastroesophageal_Junction_V15")
Summary(Cbelum)

test$tissue <- as.factor(test$tissue)
head(test)

test$tissue <- str_sub(test$tissue, end = -5)

# ggplot code
p<-ggplot(test, aes(x=tissue, y=HML2_TPM)) + geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
  labs(title="Total HML-2 expression in GTEx",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "Total HML-2 expression in GTEx.png", width = 7, height =7 )

#filter TPM expression to >= 1 and graph
TPM_filter = test[test$HML2_TPM >= 1, ]

p<-ggplot(TPM_filter, aes(x=tissue, y=HML2_TPM)) + geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
  labs(title="Total HML-2 expression in GTEx, TPM/sample/tissue >= 1",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "Total HML-2 expression in GTEx, TPM per sample per tissue >= 1.png", width = 7.5, height = 7)
