#relationship between hardy score and HML-2 expression
#options(warn=2)
#options(warn=0, error=NULL)
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
library(purrr)


dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
raw_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
TPM_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures"
metadata = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
tissuemeta = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Full GTEx_HML2_Expression 3.nosync"
Demographics = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
Hardy = paste(figures,"Demographics/Hardy",sep="/")
setwd(dir)
getwd()

#load in meta data from Aidan
meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","DTHHRDY")]
colnames(meta_mod) = c("SUBJID","DTHHRDY")
head(meta_mod)

#load in expression data and save as a list
HML2_filelist = list.files(path = TPM_counts, pattern = "*_TPM_HML2.csv")
HML2_filelist
HML2_filelist <- HML2_filelist[-23]
HML2_datalist = lapply(paste(TPM_counts, HML2_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
HML2_datalist
HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
HML2_filelist_edited
names(HML2_datalist) = HML2_filelist_edited
head(HML2_datalist)

HML2_datalist_T = list()

for (i in (1:length(HML2_datalist))) {
  HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- HML2_datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
  HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
  HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
  HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
  HML2_DF$X
  HML2_DF$difference = NULL
  head(HML2_DF)
  colnames(HML2_DF)
  HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
  HML2_DF_SampleSum$X
  HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
  
  #grab sum row and colnames
  HML2_DF_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
  rownames(HML2_DF_SampleSum) = HML2_DF_SampleSum$X
  HML2_DF_SampleSum$X = NULL
  colnames(HML2_DF_SampleSum) <- gsub("\\.", "-",colnames(HML2_DF_SampleSum))
  head(HML2_DF_SampleSum)
  
  #transpose and add new column with tissue id
  ERV_Counts_modT = as.data.frame(t(HML2_DF_SampleSum))
  ERV_Counts_modT$tissue = HML2_Tissue
  ERV_Counts_modT$SAMPID <- rownames(ERV_Counts_modT)
  head(ERV_Counts_modT)
  
  #read in run table to filter and merge with meta data
  GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
  GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
  GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == HML2_Tissue)
  string = GTEx_sample_ID_tissue$SAMPID
  subjid = sub("^(.*?-.*?)-.*", "\\1", string)
  head(subjid)
  length(subjid)
  
  GTEx_sample_ID_tissue$SUBJID = subjid

  meta_mod$SUBJID
  meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
  head(meta_updated)
  dim(meta_updated)
  
  
  ERV_Counts_modT_hardy = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
  head(ERV_Counts_modT_hardy)
  
  HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT_hardy # add it to your list
}

HML2_datalist_T
big_data = do.call(rbind, HML2_datalist_T)
head(big_data)
big_data$SUBJID = NULL
big_data$body_site = NULL
big_data$Sample_Name = NULL
test = big_data[c("SAMPID","HML2_Sum", "tissue", "DTHHRDY")]
head(test)
table(test$DTHHRDY)
 

write.csv(test,file=paste(Demographics,"HML2_expression_HRDY_based.csv",sep="/"))
rownames(test) = NULL

#remove rows with NA under DTHHRDY
test_NAOmit = test[complete.cases(test$DTHHRDY),]
dim(test)

dim(test_NAOmit)


test_NAOmit$tissue <- as.factor(test_NAOmit$tissue)
test_NAOmit$DTHHRDY <- as.factor(test_NAOmit$DTHHRDY)
head(test_NAOmit)

table(test_NAOmit$DTHHRDY)
 

#get number of donors within each Hardy score for each tissue
tissue_list = unique(test_NAOmit$tissue)
HML2_tissue_list = list()
for (i in (1:length(tissue_list))) {
  HML2_Tissue = as.character(tissue_list[[i]])
  HML2_Tissue
  list_tissue = test_NAOmit[test_NAOmit$tissue %in% HML2_Tissue,]
  df_list_tissue = as.data.frame(table(list_tissue$DTHHRDY))
  head(df_list_tissue)
  df_list_tissue$pct = (df_list_tissue$Freq/sum(df_list_tissue$Freq))*100
  colnames(df_list_tissue) = c("DTHHRDY",paste(HML2_Tissue,'sample_size',sep="_"),paste(HML2_Tissue,'_percent',sep="_"))
  HML2_tissue_list[[i]] = df_list_tissue
}
HML2_tissue_list
big_data_tissue = do.call(cbind, HML2_tissue_list)
head(big_data_tissue)
rownames(big_data_tissue)=big_data_tissue$DTHHRDY
drop <- c("DTHHRDY")
big_data_tissue_final = big_data_tissue[ , !(names(big_data_tissue) %in% drop)]
write.csv(big_data_tissue_final,paste(Demographics,"Hardy_score_distribution_for_HML2_expression_differences_07212021.csv",sep="/"))

# ggplot code
p<-ggplot(test_NAOmit, aes(x=tissue, y=HML2_Sum, fill = DTHHRDY)) + geom_boxplot(coef = 6, outlier.shape=NA)  +     
  labs(title="Differences in HML-2 expression in regards to Hardy score in GTEx",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path=Hardy, filename = "Differences in HML-2 expression in regards to Hardy score in GTEx.png")

#for zoomed in boxplot
tissues_of_interest = c("Brain_Cerebellar_Hemisphere", "Brain_Cerebellum","Brain_Cortex","Lung", "Nerve_Tibial","Thyroid","Prostate", "Testis")
SpecTiss = test_NAOmit[test_NAOmit$tissue %in% tissues_of_interest,]


#zoom in on interesting tissues
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
SpecTiss_fil <- SpecTiss %>% subset(DTHHRDY != 0)

p<-ggplot(SpecTiss_fil, aes(x=tissue, y=HML2_Sum, fill = DTHHRDY)) + geom_boxplot(position=position_dodge(width=.85))+
  labs(title="Hardy Score",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path=Hardy, filename = "Differences in HML-2 expression in regards to hardy score in GTEx, zoom in.pdf")


#with outlier names
rownames(SpecTiss)=SpecTiss$Run
SpecTiss2 = SpecTiss %>%
  tibble::rownames_to_column(var="outlier") %>%
  group_by(tissue) %>%
  mutate(is_outlier=ifelse(is_outlier(HML2_Sum), HML2_Sum, as.numeric(NA))) 
SpecTiss2$outlier[which(is.na(SpecTiss2$is_outlier))] <- as.numeric(NA)
ggplot(SpecTiss2, aes(x = factor(tissue), y = HML2_Sum)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE, hjust = -0.9, size=2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Differences in HML-2 expression in regards to age in GTEx, zoom in with outlier names",x="Tissue", y = "sample size")
ggsave(path=Hardy, filename = "Differences in HML-2 expression in regards to age in GTEx, zoom in with outlier names.png")

#reshape distribution data to plot
plot =  big_data_tissue_final[,colnames(big_data_tissue_final)%like%  "_percent"]
colnames(plot) = sub("__percent", "", colnames(plot))
plot$HRDY = rownames(plot)
sex_distr = reshape2::melt(plot)
colnames(sex_distr) = c("HRDY","tissue","freq")

p = ggplot(sex_distr, aes(fill=HRDY, y=freq, x=tissue)) + 
  geom_bar(position="stack", stat="identity") + labs(title="Hardy score range sample size in GTEx",x="Tissue", y = "sample size(%)")
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path=Hardy, filename = "Hardy score range sample size in GTEx.png")
