#plot HML-2 based on age
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
library(affy)
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

#variables
dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
raw_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
TPM_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures"
metadata = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
tissuemeta = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Full GTEx_HML2_Expression 3.nosync"
Demographics = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
age = paste(figures,"Demographics/Age",sep="/")


#load in meta data
meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE")]
colnames(meta_mod) = c("SUBJID","age")
head(meta_mod)

#load in expression data and save as a list
HML2_filelist = list.files(path = TPM_counts, pattern = "*_TPM_HML2.csv")
HML2_filelist
HML2_filelist <- HML2_filelist[-23]
HML2_datalist = lapply(paste(TPM_counts,HML2_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
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
  
  ERV_Counts_modT_age = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
  head(ERV_Counts_modT_age)
  
  HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT_age
}

HML2_datalist_T
big_data = do.call(rbind, HML2_datalist_T)
head(big_data)
big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                      header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
test = big_data[c("SAMPID","HML2_Sum", "tissue", "age")]
test$age <- as.numeric(test$age)
head(test)
range(test$age)
test$range = ifelse(test$age >= 20 & test$age <= 35, "20-35", ifelse(test$age >=36 & test$age <= 51,"36-51",ifelse(test$age >= 52 & test$age <= 70, "52-70", "NA")))
head(test)

write.csv(test,file=paste(Demographics,"HML2_expression_age_based.csv",sep="/"))

rownames(test) = NULL
test$tissue <- as.factor(test$tissue)
test$range <- as.factor(test$range)
head(test)

table(test$range)
# 20-35  36-51   52-70 

test = drop_na(test)

#for zoomed in boxplot
tissues_of_interest = c("Brain_Cerebellar_Hemisphere", "Brain_Cerebellum","Brain_Cortex","Lung", "Nerve_Tibial","Thyroid","Prostate", "Testis", "Brain_Spinal_cord")
SpecTiss = test[test$tissue %in% tissues_of_interest,]


#get number of doners within each age range for each tissue
tissue_list = unique(test$tissue)
HML2_tissue_list = list()
for (i in (1:length(tissue_list))) {
  HML2_Tissue = as.character(tissue_list[[i]])
  HML2_Tissue
  list_tissue = test[test$tissue %in% HML2_Tissue,]
  df_list_tissue = as.data.frame(table(list_tissue$range))
  head(df_list_tissue)
  df_list_tissue$pct =(df_list_tissue$Freq/sum(df_list_tissue$Freq))*100
  colnames(df_list_tissue) = c("range",paste(HML2_Tissue,'sample_size',sep="_"),paste(HML2_Tissue,'percent',sep="_"))
  HML2_tissue_list[[i]] = df_list_tissue
}
HML2_tissue_list
big_data_tissue = do.call(cbind, HML2_tissue_list)
head(big_data_tissue)
rownames(big_data_tissue)=big_data_tissue$range
drop <- c("range")
big_data_tissue_final = big_data_tissue[ , !(names(big_data_tissue) %in% drop)]
write.csv(big_data_tissue_final,paste(Demographics,"age_range_distribution_for_HML2_expression_differences_07212021.csv",sep="/"))

# ggplot code
p<-ggplot(test, aes(x=tissue, y=HML2_Sum, fill = range)) + geom_boxplot(coef = 6, outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
  labs(title="Age",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path = age, file ="Differences in HML-2 expression in regards to age in GTEx.pdf", width = 8, height = 6)

#zoom in on interesting tissues
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

p<-ggplot(SpecTiss, aes(x=tissue, y=HML2_Sum, fill = range)) + geom_boxplot(position= position_dodge(width = .85))+
  labs(title="Age",x="Tissue", y = "HML2 TPM (sum/sample)")
p + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(path = age, file = "Differences in HML-2 expression in regards to age in GTEx, zoom in_nogrid.pdf")

#geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025)
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
ggsave(path = age, file = "Differences in HML-2 expression in regards to age in GTEx, zoom in with outlier names.png")

#reshape distribution data to plot
plot = big_data_tissue_final[,colnames(big_data_tissue_final)%like% "_percent"]
colnames(plot) = sub("_percent", "", colnames(plot))
plot$range = rownames(plot)
age_distr = reshape2::melt(plot)
colnames(age_distr) = c("range","tissue","freq")

p = ggplot(age_distr, aes(fill=range, y=freq, x=tissue)) + 
  geom_bar(position="stack", stat="identity") + labs(title="Age range sample size in GTEx",x="Tissue", y = "sample size (%)")
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path = age, file = "Age range sample size in GTEx.png")
