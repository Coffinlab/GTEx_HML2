#separate HML-2 expression based on biological sex
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
PCA=paste(figures,"PCA",sep="/")
#PCA_data = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/PCA"
Demographics = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
sex = paste(figures,"Demographics/Sex",sep="/")
setwd(dir)
getwd()


#filter males and females based on y chrom
#load in individual expression files for each tissue and assign to a list
filelist = list.files(path =raw_counts, pattern = "*V15_Telescope_output.csv") #save data as a list
filelist
datalist = lapply(paste(raw_counts,filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
datalist

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "_V15_Telescope_output.csv"))
filelist_edited
names(datalist) = filelist_edited
head(datalist)

#convert raw counts to TPM for whole genome counts
for (i in (1:length(datalist))) {
  Tissue = stringr::str_replace((names(datalist)[[i]]), '_Telescope_output.csv', '')
  Tissue 
  df <- datalist[[i]]
  head(df)
  
  #convert to TPM and filter poorly expressed genes (TPM < 1 in at least half of the samples)
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
  hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(gene_id), summarise, difference=(sum(difference)/1000)) #with Telescope
  hg38_Annotation_df_lengthSum[1:5,]
  
  # Divide the read counts by the length of each gene (transcript?) in kilobases. This gives you reads per kilobase (RPK).
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
  write.csv(Counts_HML2_TPM, file=paste(TPM_counts,paste(Tissue, "TPM_total.csv", sep="_"),sep="/"))
  
  #log2(counts+1) transform 
  Counts_HML2_TPM_Filtered_Log2Trans = log2(Counts_HML2_TPM + 1)
  head(Counts_HML2_TPM_Filtered_Log2Trans)
  rownames(Counts_HML2_TPM_Filtered_Log2Trans)
  write.csv(Counts_HML2_TPM_Filtered_Log2Trans, file=paste(TPM_counts,paste(Tissue, "TPM_Log2_total.csv", sep="_"),sep="/"))
}

#plot PCA to make sure they separate properly
counts_filelist = list.files(path = TPM_counts,pattern = "*TPM_Log2_total.csv") #save data as a list
counts_filelist
counts_datalist = lapply(paste(TPM_counts,counts_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
counts_datalist

#name each data frame within the list by the tissue specifier within the file name
counts_filelist_edited = as.character(strsplit(counts_filelist, "_V15_TPM_Log2_total.csv"))
counts_filelist_edited
names(counts_datalist) = counts_filelist_edited
head(counts_datalist)

#load in meta data
meta = read.csv(paste(metadata,"SUBJID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX")]
colnames(meta_mod) = c("SUBJID","age","sex")
meta_mod$sex <- sub("1", "male", meta_mod$sex)
meta_mod$sex <- sub("2", "female", meta_mod$sex)
head(meta_mod)

for (i in (1:length(counts_datalist))) {
  Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_V15', '')
  Tissue
  #read in counts data
  df <- counts_datalist[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  colnames(df) <- gsub("\\.", "-",colnames(df))
  head(df)
  
  #read in run table to filter and merge with meta data
  GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
  GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
  GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == Tissue)
  string = GTEx_sample_ID_tissue$SAMPID
  
  subjid = sub("^(.*?-.*?)-.*", "\\1", string)
  head(subjid)
  length(subjid)
  
  GTEx_sample_ID_tissue$SUBJID = subjid
  
  meta_mod$SUBJID
  meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
  head(meta_updated)
  dim(meta_updated)
  
  
  k <- which(colnames(df) %in% meta_updated$SAMPID)
  df <- df[,k]
  meta_filtered <- meta_updated[match(colnames(df),meta_updated$SAMPID),]
  rownames(meta_filtered) <- colnames(df)
  dim(df)
  
  dim(meta_filtered)
  
  meta_filtered$sex = as.factor(meta_filtered$sex)
  
  obj = ExpressionSet(as.matrix(df))
  phenoData(obj) = AnnotatedDataFrame(meta_filtered)
  obj
  
  #grab rows from counts matrices that correspond to the Y chrom
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
  t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
  table(t2g$chromosome_name)
  # y = 522
  my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
  head(my_ids)
  my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
  head(my_ids_final)
  table(my_ids_final$chromosome_name)
  # y = 70
  
  #add featureData to eset object for Y-specific 
  storageMode(obj) <- "environment"
  featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
  storageMode(obj) <- "lockedEnvironment"
  print(obj)
  
  #extract expression data from counts matrix that originate from the Y chromosome
  listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
  head(listOfIDs)
  Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
  head(Male_exprs)
  dim(listOfIDs)
  
  dim(Male_exprs)
  
  
  #plot PCA for y chromosome based gene expression
  mypath <- file.path(paste(PCA,paste("PCA_V15_", Tissue, ".png", sep = ""), sep ="/"))
  png(file=mypath)
  plotPCA(Male_exprs, groups = as.numeric(pData(obj)$sex), groupnames = levels(pData(obj)$sex), pcs = c(1, 2), main = paste(Tissue,"Male vs Female Y chrom only",sep = " "), )
  dev.off()
}


#Adipose Subcutaneous, Adrenal Gland, Fibroblasts, colon, esophagus mucosa, esophagus muscularis , heart atrial, lung, pancreas, pituitary, skin nse, skin se,small intestine, blood?, breast, stomach?   
#the tissues that have outliers are: Adipose Subcu, Breast Mammary, Cells Fibroblasts, Colon Transverse,
#  Esophagus Mucosa, Heart Atrial Appendage, Lung, Pancreas, Skin NSE, Skin SE, Whole Blood
outlier_list = counts_datalist[c("Adipose_Subcutaneous", "Adrenal_Gland", "Breast_Mammary_Tissue","Cells_Cultured_fibroblasts","Colon_Transverse","Esophagus_Mucosa", "Esophagus_Muscularis",
                                 "Heart_Atrial_Appendage","Heart_Left_Ventricle", "Lung","Pancreas", "Pituitary", "Stomach", "Skin_Not_Sun_Exposed","Skin_Sun_Exposed", "Small_Intestine_Terminal_Ileum", "Whole_Blood")]
outlier_list
names(outlier_list)

for (i in (1:length(outlier_list))) {
  Tissue = stringr::str_replace((names(outlier_list)[[i]]), 'V15_TPM_Log2_total.csv', '')
  Tissue
  #read in counts data
  df <- outlier_list[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  colnames(df) <- gsub("\\.", "-",colnames(df))
  head(df)
  
  GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
  GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
  GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == Tissue)
  string = GTEx_sample_ID_tissue$SAMPID
  
  subjid = sub("^(.*?-.*?)-.*", "\\1", string)
  head(subjid)
  length(subjid)
  
  GTEx_sample_ID_tissue$SUBJID = subjid
  
  meta_mod$SUBJID
  meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
  head(meta_updated)
  dim(meta_updated)
  
  
  k <- which(colnames(df) %in% meta_updated$SAMPID)
  df <- df[,k]
  meta_filtered <- meta_updated[match(colnames(df),meta_updated$SAMPID),]
  rownames(meta_filtered) <- colnames(df)
  dim(df)
  
  dim(meta_updated)
  
  dim(meta_filtered)
  
  
  obj = ExpressionSet(as.matrix(df))
  phenoData(obj) = AnnotatedDataFrame(meta_filtered)
  obj
  
  #grab rows from counts matrices that correspond to the Y chrom
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
  t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
  table(t2g$chromosome_name)
  # y = 522
  my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
  head(my_ids)
  my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
  head(my_ids_final)
  table(my_ids_final$chromosome_name)
  # y = 70
  
  #add featureData to eset object for Y-specific 
  storageMode(obj) <- "environment"
  featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
  storageMode(obj) <- "lockedEnvironment"
  print(obj)
  
  #extract expression data from counts matrix that originate from the Y chromosome
  listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
  head(listOfIDs)
  Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
  head(Male_exprs)
  dim(listOfIDs)
  
  dim(Male_exprs)
  
  
  #get PC data
  stuff = prcomp(Male_exprs)
  test = stuff$rotation
  dim(test)
  
  head(test)
  PC = as.data.frame(test[, -c(3:70)])
  PC$SAMPID = rownames(PC)
  PC_1 = PC[c("SAMPID","PC1")]
  ID_Sex = pData(obj)[c("SAMPID","sex")]
  PC_Sex = join(ID_Sex,PC_1,by = "SAMPID")
  write.csv(PC_Sex, file = paste(PCA,paste(Tissue,"PC1_matrix_for_outlier_removal.csv",sep="_"),sep="/"))
}

#remove outliers
#remove outliers , then redo PCA plot
#outlier run IDs: SRR1317387, SRR1311266, SRR1366790, SRR1436085, SRR1420413, SRR1470313, SRR1479263, SRR1361860, SRR1445064, SRR1436763, SRR1356577 
#outlier sample IDs: GTEX-11ILO-0226-SM-5N9D3, GTEX-1J8EW-1826-SM-C1YQW, GTEX-11ILO-2226-SM-5A5L1, GTEX-11ILO-0008-SM-5QGR9, GTEX-11ILO-2026-SM-5N9CQ, GTEX-11ILO-1026-SM-5A5LP, GTEX-1QCLY-1426-SM-E76P8, GTEX-1I1GT-1526-SM-ARL7K, GTEX-11ILO-2526-SM-5A5LQ, GTEX-11ILO-0726-SM-5HL5I, GTEX-11ILO-1526-SM-5A5KZ"
drop = c("GTEX-11ILO-0226-SM-5N9D3", "GTEX-1J8EW-1826-SM-C1YQW", "GTEX-11ILO-2226-SM-5A5L1", "GTEX-11ILO-0008-SM-5QGR9", "GTEX-11ILO-2026-SM-5N9CQ", "GTEX-11ILO-1026-SM-5A5LP", "GTEX-1QCLY-1426-SM-E76P8", "GTEX-11ILO-0126-SM-5A5LN", "GTEX-11ILO-0626-SM-5A5LO", 
         "GTEX-1I1GT-1526-SM-ARL7K", "GTEX-11ILO-2526-SM-5A5LQ", "GTEX-11ILO-0726-SM-5HL5I", "GTEX-11ILO-1526-SM-5A5KZ", "GTEX-13NYB-0226-SM-5N9G4", "GTEX-11ILO-0006-SM-5LZWE", "GTEX-1GN1W-2126-SM-9JGI7", "GTEX-1PDJ9-1826-SM-E9U66", "GTEX-1PDJ9-0326-SM-DTX95", "GTEX-1PPGY-2926-SM-DPRZF")
# remove outliers and remake eset obj, then re-plot PCA
counts_filelist = list.files(path = TPM_counts, pattern = "*TPM_Log2_total.csv") #save data as a list
counts_filelist
counts_datalist = lapply(paste(TPM_counts,counts_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
counts_datalist

#name each data frame within the list by the tissue specifier within the file name
counts_filelist_edited = as.character(strsplit(counts_filelist, "_V15_TPM_Log2_total.csv"))
counts_filelist_edited
names(counts_datalist) = counts_filelist_edited
head(counts_datalist)

for (i in (1:length(counts_datalist))) {
  Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_V15', '')
  Tissue
  #read in counts data
  df <- counts_datalist[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  colnames(df) <- gsub("\\.", "-",colnames(df))
  head(df)
  df_outliersRemoved = df[,!(names(df) %in% drop)]
  
  #read in run table to filter and merge with meta data
  GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
  GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
  GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
  GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == Tissue)
  string = GTEx_sample_ID_tissue$SAMPID
  
  subjid = sub("^(.*?-.*?)-.*", "\\1", string)
  head(subjid)
  length(subjid)
  
  GTEx_sample_ID_tissue$SUBJID = subjid
  
  meta_mod$SUBJID
  meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
  head(meta_updated)
  dim(meta_updated)
  
  
  k <- which(colnames(df_outliersRemoved) %in% meta_updated$SAMPID)
  df_outliersRemoved <- df_outliersRemoved[,k]
  meta_filtered <- meta_updated[match(colnames(df_outliersRemoved),meta_updated$SAMPID),]
  rownames(meta_filtered) <- colnames(df_outliersRemoved)
  dim(df_outliersRemoved)
  
  dim(meta_updated)
  
  dim(meta_filtered)
  
  meta_filtered$sex = as.factor(meta_filtered$sex)
  
  obj = ExpressionSet(as.matrix(df_outliersRemoved))
  phenoData(obj) = AnnotatedDataFrame(meta_filtered)
  obj
  
  #grab rows from counts matrices that correspond to the Y chrom
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
  t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
  table(t2g$chromosome_name)
  # y = 522
  my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
  head(my_ids)
  my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
  head(my_ids_final)
  table(my_ids_final$chromosome_name)
  # y = 70
  
  #add featureData to eset object for Y-specific 
  storageMode(obj) <- "environment"
  featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
  storageMode(obj) <- "lockedEnvironment"
  print(obj)
  
  #extract expression data from counts matrix that originate from the Y chromosome
  listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
  head(listOfIDs)
  Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
  head(Male_exprs)
  dim(listOfIDs)
  
  dim(Male_exprs)
  
  
  #plot PCA for y chromosome based gene expression
  mypath <- file.path(paste(PCA,paste("PCA_V15_outliers_removed", Tissue, ".png", sep = ""), sep ="/"))
  png(file=mypath)
  plotPCA(Male_exprs, groups = as.numeric(pData(obj)$sex), groupnames = levels(pData(obj)$sex), pcs = c(1, 2), main = paste(Tissue,"Male vs Female Y chrom only",sep = " "))
  dev.off()
}

#plot expression
#now that we know which runs to remove, I'll remove those runs and then plot the expression data
HML2_filelist = list.files(path=TPM_counts, pattern = "*_TPM_HML2.csv") #save data as a list
HML2_filelist
HML2_filelist <- HML2_filelist[-23]
HML2_datalist = lapply(paste(TPM_counts, HML2_filelist, sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
HML2_datalist
HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
HML2_filelist_edited
names(HML2_datalist) = HML2_filelist_edited
head(HML2_datalist)

HML2_datalist_T = list()

for (i in (1:length(HML2_datalist))) {
  HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_V15', '')
  HML2_Tissue
  HML2_df <- HML2_datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_df_1 = HML2_df[,!(names(HML2_df) %in% drop)]
  HML2_DF = HML2_df_1[!grepl("GAPDH", HML2_df_1$X),]
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
  
  ERV_Counts_modT_sex = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
  head(ERV_Counts_modT_sex)
  
  HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT_sex # add it to your list
}

HML2_datalist_T
big_data = do.call(rbind, HML2_datalist_T)
head(big_data)
big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                      header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
test = big_data[c("SAMPID","HML2_Sum", "tissue", "sex")]
head(test)
HML2_outliersRemoved = test[!test$SAMPID %in% drop, ]
head(HML2_outliersRemoved)
dim(HML2_outliersRemoved)

length(HML2_outliersRemoved$SAMPID)


write.csv(HML2_outliersRemoved,file=paste(Demographics,"HML2_outliers_removed_for_sex_differences.csv",sep="/"))
rownames(HML2_outliersRemoved) = NULL

HML2_outliersRemoved$tissue <- as.factor(HML2_outliersRemoved$tissue)
HML2_outliersRemoved$sex <- as.factor(HML2_outliersRemoved$sex)
head(HML2_outliersRemoved)

HML2_outliersRemoved = drop_na(HML2_outliersRemoved)

write.csv(table(HML2_outliersRemoved$tissue), "Tissue_Distribution_outlierDropped.csv")

#checked the number of samples per tissue to make sure the samples were dropped. 
# looks like they were removed. continuing.
table(HML2_outliersRemoved$sex)


#for zoomed in boxplot
tissues_of_interest = c("Brain_Cerebellum","Brain_Cortex","Breast_Mammary_Tissue", "Lung", "Nerve_Tibial","Thyroid")
SpecTiss = HML2_outliersRemoved[HML2_outliersRemoved$tissue %in% tissues_of_interest,]

#get number of males and females for each tissue
tissue_list = unique(HML2_outliersRemoved$tissue)
HML2_tissue_list = list()
for (i in (1:length(tissue_list))) {
  HML2_Tissue = as.character(tissue_list[[i]])
  HML2_Tissue
  list_tissue = HML2_outliersRemoved[HML2_outliersRemoved$tissue %in% HML2_Tissue,]
  df_list_tissue = as.data.frame(table(list_tissue$sex))
  df_list_tissue$Freq_pct = (df_list_tissue$Freq / sum(df_list_tissue$Freq))*100
  head(df_list_tissue)
  colnames(df_list_tissue) = c("sex",paste(HML2_Tissue,'sample_size',sep="_"), paste(HML2_Tissue,'percent',sep="_"))
  HML2_tissue_list[[i]] = df_list_tissue
}
HML2_tissue_list
big_data_tissue = do.call(cbind, HML2_tissue_list)
head(big_data_tissue)
rownames(big_data_tissue)=big_data_tissue$sex
drop <- c("sex")
big_data_tissue_final = big_data_tissue[ , !(names(big_data_tissue) %in% drop)]
write.csv(big_data_tissue_final,paste(Demographics,"sex_distribution_for_sex_differences_07212021.csv",sep="/"))

# ggplot code
pdf(paste(sex,"sex differences HML2 expression.pdf",sep="/"), width = 11, height = 8)
p<-ggplot(test, aes(x=tissue, y=HML2_Sum, fill = sex)) + geom_boxplot(coef = 6, outlier.shape=NA)+
  geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
  labs(title="Biological Sex",x="Tissue", y = "HML2 TPM (sum/sample)")
p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#zoom in on interesting tissues
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#with outlier names
library(ggrepel)
rownames(SpecTiss)=SpecTiss$Run
SpecTiss2 = SpecTiss %>%
  tibble::rownames_to_column(var="outlier") %>%
  group_by(tissue) %>%
  mutate(is_outlier=ifelse(is_outlier(HML2_Sum), HML2_Sum, as.numeric(NA))) 
SpecTiss2$outlier[which(is.na(SpecTiss2$is_outlier))] <- as.numeric(NA)
ggplot(SpecTiss2, aes(x = factor(tissue), y = HML2_Sum)) +
  geom_boxplot() + 
  geom_text_repel(aes(label = outlier), na.rm = TRUE, hjust = -1, size=2,point.padding = NA, segment.color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Differences in HML-2 expression in regards to sex in GTEx, zoom in with outlier names",x="Tissue", y = "sample size")
ggsave(path= sex, file = "sex differences HML2 expression outlier names.png")

p <-ggplot(SpecTiss, aes(x = factor(tissue), y = HML2_Sum, fill = sex)) +
  geom_boxplot(position = position_dodge((width=1))) +
  labs(title="Biological Sex",x="Tissue", y = "HML2 TPM (sum/sample)")
p + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave(path= sex, file = "sex differences HML2 expression zoomwithnsnobackground.pdf")

#reshape distribution data to plot
test = big_data_tissue_final[,colnames(big_data_tissue_final) %like% "_percent"]
colnames(test) = sub("_percent", "", colnames(test))
test$sex = rownames(test)
sex_distr = reshape2::melt(test)
colnames(sex_distr) = c("sex","tissue","pct")

p = ggplot(sex_distr, aes(fill=sex, y=pct, x=tissue)) + 
  geom_bar(position="stack", stat="identity") + labs(title="Biological sex sample size in GTEx",x="Tissue", y = "sample size (%)")
p+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(path = sex, file = "Biological sex sample size in GTEx.png")
