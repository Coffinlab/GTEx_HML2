#graph heatmap of individual HML-2 proviruses
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

#variables
dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
TPM_counts =  "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures =  "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures/"
heatmap = 
setwd(dir)
getwd()

#graph heatmap of individual HML-2 proviruses
#load in individual expression files for each tissue and assign to a list
filelist = list.files(path = TPM_counts, pattern = "*TPM_HML2.csv") #save data as a list
filelist
datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
datalist

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "_TPM_HML2.csv"))
filelist_edited
names(datalist) = filelist_edited

HML2_datalist_AvgStDev = list()

#calculate average and stdev for each provirus
for (i in (1:length(datalist))) {
  HML2_Tissue = stringr::str_replace((names(datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_df$X
  HML2_df$difference = NULL
  head(HML2_df)
  
  HML2_df$Average = rowMeans(HML2_df[,c(-1)])
  HML2_df$StDev = rowSds(as.matrix(HML2_df[,c(-1)]))
  HML2_DF_AvgStDev <- HML2_df[, c("X", "Average", "StDev")]
  colnames(HML2_DF_AvgStDev) = c("provirus", "Average", "StDev")
  head(HML2_DF_AvgStDev)
  
  #assign provirus names as the row names and remove the old column containing just provirus names
  rownames(HML2_DF_AvgStDev) = HML2_DF_AvgStDev$provirus
  HML2_DF_AvgStDev$provirus = NULL
  head(HML2_DF_AvgStDev)
  
  HML2_datalist_AvgStDev[[HML2_Tissue]] <- HML2_DF_AvgStDev # add it to your list
}
HML2_datalist_AvgStDev
big_data = do.call(rbind, HML2_datalist_AvgStDev)
head(big_data)
big_data$test = rownames(big_data)
big_data$provirus = gsub("^.*?\\.","", big_data$test)
big_data$tissue = gsub("[.][\\s\\S]*$", "", big_data$test, perl = T)
rownames(big_data)  = c()
big_data$test = NULL
head(big_data)
big_data_transf = big_data[,c(3,1,2,4)]
head(big_data_transf)
big_data_transf_mod = big_data_transf[big_data_transf$Average >= 1, ]
for_VLP_analysis = reshape2::dcast(data = big_data_transf_mod,formula = provirus~tissue,fun.aggregate = sum,value.var = "Average")

site_average <- for_VLP_analysis[-c(41:42),]
site_average <- site_average[-c(10,22,39),]
site_average[38,2:55] <- colMeans(site_average[,2:55])
site_average[38,1] <- "Average"
write.csv(site_average, "site_average.csv")

site_sum <- for_VLP_analysis[-c(41:42),]
site_sum <- site_sum[-c(10,22,39),]
site_sum[38,2:55] <- colSums(site_sum[,2:55])
site_sum[38,1] <- "Sum"
write.csv(site_sum, "site_sum.csv")

#save data 
write.csv(big_data_transf, 'HML2_individual_expression_02042022.csv')
write.csv(for_VLP_analysis, "HML2_individual_expression_wide_form_02042022.csv")

#get data for heatmap
HML2_heatmap = big_data_transf[, c("provirus", "Average", "tissue")]
head(HML2_heatmap)
HML2_heatmap_1 = reshape2::dcast(HML2_heatmap, provirus ~ tissue, value.var = "Average", fun.aggregate = mean)
rownames(HML2_heatmap_1) = HML2_heatmap_1$provirus
HML2_heatmap_1$provirus = NULL
rownames(HML2_heatmap_1)
HML2_heatmap_2 = HML2_heatmap_1[c("ACTB","GAPDH"),]
head(HML2_heatmap_2)
rownames(HML2_heatmap_2)
remove = c("ACTB","GAPDH")
HML2_heatmap_3 = HML2_heatmap_1[!row.names(HML2_heatmap_1)%in%remove,]
head(HML2_heatmap_3)
rownames(HML2_heatmap_3)
HML2_heatmap_final=rbind(HML2_heatmap_3,HML2_heatmap_2)
head(HML2_heatmap_final)
rownames(HML2_heatmap_final)
HML2_heatmap_final <- HML2_heatmap_final[-c(19,42,83),]
write.csv(HML2_heatmap_final, 'HML2_individual_expression_heatmap_protrim.csv')

#reorganize provirus order by ascending values
order = c("1p31.1a", "1p31.1b", "1p34.3", "1p36.21a", "1p36.21c", "1q21.3", "1q22", "1q23.3", "1q24.1", "1q32.2", "1q43", "2q21.1", "3p12.3", "3p25.3", "3q12.3", "3q13.2", "3q21.2", "3q24", "3q27.2", "4p16.1a", "4p16.1b", "4p16.3a", "4p16.3b", "4q13.2", "4q32.1", "4q32.3", "4q35.2", "5p12", "5p13.3", "5q33.2", "5q33.3", "6p11.2", "6p21.1", "6p22.1", "6q14.1", "6q25.1", "7p22.1a", "7p22.1b", "7q11.21", "7q22.2", "7q34", "8p22", "8p23.1a", "8p23.1b", "8p23.1c", "8p23.1d", "8q11.1", "8q24.3a", "8q24.3b","9q34.11", "9q34.3", "10p12.1", "10p14", "10q24.2", "11p15.4", "11q12.1", "11q12.3", "11q22.1", "11q23.3", "12p11.1", "12q13.2", "12q14.1", "12q24.11", "12q24.33", "14q11.2", "14q32.33", "15q25.2", "16p11.2", "16p13.3", "19p12a", "19p12b", "19p12c", "19p12d", "19p12e", "19p13.3", "19q11", "19q13.12b", "19q13.41", "19q13.42", "20q11.22", "22q11.21", "22q11.23", "Xq12", "Xq21.33", "Xq28a", "Xq28b", "Yp11.2", "Yq11.23a", "Yq11.23b")
HML2_heatmap_3$provirus = rownames(HML2_heatmap_3)
HML2_heatmap_4 = HML2_heatmap_3 %>%
  dplyr::slice(match(order,provirus))
rownames(HML2_heatmap_4)=HML2_heatmap_4$provirus
HML2_heatmap_4$provirus = NULL
HML2_heatmap_3$provirus = NULL

# plot a basic heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "blue"))(paletteLength)

pdf("Average individual HML-2 expression in GTEx, with housekeeping genes.pdf")
map = pheatmap(as.matrix(HML2_heatmap_final), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
               main = "Average individual HML-2 expression in GTEx, with housekeeping genes", fontsize = 8, fontsize_row = 6, 
               fontsize_col = 6)
dev.off()

pdf("Average individual HML-2 expression in GTEx.pdf")
map = pheatmap(as.matrix(HML2_heatmap_4), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
               main = "Average individual HML-2 expression in GTEx",breaks = seq(0,8,by=0.2),fontsize = 8, fontsize_row = 6, 
               fontsize_col = 6)
dev.off()

order = c("1p31.1a", "1p31.1b", "1p34.3", "1p36.21a", "1p36.21c", "1q21.3", "1q22", "1q23.3", "1q24.1", "1q32.2", "1q43", "2q21.1", "3p12.3", "3p25.3", "3q12.3", "3q13.2", "3q21.2", "3q24", "3q27.2", "4p16.1a", "4p16.1b", "4p16.3a", "4p16.3b", "4q13.2", "4q32.1", "4q32.3", "4q35.2", "5p12", "5p13.3", "5q33.2", "5q33.3", "6p11.2", "6p21.1", "6p22.1", "6q14.1", "6q25.1", "7p22.1a", "7p22.1b", "7q11.21", "7q22.2", "7q34", "8p22", "8p23.1a", "8p23.1b", "8p23.1c", "8p23.1d", "8q11.1", "8q24.3a", "8q24.3b", "9q34.11", "9q34.3", "10p12.1", "10p14", "10q24.2", "11p15.4", "11q12.1", "11q12.3", "11q22.1", "11q23.3", "12p11.1", "12q13.2", "12q14.1", "12q24.11", "12q24.33", "14q11.2", "14q32.33", "15q25.2", "16p11.2", "16p13.3", "19p12a", "19p12b", "19p12c", "19p12d", "19p12e", "19p13.3", "19q11", "19q13.12b", "19q13.41", "19q13.42", "20q11.22", "22q11.21", "22q11.23", "Xq12", "Xq21.33", "Xq28a", "Xq28b", "Yp11.2", "Yq11.23a", "Yq11.23b")
for_VLP_analysis_mod = for_VLP_analysis %>%
  dplyr::slice(match(order,provirus))
rownames(for_VLP_analysis_mod)=for_VLP_analysis_mod$provirus
for_VLP_analysis_mod$provirus = NULL
for_VLP_analysis_final = for_VLP_analysis_mod[!row.names(for_VLP_analysis_mod)%in%remove,]


write.csv(for_VLP_analysis_final, "average individual HML-2 expression in GTEx, average >= 1.csv")

pdf("Average expression of individual HML-2 proviruses in GTEx, Average >= 1 row.pdf")
map = pheatmap(as.matrix(for_VLP_analysis_final), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
               main = "Average expression of individual HML-2 proviruses in GTEx, Average >= 1",breaks = seq(0,8,by=0.2),show_colnames = TRUE, fontsize = 8, fontsize_row = 6)
dev.off()
