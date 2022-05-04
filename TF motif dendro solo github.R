library(readr)
library(tidyverse)
library(dplyr)
library(pvclust)

setwd("/System/Volumes/Data/Users/aidanburn/Library/CloudStorage/Box-Box/Will Rotation/Will's Rotation Files/PANDA/PANDA_Data")


HML2_MotifScan <- read_tsv("~/Documents/Coffin Lab/fimo_out/fimo.tsv")
colnames(HML2_MotifScan) <- gsub(x = colnames(HML2_MotifScan), "\\-", "_")
updatefimo <- read_tsv(list.files("../PANDA_Data/Motif Data/HervHumanTF_fimo_out/12q14.1", pattern = "fimo.tsv", full.names = T))
colnames(updatefimo) <- gsub(x = colnames(updatefimo), "\\-", "_")
HML2_MotifScan <- rbind(HML2_MotifScan,updatefimo)
updatefimo2 <- read_tsv(list.files("../PANDA_Data/Motif Data/HervHumanTF_fimo_out/22q11.23HS", pattern = "fimo22hs.tsv", full.names = T))
colnames(updatefimo2) <- gsub(x = colnames(updatefimo2), "\\-", "_")
HML2_MotifScan <- rbind(HML2_MotifScan,updatefimo2)
table(HML2_MotifScan$sequence_name)
SingleMotif <- HML2_MotifScan%>% filter(sequence_name == "20q11.22-5-LTR")
SingleMotiftrim <- SingleMotif %>% filter(start < 250)
SingleMotiftrimafter <- SingleMotif %>% filter(start > 842)
SingleMotifedit <- rbind(SingleMotiftrim,SingleMotiftrimafter)
SingleMotif <- HML2_MotifScan%>% filter(sequence_name != "20q11.22-5-LTR")
SingleMotif <- rbind(SingleMotif, SingleMotifedit)
HML2_MotifScan = SingleMotif
HML2_MotifScan <- HML2_MotifScan %>%
  filter(q_value < 0.05)
binMatrix_MotifScan <- select(HML2_MotifScan, motif_alt_id, sequence_name)
binMatrix_MotifScan <- table(binMatrix_MotifScan)
binMatrix_MotifScan <- as.data.frame.matrix(binMatrix_MotifScan)
binMatrix_MotifScan <- as.matrix(t(binMatrix_MotifScan))
binMatrix_MotifScan

binMatrix_MotifScan_df <- as.data.frame(binMatrix_MotifScan)
binMatrix_MotifScan <- binMatrix_MotifScan_df[-c(30,33,44,45,56,57,62),]
t_binMatrix_MotifScan <- t(binMatrix_MotifScan)
ltr.pv <-pvclust(t_binMatrix_MotifScan, nboot = 1000, method.hclust="ward.D2")
pdf("HML2_LTRMotif_Cluster_2LTR_boot_procut_12qadd22add.pdf", width = 15)
plot(ltr.pv, print.pv =1, print.num=FALSE, col.pv=c(si=1, au=1, bp=1))
dev.off()

                                               