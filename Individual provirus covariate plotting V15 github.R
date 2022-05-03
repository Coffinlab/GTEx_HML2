dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
raw_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
TPM_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures"
metadata = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
tissuemeta = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Full GTEx_HML2_Expression 3.nosync"
PCA=paste(figures,"PCA",sep="/")
#PCA_data = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/PCA"
Demographics = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
setwd(TPM_counts)

library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape)


#Limma for Individual Sex
#Read in raw Counts for Breast
setwd("../V15 Telescope")
library(edgeR)
counts <-read.csv("Breast_Mammary_Tissue_V15_Telescope_output.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(counts)
#Clean up Row/Column names and remove unnecessary columns
rownames(counts) <- counts[,1]
counts <- counts[,-c(1,69,240,295,321,194,212)]
counts <- counts[-c(26086:26087),]
colnames(counts) <- gsub("\\.", "-",colnames(counts))

#Remove Y chromosome outliers
colnum <- which(colnames(counts) == "GTEX-11ILO-2226-SM-5A5L1")
counts <- counts[,-colnum]

# #Generate Meta file
 meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
 meta_mod = meta[,c("SUBJID","AGE", "SEX")]
 colnames(meta_mod) = c("SUBJID","age", "sex")
 head(meta_mod)


 GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
 GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
 GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
 GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
 GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
 GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == "Breast_Mammary_Tissue")
 string = GTEx_sample_ID_tissue$SAMPID
 subjid = sub("^([^.]+.[^.]+).*", "\\1", string)
 subjid = gsub("\\.", "-", string)
#colnames(df) = subjid
subjid = sub("^(.*?-.*?)-.*", "\\1", string)
head(subjid)
length(subjid)
# 385
GTEx_sample_ID_tissue$SUBJID = subjid
#SRA_run_ID_to_merge$SUBJID = subjid
meta_mod$SUBJID
meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
head(meta_updated)
dim(meta_updated)
meta_updated$sex[meta_updated$sex == 1] <- "Male"
meta_updated$sex[meta_updated$sex == 2] <- "Female"
meta_updated$squish <- paste(meta_updated$SAMPID, meta_updated$sex, sep = "_")

colnames(counts) <- meta_updated$squish[match(names(counts),meta_updated$SAMPID)]


d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

cutoff <- 20
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

group <- sub(".*_", "", colnames(counts))
group
df=data.frame(group)
df$col = ifelse(group == "Male","red","blue")
df

plotMDS(d, col = c(df$col), top = 500)

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

#Redo MDS and voom with just HERVs
grep("HML-2", rownames(d$counts))

fit <- lmFit(y, mm)
head(coef(fit))


contr <- makeContrasts(groupFemale - groupMale, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])


#Individual Sex

HML2_Tissue = "Breast_Mammary_Tissue"
HML2_Tissue
HML2_df <- read.csv(paste(HML2_Tissue, "TPM_Log2_HML2.csv", sep = "_"),header=TRUE, sep =",", stringsAsFactors=FALSE)
head(HML2_df)[1:5,1:5]
HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
HML2_DF$X
HML2_DF$difference = NULL
head(HML2_DF)
colnames(HML2_DF)
colnames(HML2_DF) <- gsub("\\.", "-",colnames(HML2_DF))
colnames(HML2_DF)
ERV_Counts_modT = as.data.frame(t(HML2_DF))
#ERV_Counts_modT$tissue = HML2_Tissue
ERV_Counts_modT$SAMPID <- rownames(ERV_Counts_modT)
head(ERV_Counts_modT)
colnames(ERV_Counts_modT) <- ERV_Counts_modT["X",]
ERV_Counts_modT <- dplyr::rename(ERV_Counts_modT, SAMPID = "X")
ERV_Counts_modT <- ERV_Counts_modT[-1,]
head(ERV_Counts_modT)

meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX")]
colnames(meta_mod) = c("SUBJID","age", "sex")
head(meta_mod)

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

ERV_Counts_modT_pheno = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
head(ERV_Counts_modT_pheno)

ERV_Counts_modT_trim_sex = subset(ERV_Counts_modT_pheno, select=-c(SMTSD,SAMPID,age))
head(ERV_Counts_modT_trim_sex)

ERV_Counts_modT_trim_sex_filt <- filter_all(ERV_Counts_modT_trim_sex, any_vars(. > 1))

ERV_Counts_modT_melt_sex = reshape2::melt(ERV_Counts_modT_trim_sex, c("SUBJID", "sex"))


ERV_Counts_modT_melt_sex$value <- as.numeric(ERV_Counts_modT_melt_sex$value)
ERV_Counts_modT_melt_sex$sex <- as.factor(ERV_Counts_modT_melt_sex$sex)
head(ERV_Counts_modT_melt_sex)

ERV_Counts_modT_melt_sex_filt <- ERV_Counts_modT_melt_sex %>% dplyr::filter(value>1)

levels(ERV_Counts_modT_melt_sex_filt$sex) <- c(levels(ERV_Counts_modT_melt_sex_filt$sex), "Male")
ERV_Counts_modT_melt_sex_filt$sex[ERV_Counts_modT_melt_sex_filt$sex == 1] <- "Male"
levels(ERV_Counts_modT_melt_sex_filt$sex) <- c(levels(ERV_Counts_modT_melt_sex_filt$sex), "Female")
ERV_Counts_modT_melt_sex_filt$sex[ERV_Counts_modT_melt_sex_filt$sex == 2] <- "Female"

ERV_Counts_modT_melt_sex_filt <- na.omit(ERV_Counts_modT_melt_sex_filt)



trial <- ERV_Counts_modT_melt_sex_filt[ERV_Counts_modT_melt_sex_filt$variable %in% names(which(table(ERV_Counts_modT_melt_sex_filt$variable) > 5)), ]


my_comparisons <- list(c("Male", "Female"))
eliminated_sex$sex <- factor(eliminated_sex$sex, levels =c("Male", "Female"))


pdf("breast_boxout.pdf", width = 8)
myplot <- ggplot(trial, aes(x=variable, y=value, fill = sex)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Biological Sex in Breast", x="Provirus", y = "TPM")
print(myplot) 
dev.off()

#Limma for Individual Age
#Read in raw Counts for Spinal Cord
setwd("../V15 Telescope")
library(edgeR)
counts <-read.csv("Brain_Spinal_Cord_V15_Telescope_output.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(counts)
#Clean up Row/Column names and remove unnecessary columns
rownames(counts) <- counts[,1]
counts <- counts[,-(1)]
counts <- counts[-c(26086:26087),]
colnames(counts) <- gsub("\\.", "-",colnames(counts))


# #Generate Meta file
meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX")]
colnames(meta_mod) = c("SUBJID","age", "sex")
head(meta_mod)


GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == "Brain_Spinal_cord")
string = GTEx_sample_ID_tissue$SAMPID
subjid = sub("^([^.]+.[^.]+).*", "\\1", string)
subjid = gsub("\\.", "-", string)
#colnames(df) = subjid
subjid = sub("^(.*?-.*?)-.*", "\\1", string)
head(subjid)
length(subjid)
# 385
GTEx_sample_ID_tissue$SUBJID = subjid
#SRA_run_ID_to_merge$SUBJID = subjid
meta_mod$SUBJID
meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
head(meta_updated)
dim(meta_updated)


meta_updated_trim = subset(meta_updated, select=-c(SMTSD,SUBJID,sex))
head(meta_updated_trim)


young = c(20:35)
middle = c(36:51)
older = c(52:70)
meta_updated_mute <- meta_updated_trim %>%
  mutate(age = case_when(age %in% young ~ "20-35",
                         age %in% middle ~ "36-51",
                         age %in% older ~ "52-70"))

meta_updated_mute$age[meta_updated_mute$age == "20-35"] <- "Low"
meta_updated_mute$age[meta_updated_mute$age == "36-51"] <- "Middle"
meta_updated_mute$age[meta_updated_mute$age == "52-70"] <- "High"
meta_updated_mute$squish <- paste(meta_updated_mute$SAMPID, meta_updated_mute$age, sep = "_")


colnames(counts) <- meta_updated_mute$squish[match(names(counts),meta_updated$SAMPID)]



d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

group <- sub(".*_", "", colnames(counts))
group
df=data.frame(group)
df$col <- ifelse(df$group != "Low" & df$group != "Middle","red",
                 ifelse(df$group != "Low" & df$group != "High","blue","green"))

plotMDS(d, col = df$col, top = 500)

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))

#low and high contrasts
contr <- makeContrasts(groupLow - groupHigh, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "spinalcordageretoplowhigh.csv")

#low and middle contrasts
contr <- makeContrasts(groupLow - groupMiddle, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "spinalcordageretoplowmid.csv")

#middle and high contrast
contr <- makeContrasts(groupMiddle - groupHigh, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "spinalcordageretopmidhigh.csv")

#Individual Age 

HML2_Tissue = "Brain_Spinal_Cord"
HML2_Tissue
HML2_df <- read.csv(paste(HML2_Tissue, "TPM_HML2.csv", sep = "_"),header=TRUE, sep =",", stringsAsFactors=FALSE)
head(HML2_df)[1:5,1:5]
HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
HML2_DF$X
HML2_DF$difference = NULL
head(HML2_DF)
colnames(HML2_DF)
colnames(HML2_DF) <- gsub("\\.", "-",colnames(HML2_DF))
ERV_Counts_modT = as.data.frame(t(HML2_DF))
#ERV_Counts_modT$tissue = HML2_Tissue
ERV_Counts_modT$SAMPID <- rownames(ERV_Counts_modT)
head(ERV_Counts_modT)
colnames(ERV_Counts_modT) <- ERV_Counts_modT["X",]
ERV_Counts_modT <- ERV_Counts_modT[-1,]
ERV_Counts_modT <- dplyr::rename(ERV_Counts_modT, "SAMPID" = "X")
head(ERV_Counts_modT)

meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX")]
colnames(meta_mod) = c("SUBJID","age", "sex")
head(meta_mod)

GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == "Brain_Spinal_cord")
string = GTEx_sample_ID_tissue$SAMPID

subjid = sub("^(.*?-.*?)-.*", "\\1", string)
head(subjid)
length(subjid)

GTEx_sample_ID_tissue$SUBJID = subjid

meta_mod$SUBJID
meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
head(meta_updated)
dim(meta_updated)

ERV_Counts_modT_pheno = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
head(ERV_Counts_modT_pheno)

ERV_Counts_modT_trim_age = subset(ERV_Counts_modT_pheno, select=-c(SMTSD,SAMPID,sex))
head(ERV_Counts_modT_trim_age)

ERV_Counts_modT_melt_age = reshape2::melt(ERV_Counts_modT_trim_age, c("SUBJID", "age"))


young = c(20:35)
middle = c(36:51)
older = c(52:70)
ERV_Counts_modT_melt_age <- ERV_Counts_modT_melt_age %>%
  mutate(age = case_when(age %in% young ~ "20-35",
                         age %in% middle ~ "36-51",
                         age %in% older ~ "52-70"))


ERV_Counts_modT_melt_age$value <- as.numeric(ERV_Counts_modT_melt_age$value)

ERV_Counts_modT_melt_age_filt <- ERV_Counts_modT_melt_age %>% dplyr::filter(value>2)
ERV_Counts_modT_melt_age_filt <- na.omit(ERV_Counts_modT_melt_age_filt)

trial <- ERV_Counts_modT_melt_age_filt[ERV_Counts_modT_melt_age_filt$variable %in% names(which(table(ERV_Counts_modT_melt_age_filt$variable) > 5)), ]



pdf("spinal_cord_box_age.pdf", width = 8)
myplot <- ggplot(trial, aes(x=variable, y=value, fill = age)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Age in Spinal Cord",x="Provirus", y = "TPM")
print(myplot) 
dev.off()

stat_compare_means(method = "t.test", aes(label = ..p.signif..)) 

#Limma for Individual Hardy Score
#Read in raw Counts for Cerebellum
setwd("../V15 Telescope")
library(edgeR)
counts <-read.csv("Brain_Cerebellum_V15_Telescope_output.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(counts)
#Clean up Row/Column names and remove unnecessary columns
rownames(counts) <- counts[,1]
counts <- counts[,-(1)]
counts <- counts[,-c(7,27,109,111,113,116,136,143,149,157,184,185,194)]
counts <- counts[,-c(121,123,125)]
counts <- counts[-c(26086:26087),]
colnames(counts) <- gsub("\\.", "-",colnames(counts))


# #Generate Meta file
meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX","DTHHRDY")]
colnames(meta_mod) = c("SUBJID","age", "sex","hardy")
head(meta_mod)


GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == "Brain_Cerebellum")
string = GTEx_sample_ID_tissue$SAMPID
subjid = sub("^([^.]+.[^.]+).*", "\\1", string)
subjid = gsub("\\.", "-", string)

subjid = sub("^(.*?-.*?)-.*", "\\1", string)
head(subjid)
length(subjid)

GTEx_sample_ID_tissue$SUBJID = subjid

meta_mod$SUBJID
meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
head(meta_updated)
dim(meta_updated)


meta_updated_trim = subset(meta_updated, select=-c(SMTSD,SUBJID,age,sex))
head(meta_updated_trim)

meta_updated_trim <- subset(meta_updated_trim, hardy != 0)



meta_updated_trim$squish <- paste(meta_updated_trim$SAMPID, meta_updated_trim$hardy, sep = "_")


colnames(counts) <- meta_updated_trim$squish[match(names(counts),meta_updated_trim$SAMPID)]



d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)

group <- sub(".*_", "", colnames(counts))
group
df=data.frame(group)
df$col <- ifelse(df$group != "1" & df$group != "2" & df$group != "3","red",
                 ifelse(df$group != "1" & df$group != "2" & df$group != "4","blue",
                        ifelse(df$group !=1 & df$group != "3" & df$group != "4","green","yellow")))

plotMDS(d, col = df$col, top =500, gene.selection = "common")



mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

fit <- lmFit(y, mm)
head(coef(fit))



#1 and 2 contrasts
contr <- makeContrasts(group1 - group2, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "cerebellumhardyretop1-2.csv")

#1 and 4 contrast
contr <- makeContrasts(group1 - group4, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "cerebellumhardyretop1-4.csv")

#2 and 4 contrast
contr <- makeContrasts(group2 - group4, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "cerebellumhardyretop2-4.csv")

#2 and 3 contrast
contr <- makeContrasts(group2 - group3, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "cerebellumhardyretop2-3.csv")

#3 and 4 contrast
contr <- makeContrasts(group3 - group4, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
hmlindic <- grep("HML-2", rownames(tmp))
top.tablehml <- topTable(tmp[hmlindic,])
write.csv(top.tablehml, "cerebellumhardyretop3-4.csv")



#Individual Hardy Score



HML2_Tissue = "Brain_Cerebellum"
HML2_Tissue
HML2_df <- read.csv(paste(HML2_Tissue, "TPM_HML2.csv", sep = "_"),header=TRUE, sep =",", stringsAsFactors=FALSE)
head(HML2_df)[1:5,1:5]
HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
HML2_DF$X
HML2_DF$difference = NULL
head(HML2_DF)
colnames(HML2_DF)
colnames(HML2_DF) <- gsub("\\.", "-",colnames(HML2_DF))
ERV_Counts_modT = as.data.frame(t(HML2_DF))
ERV_Counts_modT$SAMPID <- rownames(ERV_Counts_modT)
head(ERV_Counts_modT)
colnames(ERV_Counts_modT) <- ERV_Counts_modT["X",]
ERV_Counts_modT <- ERV_Counts_modT[-1,]
ERV_Counts_modT <- dplyr::rename(ERV_Counts_modT, "SAMPID" = "X")
head(ERV_Counts_modT)

meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX", "DTHHRDY")]
colnames(meta_mod) = c("SUBJID","age", "sex","hardy")
head(meta_mod)

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

ERV_Counts_modT_pheno = join(ERV_Counts_modT, meta_updated, by = "SAMPID", type = "left")
head(ERV_Counts_modT_pheno)

ERV_Counts_modT_trim_hardy = subset(ERV_Counts_modT_pheno, select=-c(SMTSD,SAMPID,age,sex))
head(ERV_Counts_modT_trim_hardy)

ERV_Counts_modT_melt_hardy = reshape2::melt(ERV_Counts_modT_trim_hardy, c("SUBJID", "hardy"))

ERV_Counts_modT_melt_hardy_filt <- ERV_Counts_modT_melt_hardy %>% dplyr::filter(value>2)
ERV_Counts_modT_melt_hardy_filt <- na.omit(ERV_Counts_modT_melt_hardy_filt)



ERV_Counts_modT_melt_hardy_filt$hardy <- as.factor(ERV_Counts_modT_melt_hardy_filt$hardy)
ERV_Counts_modT_melt_hardy_filt$value <- as.numeric(ERV_Counts_modT_melt_hardy_filt$value)

trial <- ERV_Counts_modT_melt_hardy_filt[ERV_Counts_modT_melt_hardy_filt$variable %in% names(which(table(ERV_Counts_modT_melt_hardy_filt$variable) > 5)), ]

table(trial$variable)

pdf("Cerebellum_hardy.pdf", width = 8)
myplot <- ggplot(trial, aes(x=variable, y=value, fill = hardy)) +
  geom_boxplot(position=position_dodge(width=0.75)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Hardy Score in Cerebellum",x="Provirus", y = "TPM")
print(myplot) 
dev.off()

png("Cerebellum_box_age.png", width = 1000)
myplot <- ggplot(ERV_Counts_modT_melt_age_filt, aes(x=variable, y=value, color = age)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title= "Age in Spinal Cord",x="Provirus", y = "TPM")
print(myplot) 
dev.off()
