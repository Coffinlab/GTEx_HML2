
#Variables
getwd()
setwd("/Users/aidanburn/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Gtex V8")
library(plyr)
library(stringi)
library(dplyr)
library(stringr)
library(readr)

#read telescope output in
Sample_Type="*_report"
File_Type="*-telescope_report.tsv"
Samplelist = list.files(pattern = Sample_Type)
for (i in 1:length(Samplelist)) {
filelist <- list.files(recursive = T, pattern = File_Type, full.names = T) #save data as a list
}
filelist

datalist = lapply(filelist, read.table, header=FALSE, sep ="\t", stringsAsFactors=FALSE, quote="\"")
datalist

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "-telescope_report.tsv"))
filelist_edited
filelist_edited = as.character(stri_sub(filelist_edited,22 ))
filelist_edited

names(datalist) = filelist_edited
head(datalist)


#Move each tissue sample to correct dataframe by match with sample file
V8ds <- read.delim("../GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
#Pull out SampleID and Tissue Site
V8ds_trim <- select(V8ds, SAMPID, SMTSD)
V8ds_trim$SMTSD <- gsub("\\s*\\([^\\)]+\\)", "", V8ds_trim$SMTSD)

#Run the first time working with dataset to generate unique tissue list 
for (i in (1:length(filelist_edited))) {
  c <-str_which(V8ds_trim$SAMPID, filelist_edited[i] )
  filelist_edited[i] = as.character(paste(V8ds_trim[c,2] , i, sep = "_"))
}


Tissue_list <- sort(filelist_edited)
Tissue_list_unique <- gsub("\\_.*","",Tissue_list)
Tissue_list_unique <- unique(Tissue_list_unique)
Tissue_list
Tissue_list_unique

write.csv(Tissue_list_unique, "Tissue_listref.csv")

#Create V8 list and loop through each tissue adding samples that we have to the V8 sample list for each bodysite
V8_samples <- list()
Tissue_list_unique <- read.delim("Tissue_listref.txt", header = FALSE)

  for (t in Tissue_list_unique$V1) {
    print(t)
    body_site <- V8ds_trim[grep(t, V8ds_trim$SMTSD),]
    body_site1 <- datalist[names(datalist) %in% body_site$SAMPID]
    V8_samples[[t]] <- body_site1
    
  }
 

#load in row names
gene_names = read.table("../all_ids", header=FALSE)
colnames(gene_names) = c("transcript")
head(gene_names)
dim(gene_names)

#Split V8_Samples to individual tissue file
for (i in (1:length(V8_samples))) {
  print(i)
  df <- V8_samples[[i]]
  #make empty data frame and append final counts data via matching to the first column
  Telescope_Counts = as.data.frame(gene_names)
  for (a in names(df)) {
    tempdf <- as.data.frame(df[a])
    colnames(tempdf) <- tempdf[1,]
    tempdf <- tempdf[-1,]
    Telescope_Counts <- join(Telescope_Counts, tempdf[,c(1,3)], by = "transcript", type = "left")
    colnames(Telescope_Counts)[colnames(Telescope_Counts) == 'final_count'] <- a
  }
  head(Telescope_Counts)
  Telescope_Counts[is.na(Telescope_Counts)] <- 0
  head(Telescope_Counts)
  rownames(Telescope_Counts) = Telescope_Counts$transcript
  Telescope_Counts[,1] = NULL
  head(Telescope_Counts)
  
  write.csv(Telescope_Counts, file = paste(names(V8_samples[i]), 'V8_Telescope_output.csv', sep="_"))
  
}  

#Rename V7 samples and Join with V8
#Load in the datalist of V7 samples
tissuemeta = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Full GTEx_HML2_Expression 3.nosync"
setwd(tissuemeta)
File_Type="*_Telescope_output.csv"
Telescopelist = list.files(tissuemeta, pattern = File_Type)
V7_datalist = lapply(Telescopelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)

Telescopelist_edited = as.character(strsplit(Telescopelist, "_Telescope_output.csv"))
Telescopelist_edited
names(V7_datalist) = Telescopelist_edited
head(V7_datalist)

#Replace V7 SRARunID with V8 GTEx Sample Name
for (i in names(V7_datalist)) {
  SRA_run_ID = read.delim(paste(tissuemeta,paste(i,"SraRunTable.txt",sep="_"),sep="/"),header=TRUE,sep="\t")
  meta_run_id <- SRA_run_ID[,c("Run","Sample_Name")]
  head(meta_run_id)
  names(V7_datalist[[i]]) <- meta_run_id$Sample_Name[match(names(V7_datalist[[i]]), meta_run_id$Run)]
  }
#Save new V7 files
for (i in names(V7_datalist)) {
  write.csv(V7_datalist[[i]], file = paste(names(V7_datalist[i]), "V7_SUBJID_Telescope_output.csv", sep = "_"))
}

#Join V7 and V8 files into V15 files
File_Type="*V15_Telescope_output.csv"
Telescopelist = list.files(pattern = File_Type)
for (i in Telescopelist) {
  tissue <- sub("_V15.*", "", i)
  print(tissue)
  
  #V15_Tissue <- read.csv(paste(tissue, "V15_Telescope_output.csv", sep = "_"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
  #V15_Tissue$X = NULL
  #names(V15_Tissue)[names(V15_Tissue)=="NA."] <- "X"
  V8_Tissue <- read.csv(paste(tissue, "V8_Telescope_output.csv", sep="_"), header=TRUE, sep =",", stringsAsFactors=FALSE)
  V7_Tissue <- read.csv(paste(tissue, "V7_SUBJID_Telescope_output.csv", sep="_"), header=TRUE, sep = ",", stringsAsFactors = FALSE)
  V7_Tissue$X = NULL
  names(V7_Tissue)[names(V7_Tissue)=="NA."] <- "X"
  V15_Tissue <- right_join(V7_Tissue, V8_Tissue, by = "X")
  write.csv(V15_Tissue, paste(tissue, "V15_Telescope_output.csv", sep="_"), row.names = FALSE)
}

