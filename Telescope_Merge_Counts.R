#Variables
getwd()
setwd("/Users/far122/Desktop/GTEx_10252019/kidney")
library(plyr)

#read telescope output in
File_Type="*-telescope_report.tsv"
filelist = list.files(pattern = File_Type) #save data as a list
filelist
datalist = lapply(filelist, read.table, header=TRUE, sep ="\t", stringsAsFactors=FALSE)
datalist

#load in row names
gene_names = read.table("../all_ids", header=FALSE)
colnames(gene_names) = c("transcript")
head(gene_names)
dim(gene_names)

#make empty data frame and append final counts data via matching to the first column
Telescope_Counts = as.data.frame(gene_names)

for (i in (1:length(datalist))) {
  print(i)
  df <- datalist[[i]]
  Telescope_Counts <- join(Telescope_Counts, df[,c(1,3)], by = "transcript", type = "left")
}  
  
#Telescope_Counts = cbind(merge(Telescope_Counts, datalist[[i]][,c("V1", "V3")], by.x = "Gene_names", by.y = "V1"))
 # add names to columns while binding
head(Telescope_Counts)
Telescope_Counts[is.na(Telescope_Counts)] <- 0
colnames(Telescope_Counts) <- c("Gene_names", unlist(lapply(strsplit(filelist, "-telescope_report.tsv"), function(x) x[[1]])) )
head(Telescope_Counts)
rownames(Telescope_Counts) = Telescope_Counts$Gene_names
Telescope_Counts[,1] = NULL
head(Telescope_Counts)

write.csv(Telescope_Counts, file = paste('8cell_1_Telescope_output', '.csv', sep="_"))
