##subset SNPs
library(data.table)
path<-setwd("/scratch/users/s//h/shifang/QTL/RAW")
fileNames = list.files(path=path,pattern="*.txt", full.names = TRUE)
for(i in 1:length(fileNames)){
  data<-fread(fileNames[i],header=T)
  data<-na.omit(data)
  fileNames[i] <- sub('/scratch/users/s//h/shifang/QTL/RAW/', '', fileNames[i])
  fileNames[i] <- gsub(" ", "", paste('/scratch/users/s/h/shifang/QTL/RAW/sig/clump/new/', fileNames[i]))
  data<-fread(fileNames[i],header=F)
  data<-data[data$V4 %in% data$SNP,]
  fileNames[i] <- sub('/scratch/users/s/h/shifang/QTL/RAW/sig/clump/new/', '', fileNames[i])
  fileNames[i] <- gsub(" ", "", paste('/scratch/users/s/h/shifang/QTL/RAW/sig/clump/new/final/', fileNames[i]))
  write.table(data,paste(fileNames[i],"txt",sep = "."), quote = FALSE,row.names = FALSE)
}