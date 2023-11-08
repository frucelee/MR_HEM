##Colocalization analysis
library(data.table)
library(magrittr) 
library(dplyr)  
library(coloc)
HEM_gwas<-fread("33888516-GCST90014033-EFO_0009552.h.tsv",header=T)
path<-setwd("coloc")
fileNames = list.files(path=path,pattern=".txt", full.names = TRUE)
coloc_sum <- matrix(nrow=175,
                     ncol = 7)
coloc_sum<-data.frame(coloc_sum)
for (j in c(1:175)){
  pQTL_subset<-fread(fileNames[j],header=F)
  head(pQTL_subset)
  data<-pQTL_subset
  pQTL_subset$V1 <- sub('chr', '', pQTL_subset$V1)
  data %<>% dplyr::select(V2, V7, V10, V5,V6, V12, V11, V8, V4)
  colnames(data)=c("position", "beta", "SE", "ALT", "REF", "MAF", "N", "pvalues", "rsid")
  #data6<-data[,c(9,8)]
  # give varbeta
  data$varbeta <- (data$SE)^2
  # limit rsid to one
  data$rsid <- unlist(lapply(strsplit(data$rsid, ","), function(x) x[1]))
  data<-na.omit(data)
  data=data[!duplicated(data$rsid),]
  #
  data3<-subset(HEM_gwas,Chromozes==pQTL_subset$V1[[1]])
  data3<-subset(data3,pos>=min(pQTL_subset$V2)&pos<=max(pQTL_subset$V2))
  data3 %<>% dplyr::select("hm_chrom", "hm_pos", "hm_beta", "standard_error", "hm_effect_allele", "hm_other_allele", "hm_effect_allele_frequency", "p_value", "hm_rsid") 
  colnames(data3)=c("chr", "position", "beta", "SE", "ALT", "REF", "MAF", "pvalues", "rsid")
  #data7<-data3[,c(9,8)]
  # give varbeta
  data3$varbeta <- (data3$SE)^2
  # limit rsid to one
  data3$rsid <- unlist(lapply(strsplit(data3$rsid, ","), function(x) x[1]))
  data3=data3[!duplicated(data3$rsid),]
  #devtools::install_github("chr1swallace/coloc",force = TRUE)
  input <- merge(data, data3, by="rsid", all=FALSE, suffixes=c("_eqtl","_gwas"))
  library("coloc")
  result <- coloc.abf(dataset1=list(pvalues=input$pvalues_eqtl, snp=input$rsid,type="quant", N=35559,MAF=input$MAF_eqtl), dataset2=list(pvalues=input$pvalues_gwas,MAF=input$MAF_gwas, snp=input$rsid, type="cc", s=0.23, N=944133))
  library(dplyr)
  dd<-data.frame(t(data.frame(print(result[[1]]))))
  dd$ID<-fileNames[j]
  coloc_sum[j,]<-dd[1, ]
}
write.csv(coloc_sum,"HEM_coloc.csv", quote = FALSE,row.names = FALSE)

##you can view the results of colocalization analysis using locuscomparer based on above
#library(locuscomparer)
#colnames(data7)=c("rsid","pval")
#colnames(data6)=c("rsid","pval")
#data6$pval<-as.numeric(data6$pval)
#data7$pval<-as.numeric(data7$pval)
#locuscompare(in_fn1 = data6, in_fn2 = data7, title1 = 'pQTL', title2 = 'GWAS',snp="rs5242580")
