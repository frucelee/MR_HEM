library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table)
#####################
# Twosample MR  (proteins(cis-pQTL) to HEM)
#####################
df_total1 = data.frame()
df_total10 = data.frame()
chd_out_dat <- fread("/scratch/users/s/h/shifang/QTL/subset_SNPs_HEM.tsv")

args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]

#exposure
data<-fread(exp_path,header=T)
bmi_exp_dat <- data[,c(4,1,2,6,5,7,11,8,10,1)]
#rm(data)
colnames(bmi_exp_dat)=c("SNP","Chromozes","pos","other_allele.exposure","effect_allele.exposure","beta.exposure","samplesize.exposure","pval.exposure","se.exposure","eaf.exposure")
bmi_exp_dat$id.exposure<-c(rep("Ho",dim(bmi_exp_dat)[1]))
bmi_exp_dat$exposure<-c(rep("Ho",dim(bmi_exp_dat)[1]))
bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","samplesize.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","id.exposure","eaf.exposure","exposure")]
bmi_exp_dat<-na.omit(bmi_exp_dat)

data0 <- fread("/scratch/users/s/h/shifang/QTL/MR/cis_pQTL_meta.csv")
data0$chromosome_name <- sub("^", "chr", data0$chromosome_name )

data0 <-data0[data0$gene2%in% exp_path,]
cc<-data0[1,3]
data<-subset(data,V1==cc[[1]])
#data3$V6<-10^((-1)*data3$V6)
dd<-data0[1,4]
ee<-data0[1,5]
data<-subset(data,V2>=(dd[[1]]-500000)&V2<=(ee[[1]]+500000))  ##cis-PQTL (±500 kb of gene body)
bmi_exp_dat<-bmi_exp_dat[bmi_exp_dat$SNP%in% data$V4,]

# read the outcome
chd_out_dat <- chd_out_dat[chd_out_dat$SNP%in% bmi_exp_dat $SNP,]
#chd_out_dat<-chd_out_dat[,-8]
dat <- harmonise_data(bmi_exp_dat,chd_out_dat)
res <- mr(dat)
if(dim(res)[1] == 0){
  print("None")
}else if(res$nsnp[1]==1){
    res1<-subset(res,method=="Wald ratio")
    res1$ID<-c(rep(exp_path,dim(res1)[1]))
    dat$samplesize.exposure<-c(rep("35559",dim(dat)[1]))
    dat$samplesize.outcome<-c(rep("944133",dim(dat)[1]))
    dat$samplesize.exposure<-as.numeric(dat$samplesize.exposure)
    dat$samplesize.outcome<-as.numeric(dat$samplesize.outcome)
    p<-directionality_test(dat)
    
    ab_odds <- generate_odds_ratios(res)
    ab_odds<-subset(ab_odds,method=="Wald ratio")
    res1$or<-ab_odds$or
    res1$or_lci95<-ab_odds$or_lci95
    res1$or_uci95<-ab_odds$or_uci95
    res1$rev<-p$correct_causal_direction
    res1$rev_p<-p$steiger_pval
    output_dir="pQTL_to_HEM"
    #exp_path <- gsub(" ", "", paste(exp_path, '.txt'))
    exp_path <- gsub('.txt','',exp_path)
    output_dir <- gsub(" ", "", paste(output_dir, exp_path))
    write.table(res1,paste(output_dir,"WR",sep = "."), quote = FALSE,row.names = FALSE)
  }else if(res$nsnp[1]==2){
    dat$samplesize.exposure<-c(rep("35559",dim(dat)[1]))
    dat$samplesize.outcome<-c(rep("944133",dim(dat)[1]))
    dat$samplesize.exposure<-as.numeric(dat$samplesize.exposure)
    dat$samplesize.outcome<-as.numeric(dat$samplesize.outcome)
    p<-directionality_test(dat)
    ab_odds <- generate_odds_ratios(res)
    ab_odds<-subset(ab_odds,method=="Inverse variance weighted")
    res1<-subset(res,method=="Inverse variance weighted")
    res1$or<-ab_odds$or
    res1$or_lci95<-ab_odds$or_lci95
    res1$or_uci95<-ab_odds$or_uci95
    het<-mr_heterogeneity(dat)
    het<-het[,c("method","Q","Q_df","Q_pval")]
    test2<-subset(het,method=="Inverse variance weighted")
    df_total10<-cbind(res1,test2)
    res_single <- mr_singlesnp(dat)
    res_single_beta <- res_single$b
    res_single_se <- res_single$se
    res_Isq <- Isq(res_single_beta, res_single_se)
    df_total10$i_squared <- res_Isq
    df_total10$ID<-c(rep(exp_path,dim(df_total10)[1]))
    df_total10$rev<-p$correct_causal_direction
    df_total10$rev_p<-p$steiger_pval
    df_total10<-subset(df_total10,nsnp<=2)
    output_dir="/scratch/users/s/h/shifang/QTL/MR/pQTL_to_HEM/"
    #exp_path <- gsub(" ", "", paste(exp_path, '.txt'))
    exp_path <- gsub('.txt','',exp_path)
    output_dir <- gsub(" ", "", paste(output_dir, exp_path))
    write.table(df_total10,paste(output_dir,"IVW",sep = "."), quote = FALSE,row.names = FALSE)
    
  }else if(res$nsnp[1]>=3) {
    dat$samplesize.exposure<-c(rep("35559",dim(dat)[1]))
    dat$samplesize.outcome<-c(rep("944133",dim(dat)[1]))
    dat$samplesize.exposure<-as.numeric(dat$samplesize.exposure)
    dat$samplesize.outcome<-as.numeric(dat$samplesize.outcome)
    p<-directionality_test(dat)
    
    res100<-subset(res,method=="Inverse variance weighted")
    res100$ID<-c(rep(exp_path,dim(res100)[1]))
    res$id.exposure<-c(rep(exp_path,dim(res)[1]))
    output_dir="pQTL_to_HEM"
    #exp_path <- gsub(" ", "", paste(exp_path, '.txt'))
    exp_path <- gsub('.txt','',exp_path)
    output_dir <- gsub(" ", "", paste(output_dir, exp_path))
    ab_odds <- generate_odds_ratios(res)
    res$ID<-c(rep(res[1,3],dim(res)[1]))
    res<-res[,c("method","nsnp","b","se","pval")]
    res$or<-ab_odds$or
    res$or_lci95<-ab_odds$or_lci95
    res$or_uci95<-ab_odds$or_uci95
    res0<-subset(res,method=="Inverse variance weighted")
    res1<-subset(res,method=="MR Egger")
    res1<-res1[,c(1,3:5)]
    res2<-subset(res,method=="Weighted median")
    res2<-res2[,c(1,3:5)]
    res3<-subset(res,method=="Simple mode")
    res3<-res3[,c(1,3:5)]
    res4<-subset(res,method=="Weighted mode")
    res4<-res4[,c(1,3:5)]
    #res0<-subset(res0,pval<=0.05)
    #Draw scatter plot
    #p1 <- mr_scatter_plot(res, dat)
    #p1[[1]]
    #other analysis
    #Heterogeneity statistics
    het<-mr_heterogeneity(dat)
    het<-het[,c("method","Q","Q_df","Q_pval")]
   #Horizontal pleiotropy
    plt<-mr_pleiotropy_test(dat)##
    plt<-plt[,c("egger_intercept","se","pval")]
    # Calculate F statistics
    # and I-squared statistics
    # to measure Instrument
    # strength for MR-Egger
    F= BXG^2/seBetaXG^2
    plt$F.statics= mean(F)
    res_single <- mr_singlesnp(dat)
    res_single_beta <- res_single$b
    res_single_se <- res_single$se
    res_Isq <- Isq(res_single_beta, res_single_se)
    test1<-subset(het,method=="MR Egger")
    test2<-subset(het,method=="Inverse variance weighted")
    df_total10<-cbind(res0,plt,test1,test2,res1,res2,res3,res4)
    df_total10$i_squared <- res_Isq
    df_total1 <- rbind(df_total1,df_total10)
    df_total1$ID<-c(rep(exp_path,dim(df_total1)[1]))
    df_total1$rev<-p$correct_causal_direction
    df_total1$rev_p<-p$steiger_pval
    write.table(df_total1,paste(output_dir,"txt",sep = "."), quote = FALSE,row.names = FALSE)
  }