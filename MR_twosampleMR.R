##MR analysis using TwosampleMR

#Causal effect of risk factors on HEM
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table)
rm(list = ls())
path<-setwd("MR")
fileNames = list.files(path=path,pattern="*.csv", full.names = TRUE)
chd_out_dat_raw<-fread("33888516-GCST90014033-EFO_0009552.h.tsv",header=T)
chd_out_dat_raw<-chd_out_dat_raw[,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","hm_effect_allele_frequency","p_value","standard_error")] 
colnames(chd_out_dat_raw)=c("SNP","Chromozes","pos","other_allele.outcome","effect_allele.outcome","beta.outcome","eaf.outcome","pval.outcome","se.outcome")
chd_out_dat_raw$outcome<-c(rep("Ho",dim(chd_out_dat_raw)[1]))
chd_out_dat_raw$id.outcome<-c(rep("Ho",dim(chd_out_dat_raw)[1]))
df_total1 = data.frame()
df_total10 = data.frame()
for(i in 1:length(fileNames)){
  bmi_exp_dat<-fread(fileNames[i],header=T)
  mhc <- bmi_exp_dat %>%
    dplyr::arrange(chr.exposure, pos.exposure) %>%
    dplyr::filter(chr.exposure == 6) %>%
    dplyr::filter(pos.exposure >= 28477797 & pos.exposure <= 33448354)
  
  # create MHC region's snp list
  mhcsnp <- mhc %>% dplyr::select(SNP)
  
  # remove MHC's snp from the full GWAS
  bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
  
  bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","eaf.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","exposure","id.exposure")]
  chd_out_dat <- chd_out_dat_raw[chd_out_dat_raw$SNP%in% bmi_exp_dat $SNP,]
   dat <- harmonise_data(
      exposure_dat = bmi_exp_dat,
      outcome_dat = chd_out_dat
    )
    res <- mr(dat)
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
    het<-mr_heterogeneity(dat)
    het<-het[,c("method","Q","Q_df","Q_pval")]
    #Horizontal pleiotropy
    plt<-mr_pleiotropy_test(dat)##
    plt<-plt[,c("egger_intercept","se","pval")]
    res_single <- mr_singlesnp(dat)
    res_single_beta <- res_single$b
    res_single_se <- res_single$se
    res_Isq <- Isq(res_single_beta, res_single_se)
    test1<-subset(het,method=="MR Egger")
    test2<-subset(het,method=="Inverse variance weighted")
    df_total10<-cbind(res0,plt,test1,test2,res1,res2,res3,res4)
    df_total10$i_squared <- res_Isq
    df_total10$outcome<-c(rep("HEM",dim(df_total10)[1]))
    df_total10$expose<-c(rep(fileNames[i],dim(df_total10)[1]))
    df_total1 <- rbind(df_total1,df_total10)
  }
write.csv(df_total1,"risk_factors_to_HEM.csv", quote = FALSE,row.names = FALSE)

#Leave one out analysis
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table)
rm(list = ls())
path<-setwd("MR")
fileNames = list.files(path=path,pattern="*.csv", full.names = TRUE)
chd_out_dat_raw<-fread("33888516-GCST90014033-EFO_0009552.h.tsv",header=T)
chd_out_dat_raw<-chd_out_dat_raw[,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","hm_effect_allele_frequency","p_value","standard_error")] 
colnames(chd_out_dat_raw)=c("SNP","Chromozes","pos","other_allele.outcome","effect_allele.outcome","beta.outcome","eaf.outcome","pval.outcome","se.outcome")
chd_out_dat_raw$outcome<-c(rep("Ho",dim(chd_out_dat_raw)[1]))
chd_out_dat_raw$id.outcome<-c(rep("Ho",dim(chd_out_dat_raw)[1]))
df_total1 = data.frame()
df_total10 = data.frame()
for(i in 1:length(fileNames)){
  bmi_exp_dat<-fread(fileNames[i],header=T)
  bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","eaf.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","exposure","id.exposure")]
  chd_out_dat <- chd_out_dat_raw[chd_out_dat_raw$SNP%in% bmi_exp_dat $SNP,]
  dat <- harmonise_data(
    exposure_dat = bmi_exp_dat,
    outcome_dat = chd_out_dat
  )
  res <- mr(dat)
  res_loo <- mr_leaveoneout(dat)
  p1 <- mr_scatter_plot(res, dat)
  #p3 <- mr_leaveoneout_plot(res_loo)
  #p3[[1]]
  res_loo<-subset(res_loo,p>0.05)
  rownames(bmi_exp_dat)<-bmi_exp_dat$SNP
  remove_name<-res_loo$SNP
  bmi_exp_dat<-bmi_exp_dat[!(rownames(bmi_exp_dat) %in% remove_name),]
  dat <- harmonise_data(
    exposure_dat = bmi_exp_dat,
    outcome_dat = chd_out_dat)
  res <- mr(dat)
  if(dim(res)[1] == 0){
    print("None")}else if(res$nsnp[1]>=1){ab_odds <- generate_odds_ratios(res)
  res$ID<-c(rep(res[1,3],dim(res)[1]))
  res<-res[,c("method","nsnp","b","se","pval")]
  res$or<-ab_odds$or
  res$or_lci95<-ab_odds$or_lci95
  res$or_uci95<-ab_odds$or_uci95
  df_total10<-cbind(res)
  df_total10$outcome<-c(rep("HEM",dim(df_total10)[1]))
  df_total10$expose<-c(rep(fileNames[i],dim(df_total10)[1]))
  df_total1 <- rbind(df_total1,df_total10)
}}
write.csv(df_total1,"MR_LOCO_risks_to_HEM.csv", quote = FALSE,row.names = FALSE)

#Reverse causal effect of HEM on risk factors
rm(list = ls())
bmi_exp_dat<-fread("HO_clumped_exposure.csv",header=T)
df_total1 = data.frame()
df_total10 = data.frame()
test300<-read.table("123.txt",header=F)
#jj in c(39:42,44,46:50,52:55,57:60,63:66,68:72,74:86,88:91,93:99,102:104,108,110:115,117:121,123:128,131:137,139:144,146:148,150,153:160,163:164,166,168:182,184:195,197:199,201,202:204,206:212,214:216,219:234,236:238,240:248)
mhc <- bmi_exp_dat %>%
    dplyr::arrange(chr.exposure, pos.exposure) %>%
    dplyr::filter(chr.exposure == 6) %>%
    dplyr::filter(pos.exposure >= 28477797 & pos.exposure <= 33448354)
  
  # create MHC region's snp list
  mhcsnp <- mhc %>% dplyr::select(SNP)
  
  # remove MHC's snp from the full GWAS
  bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
  bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","samplesize.exposure","eaf.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","exposure")]
  
  for(j in 1:dim(test300)[1]){
    chd_out_dat <- extract_outcome_data(
      snps = bmi_exp_dat$SNP,
      outcomes =test300[j,2])
    dat <- harmonise_data(
      exposure_dat = bmi_exp_dat,
      outcome_dat = chd_out_dat)
    res <- mr(dat)
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
    het<-mr_heterogeneity(dat)
    het<-het[,c("method","Q","Q_df","Q_pval")]
    #Horizontal pleiotropy
    plt<-mr_pleiotropy_test(dat)##
    plt<-plt[,c("egger_intercept","se","pval")]
    res_single <- mr_singlesnp(dat)
    res_single_beta <- res_single$b
    res_single_se <- res_single$se
    res_Isq <- Isq(res_single_beta, res_single_se)
    test1<-subset(het,method=="MR Egger")
    test2<-subset(het,method=="Inverse variance weighted")
    df_total10<-cbind(res0,plt,test1,test2,res1,res2,res3,res4)
    df_total10$i_squared <- res_Isq
    df_total10$outcome<-c(rep("HEM",dim(df_total10)[1]))
    df_total10$expose<-c(rep(test300[j,2],dim(df_total10)[1]))
    df_total1 <- rbind(df_total1,df_total10)
  }
write.csv(df_total1,"HEM_to_risk_factors.csv", quote = FALSE,row.names = FALSE)


#Leave one out analysis
rm(list = ls())
bmi_exp_dat<-fread("HO_clumped_exposure.csv",header=T)
df_total1 = data.frame()
df_total10 = data.frame()
test300<-read.table("risk_factors_ID.txt",header=F)
mhc <- bmi_exp_dat %>%
  dplyr::arrange(chr.exposure, pos.exposure) %>%
  dplyr::filter(chr.exposure == 6) %>%
  dplyr::filter(pos.exposure >= 28477797 & pos.exposure <= 33448354)

# create MHC region's snp list
mhcsnp <- mhc %>% dplyr::select(SNP)

# remove MHC's snp from the full GWAS
bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","samplesize.exposure","eaf.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","exposure")]

for(j in 1:dim(test300)[1]){
  chd_out_dat <- extract_outcome_data(
    snps = bmi_exp_dat$SNP,
    outcomes =test300[j,2])
  dat <- harmonise_data(
    exposure_dat = bmi_exp_dat,
    outcome_dat = chd_out_dat)
  res <- mr(dat)
  res_loo <- mr_leaveoneout(dat)
  p1 <- mr_scatter_plot(res, dat)
  #p3 <- mr_leaveoneout_plot(res_loo)
  #p3[[1]]
  res_loo<-subset(res_loo,p>0.05)
  rownames(bmi_exp_dat)<-bmi_exp_dat$SNP
  remove_name<-res_loo$SNP
  bmi_exp_dat<-bmi_exp_dat[!(rownames(bmi_exp_dat) %in% remove_name),]
  dat <- harmonise_data(
    exposure_dat = bmi_exp_dat,
    outcome_dat = chd_out_dat)
  res <- mr(dat)
  if(dim(res)[1] == 0){
    print("None")}else if(res$nsnp[1]>=1){
    ab_odds <- generate_odds_ratios(res)
    res$ID<-c(rep(res[1,3],dim(res)[1]))
    res<-res[,c("method","nsnp","b","se","pval")]
    res$or<-ab_odds$or
    res$or_lci95<-ab_odds$or_lci95
    res$or_uci95<-ab_odds$or_uci95
    df_total10<-cbind(res)
    df_total10$outcome<-c(rep("HEM",dim(df_total10)[1]))
    df_total10$expose<-c(rep(test300[j,1],dim(df_total10)[1]))
    df_total1 <- rbind(df_total1,df_total10)
}}
write.csv(df_total1,"LOCO_HEM_to_risk_factors.csv", quote = FALSE,row.names = FALSE)


##Replication analysis using INTERVAL cohort

##extract cis-pQTL from INTERVAL cohort via ieu open gwas project (https://gwas.mrcieu.ac.uk/datasets/)

protein<-c("prot-a-979","prot-a-192")
data0<-fread("E:\\drug\\33888516-GCST90014033-EFO_0009552.h.tsv",header=T)
df_total1 = data.frame()
for(j in 1:2){
bmi_exp_dat <- extract_instruments(outcomes = protein[[j]])
mhc <- bmi_exp_dat %>%
  dplyr::arrange(chr.exposure, pos.exposure) %>%
  dplyr::filter(chr.exposure == 6) %>%
  dplyr::filter(pos.exposure >= 28477797 & pos.exposure <= 33448354)

# create MHC region's snp list
mhcsnp <- mhc %>% dplyr::select(SNP)
# remove MHC's snp from the full GWAS
bmi_exp_dat <- bmi_exp_dat %>% dplyr::filter(!(SNP %in% unlist(mhcsnp)))
bmi_exp_dat<-bmi_exp_dat[,c("SNP","beta.exposure","se.exposure","samplesize.exposure","eaf.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure","exposure")]
bmi_exp_dat <- clump_data(bmi_exp_dat)
dim(bmi_exp_dat)
chd_out_dat<-data0[,c("hm_rsid","hm_chrom","hm_pos","hm_other_allele","hm_effect_allele","hm_beta","hm_effect_allele_frequency","p_value","standard_error")] 
head (chd_out_dat)
colnames(chd_out_dat)=c("SNP","Chromozes","pos","other_allele.outcome","effect_allele.outcome","beta.outcome","eaf.outcome","pval.outcome","se.outcome")
chd_out_dat <- chd_out_dat[chd_out_dat$SNP%in% bmi_exp_dat $SNP,]
chd_out_dat$id.outcome<-c(rep("Ho",dim(chd_out_dat)[1]))
chd_out_dat$outcome<-c(rep("Ho",dim(chd_out_dat)[1]))
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat,
  outcome_dat = chd_out_dat
)
res <- mr(dat)
res
ab_odds <- generate_odds_ratios(res)
res$or<-ab_odds$or
res$or_lci95<-ab_odds$or_lci95
res$or_uci95<-ab_odds$or_uci95
res$outcome<-c(rep("HEM",dim(res)[1]))
res$expose<-c(rep(protein[[j]],dim(res)[1]))
df_total1 <- rbind(df_total1,res)
}
write.csv(df_total1,"replicated_MR_pQTL_HEM.csv", quote = FALSE,row.names = FALSE)