##MR Rev analysis: Disease to pQTL
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(data.table)
setDTthreads(5)
#####################
df_total1 = data.frame()
df_total10 = data.frame()

##using pre-clumped GWASs datasets
args1 <- commandArgs(trailingOnly=TRUE)
exp_path1 <- args1[1]

bmi_exp_dat <-read.csv(exp_path1,header=T)
args2 <- commandArgs(trailingOnly=TRUE)
exp_path2 <- args2[2]
data<-fread(exp_path2,header=T)
chd_out_dat <- data[,c(4,1,2,6,5,7,1,8,10)]
rm(data)
colnames(chd_out_dat)=c("SNP","Chromozes","pos","other_allele.outcome","effect_allele.outcome","beta.outcome","eaf.outcome","pval.outcome","se.outcome")
chd_out_dat <- chd_out_dat[chd_out_dat$SNP%in% bmi_exp_dat $SNP,]
chd_out_dat$outcome<-c(rep("Ho",dim(chd_out_dat)[1]))
chd_out_dat$id.outcome<-c(rep("Ho",dim(chd_out_dat)[1]))
chd_out_dat<-chd_out_dat[,c("SNP","beta.outcome","se.outcome","eaf.outcome","pval.outcome","effect_allele.outcome","other_allele.outcome","outcome","Chromozes","id.outcome")]
dat <- harmonise_data(
  exposure_dat = bmi_exp_dat,
  outcome_dat = chd_out_dat
)
res <- mr(dat)
res100<-subset(res,method=="Inverse variance weighted")
res100$ID<-c(rep(exp_path2,dim(res100)[1]))
res$id.exposure<-c(rep(exp_path2,dim(res)[1]))
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
#het$outcome<-c(rep(ID0[i,1],dim(het)[1]))
#het$expose<-c(rep(ID0[jj,1],dim(het)[1]))
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
df_total10$Isq<-unIsq
df_total10$i_squared <- res_Isq
df_total1 <- rbind(df_total1,df_total10)
df_total1$ID<-c(rep(exp_path2,dim(df_total1)[1]))
df_total1$ID0<-c(rep(exp_path1,dim(df_total1)[1]))
write.table(df_total1,"mr.txt",quote = FALSE,row.names = FALSE)