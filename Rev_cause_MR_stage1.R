library(cause)
library(ieugwasr)
library(tidyverse)
library(data.table)
library(readr)
library(dplyr)
Disease1<-fread("/scratch/users/s/h/shifang/QTL/33888516-GCST90014033-EFO_0009552.h.tsv")
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
Disease2<-fread(exp_path,header=T)
Disease1<-Disease1[,c("hm_rsid","hm_beta","standard_error","hm_effect_allele","hm_other_allele","p_value")]
head(Disease1)
head(Disease2)

Disease1<-Disease1[,c("hm_rsid","hm_beta","standard_error","hm_effect_allele","hm_other_allele","p_value")]
head(Disease1)
head(Disease2)
Disease2_Disease1 <- gwas_merge(Disease1, Disease2, snp_name_cols = c("hm_rsid","SNP"),
                      beta_hat_cols = c("hm_beta","beta.outcome"),
                      se_cols = c("standard_error","se.outcome"),
                      A1_cols = c("hm_effect_allele","effect_allele.outcome"),
                      A2_cols = c("hm_other_allele","other_allele.outcome"))
## calculate nuisance parameters
# set random seed
set.seed(0)

Disease2_Disease1_varlist <- with(
  Disease2_Disease1,
  sample(snp, size=1000000, replace=FALSE)
)

Disease2_Disease1_params <- est_cause_params( Disease2_Disease1, Disease2_Disease1_varlist)
Disease20<-Disease1[Disease1$hm_rsid%in%Disease2_Disease1$snp,]
Disease20<-Disease20[,c("hm_rsid","p_value")]
colnames(Disease20)<-c("SNP","P")
write.table(Disease20, 'SNP.txt', quote = FALSE,row.names = FALSE)
save(Disease2_Disease1,Disease2_Disease1_params,file='cause_MR_stage1.Rdata')