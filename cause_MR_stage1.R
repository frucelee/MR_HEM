library(cause)
library(ieugwasr)
library(tidyverse)
library(data.table)
library(readr)
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
exp_path <- args[1]
Disease1<-fread("/scratch/users/s/h/shifang/QTL/33888516-GCST90014033-EFO_0009552.h.tsv")
Disease2<-fread(exp_path,header=T)
Disease1<-Disease1[,c("hm_rsid","hm_beta","standard_error","hm_effect_allele","hm_other_allele","p_value")]
head(Disease1)
head(Disease2)
Disease2_Disease1 <- gwas_merge(Disease2, Disease1, snp_name_cols = c("SNP","hm_rsid"),
                      beta_hat_cols = c("beta.outcome","hm_beta"),
                      se_cols = c("se.outcome","standard_error"),
                      A1_cols = c("effect_allele.outcome","hm_effect_allele"),
                      A2_cols = c("other_allele.outcome","hm_other_allele"))
## calculate nuisance parameters
# set random seed
set.seed(0)

Disease2_Disease1_varlist <- with(
  Disease2_Disease1,
  sample(snp, size=1000000, replace=FALSE)
)

Disease2_Disease1_params <- est_cause_params( Disease2_Disease1, Disease2_Disease1_varlist)
Disease20<-Disease2[Disease2$SNP%in%Disease2_Disease1$snp,]
Disease20<-Disease20[,c("SNP","pval.outcome")]
colnames(Disease20)<-c("SNP","P")
write.table(Disease20, 'SNP.txt', quote = FALSE,row.names = FALSE)
save(Disease2_Disease1,Disease2_Disease1_params,file='cause_MR_stage1.Rdata')