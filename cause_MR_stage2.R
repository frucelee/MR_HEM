library(cause)
library(ieugwasr)
library(tidyverse)
library(data.table)
library(readr)
library(dplyr)
load ('cause_MR_stage1.Rdata')
Disease2_Disease1_clump<-fread("HO_clumped.txt")

# fit CAUSE
Disease2_Disease1_res <- cause(
  X=Disease2_Disease1,
  variants = Disease2_Disease1_clump$SNP,
  param_ests = Disease2_Disease1_params,force=TRUE)

Disease2_Disease1_res$loos[[2]]
# causal model
Disease2_Disease1_res$loos[[3]]
## results
# expected log pointwise posterior density
Disease2_Disease1_res$elpd %>% 
  mutate(
    pval = pnorm(z, lower.tail=TRUE)
  )
# summary
summary(Disease2_Disease1_res, ci_size=0.95)
save(Disease2_Disease1_res,file='cause_final.Rdata')