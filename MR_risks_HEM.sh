#!/bin/bash
#
#SBATCH --job-name=MR_risks_HEM
#SBATCH --output=MR_mpi.txt
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000

##identified the causal effect using Twosample MR
Rscript MR_risks_to_HEM.R

##identified the reverse causal effect using Twosample MR
Rscript MR_HEM_to_risks.R

##identified the causal effect using CAUSE model

#Effect of risk factors (exposure) on HEM (outcome)
cat /scratch/users/s/h/shifang/QTL/MR/risk_factors_list | while read id
do

##Format data for CAUSE and calculate nuisance parameters
Rscript --vanilla cause_MR_stage1.R ${id}

#here plink was used for LD Pruning (r2 0.001; p-value 1e-03)
plink --bfile /scratch/users/s/h/shifang/QTL/REF/EUR --clump SNP.txt  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 1e-3 --clump-p2 1e-3 --threads 10
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >HO_clumped_raw.txt
sed '$d' HO_clumped_raw.txt | sed '$d' >HO_clumped.txt
rm -r MetaGWAS_CAD_clumped_r0.05.clumped *.log *.nosex
rm -r SNP.valid HO_clumped_raw.txt

#Fit CAUSE and save results
Rscript --vanilla cause_MR_stage2.R ${id}
cp cause_final.Rdata cause_${id}_to_HEM.Rdata
rm -r cause_final.Rdata
done

#Effect of HEM (exposure) on risk factors (outcome)
cat /scratch/users/s/h/shifang/QTL/MR/risk_factors_list | while read id
do

##Format data for CAUSE and calculate nuisance parameters
Rscript --vanilla Rev_cause_MR_stage1.R ${id}

#change the order between exposure and outcome
cp SNP_1.txt SNP.txt

#here plink was used for LD Pruning (r2 0.001; p-value 1e-03)
plink --bfile /scratch/users/s/h/shifang/QTL/REF/EUR --clump SNP.txt  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 1e-3 --clump-p2 1e-3 --threads 10
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >HO_clumped_raw.txt
sed '$d' HO_clumped_raw.txt | sed '$d' >HO_clumped.txt
rm -r MetaGWAS_CAD_clumped_r0.05.clumped *.log *.nosex
rm -r SNP.valid HO_clumped_raw.txt

#Fit CAUSE and save results
Rscript --vanilla cause_MR_stage2.R ${id}
cp cause_final.Rdata cause_HEM_to_${id}.Rdata
rm -r cause_final.Rdata
done