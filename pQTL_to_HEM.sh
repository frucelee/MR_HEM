#!/bin/bash
#
#SBATCH --job-name=MR
#SBATCH --output=MR_mpi.txt
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000

cat /scratch/users/s/h/shifang/QTL/RAW/protein_list | while read id
do
awk '$8<=5e-08' /scratch/users/s/h/shifang/QTL/RAW/sig/${id%}.txt>${id}.subset.txt
awk '{print $4,$8}' /scratch/users/s/h/shifang/QTL/RAW/sig/${id}.subset.txt | sed 's/ /\t/g' >${id}.txt0
awk 'BEGIN {print "SNP   P"} {print $0}' ${id}.txt0 >${id}.txt1

##LD Pruning
plink --bfile /scratch/users/s/h/shifang/QTL/REF/EUR --clump ${id}.txt1  --clump-r2 0.001 --out MetaGWAS_CAD_clumped_r0.05 --clump-kb 10000 --clump-p1 1 --clump-p2 1 --threads 10
awk '{print $1,$3}' MetaGWAS_CAD_clumped_r0.05.clumped > SNP.valid
awk '{print $2}' SNP.valid >/scratch/users/s/h/shifang/QTL/RAW/sig/clump/${id}_clumped.txt
rm -r ${id}.txt1
rm -r ${id}.txt0
rm -r MetaGWAS_CAD_clumped_r0.05.clumped *.log *.nosex
rm -r SNP.valid
mv /scratch/users/s/h/shifang/QTL/RAW/sig/clump/${id}_clumped.txt /scratch/users/s/h/shifang/QTL/RAW/sig/clump/new/${id}.txt; 
done

##subset SNPs for MR
Rscript subset_SNP_of_pQTL.R

##performing MR analysis (pQTL_disease) using TwoSampleMR.
cat pQTL_list | while read id
do
Rscript --vanilla pQTL_to_HEM.R ${id}
done
