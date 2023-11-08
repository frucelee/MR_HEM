#!/bin/bash
#
#SBATCH --job-name=MR_risks_to_protein
#SBATCH --output=MR_mpi.txt
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=20000

cat risk_list | while read id
do
cat protein_list | while read id0
do
Rscript --vanilla MR_risk_to_protein.R ${id} ${id0}
cp mr.txt ${id0}_${id}
rm mr.txt
done
done