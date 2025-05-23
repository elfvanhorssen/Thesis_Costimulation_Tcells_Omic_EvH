#!/bin/bash
#SBATCH --job-name=diffbind
#SBATCH --account=div4-ccb
#SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=40gb
#SBATCH --time=1:05:00
#SBATCH --output=diffbind_atac_all.out
#SBATCH --error=diffbind_atac_all.err
# change to working directory default for SLURM
#SBATCH -D ./

echo "SLURM_JOB_ACCOUNT="${SLURM_JOB_ACCOUNT}
echo "SLURM_JOB_PARTITION="${SLURM_JOB_PARTITION}
echo "SLURM_JOB_QOS="${SLURM_JOB_QOS}
echo "SLURM_JOBID="${SLURM_JOBID}
echo "SLURM_JOB_NODELIST"=${SLURM_JOB_NODELIST}
echo "SLURM_NNODES"=${SLURM_NNODES}
echo "SLURM_NTASKS"=${SLURM_NTASKS}
echo "SLURM_TASKS_PER_NODE"=${SLURM_TASKS_PER_NODE}
echo "working directory = "$SLURM_SUBMIT_DIR

echo "start job DiffBind `date`"

conda activate /exports/chemimm/minkang/tools/miniconda3/envs/r4-env
Rscript plotPCA.r
#Rscript diffbind_atac.R /exports/archive/chemimm/Minkang_archived/Felix_ATACseq/105042-003/DiffBind/ /exports/archive/chemimm/Minkang_archived/Felix_ATACseq/105042-003/DiffBind/sample_input.csv 
echo "Finished `date`"

