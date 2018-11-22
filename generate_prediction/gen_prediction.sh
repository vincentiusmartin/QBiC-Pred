#!/bin/bash  

#SBATCH -o vm_%A_%a.out
#SBATCH -e vm_%A_%a.err

#SBATCH --array=0-665%80
#SBATCH --mem=20G
#SBATCH --mail-type=END
#SBATCH --mail-user=vm76@duke.edu

# Directory setting -- parallelize per TF
#"/gpfs/fs0/data/gordanlab/vincentius" "/home/vincentius/Research/mutpred/Mutpred-all/preparator"
output="/data/gordanlab/vincentius/predmodel"
# (/home/vincentius/Research/mutpred/mutation-predictor/all-test2/*) #(/home/vm76/tf/*)
input=(/data/gordanlab/vincentius/pbmdata/*)
# Parameter configuration
kmer=6
chunk=32


# TO RUN: sbatch gen_prediction.sh 
#------------End of Configuration-----------------
 
filein=${input[SLURM_ARRAY_TASK_ID]}
echo $filein
echo "task_id: $SLURM_ARRAY_TASK_ID"

python olskmer.py $filein $output -k $kmer -d $chunk



