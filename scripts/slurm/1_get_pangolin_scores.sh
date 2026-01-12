#!/bin/bash
#SBATCH --chdir=/storage/group/izg5139/default/akshatha/cancer_mutation_model
#SBATCH -o /storage/group/izg5139/default/akshatha/cancer_mutation_model/logs/jobs/out_%a.log
#SBATCH -e /storage/group/izg5139/default/akshatha/cancer_mutation_model/logs/jobs/err_%a.log
#SBATCH --account=izg5139_bc
#SBATCH --partition=sla-prio
#SBATCH --mem=16G
#SBATCH --array 1-31
#SBATCH --time=5-00:00:00

source /storage/home/abn5461/work/miniforge3/bin/activate /storage/home/abn5461/work/miniforge3/envs/cancer-model 

# path to logs
log_path=logs/pangolin_scores
mkdir -p $log_path

files="data/driver_genes"
CANCER_TYPES=$(ls $files | cut -f 1 | sort | uniq)
CANCER_TYPES=$(echo $CANCER_TYPES | sed 's/.tsv//g')
echo $CANCER_TYPES

cancer_type=$(echo $CANCER_TYPES | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)
python scripts/2_get_pangolin_scores.py --cancer_type $cancer_type > $log_path/$cancer_type.log
