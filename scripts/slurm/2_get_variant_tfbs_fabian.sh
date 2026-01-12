#!/bin/bash
#SBATCH --chdir=/storage/group/izg5139/default/akshatha/cancer_mutation_model
#SBATCH -o /storage/group/izg5139/default/akshatha/cancer_mutation_model/logs/fabian.out
#SBATCH -e /storage/group/izg5139/default/akshatha/cancer_mutation_model/logs/fabian.err
#SBATCH --account=izg5139_bc
#SBATCH --partition=sla-prio
#SBATCH --mem=16G
#SBATCH --time=4-00:00:00

source /storage/home/abn5461/work/miniforge3/bin/activate /storage/home/abn5461/work/miniforge3/envs/cancer-model 

files="data/driver_genes"
CANCER_TYPES=$(ls $files | cut -f 1 | sort | uniq)
CANCER_TYPES=$(echo $CANCER_TYPES | sed 's/.tsv//g')
echo $CANCER_TYPES

for CANCER_TYPE in ${CANCER_TYPES[@]}; do
    echo $CANCER_TYPE
    python scripts/3_get_fabian_input.py --cancer_type $CANCER_TYPE
    bash scripts/3.1_run_fabian.sh $CANCER_TYPE
    python scripts/3.2_process_fabian_output.py --cancer_type $CANCER_TYPE
done
