#/bin/bash
set -x

files="data/driver_genes"
CANCER_TYPES=$(ls $files | cut -f 1 | sort | uniq)
CANCER_TYPES=$(echo $CANCER_TYPES | sed 's/.tsv//g')
echo $CANCER_TYPES

# run notebooks for each cancer type
cd notebooks
for CANCER_TYPE in ${CANCER_TYPES[@]}; do
    echo $CANCER_TYPE
    papermill 1.1_mutation_density_genic_region_by_cancer.ipynb temp.ipynb -p cancer_type $CANCER_TYPE
    papermill 2.1_cadd_score_genic_region_by_cancer.ipynb temp.ipynb -p cancer_type $CANCER_TYPE
    papermill 6.2_splicing_analysis_pangolin_by_cancer.ipynb temp.ipynb -p cancer_type $CANCER_TYPE
    papermill 7.2_TFB_effect_differential_by_cancer.ipynb temp.ipynb -p cancer_type $CANCER_TYPE
done