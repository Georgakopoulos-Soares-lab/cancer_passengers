#/bin/bash
set -x

# Generate intermediate files for processing mutation data
python scripts/0_prep_data.py

# Get mutation data for each cancer type
cancer_subtypes_list="metadata/cancer_subtypes_counts.tsv"
CANCER_TYPES=$(cut -f 1 $cancer_subtypes_list | tail -n +2 | sort | uniq)
CANCER_TYPES_ARRAY=($CANCER_TYPES)
echo ${#CANCER_TYPES_ARRAY[@]} ${CANCER_TYPES_ARRAY[@]}
for CANCER_TYPE in ${CANCER_TYPES_ARRAY[@]}; do
    python scripts/0.1_get_mutations_by_cancer_type.py $CANCER_TYPE
done

# Get genome-wide mutation data
python scripts/0.2_get_genome_wide_mutation_data.py

# Get the top 10 driver genes for each cancer type along with coordinate information
python scripts/0.3_get_top_driver_genes.py

# Get genes overlapping cytoband regions
python scripts/0.4.1_get_cytoband_regions.py
python scripts/0.4.2_get_cytoband_genes.py

# Process mutation data for each cancer type
files="data/snv_mv_indels_by_cancer_subtype"
CANCER_TYPES=$(ls $files | cut -f 1 | sort | uniq)
CANCER_TYPES=$(echo $CANCER_TYPES | sed 's/.tsv//g')
echo $CANCER_TYPES
for CANCER_TYPE in ${CANCER_TYPES[@]}; do
    echo $CANCER_TYPE

    # Annotate mutations with driver status
    python scripts/0.5_get_driver_mutation_status.py --cancer_type $CANCER_TYPE

    # Annotate mutations with various genomic features using ANNOVAR and with pathogenicity scores from CADD
    python scripts/1_annotate_mutations.py --cancer_type $CANCER_TYPE

    # Generate BED files for liftover from hg19 to hg38 (required for certain tools)
    python scripts/1.1_generate_bed_files_for_liftover.py --cancer_type $CANCER_TYPE

    # Get splicing prediction from Pangolin tool 
    # Recommended to run on HPC cluster
    # Slurm script available at scripts/slurm/1_get_pangolin_scores.sh
    python scripts/2_get_pangolin_scores.py --cancer_type $CANCER_TYPE

    # Get TF binding prediction from FABIAN-Variant tool
    # Recommended to run on HPC cluster
    # Slurm script available at scripts/slurm/2_get_variant_tfbs_fabian.sh
    python scripts/3_get_fabian_input.py --cancer_type $CANCER_TYPE
    bash scripts/3.1_run_fabian.sh $CANCER_TYPE
    python scripts/3.2_process_fabian_output.py --cancer_type $CANCER_TYPE
done
