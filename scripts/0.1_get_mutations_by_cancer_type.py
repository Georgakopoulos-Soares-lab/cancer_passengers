import os
import sys
import pandas as pd

if len(sys.argv) < 2:
    raise ValueError("Usage: python process_subtype.py <SUBTYPE>")

SUBTYPE = sys.argv[1]

# input files
CANCER_SUBTYPE_COUNTS = "metadata/cancer_subtypes_counts.tsv"
PREP_DIR = "data/intermediate_files"

# output files/dir
MUT_DATA_BY_CANCER_SUBTYPE = "data/snv_mv_indels_by_cancer_subtype"
os.makedirs(MUT_DATA_BY_CANCER_SUBTYPE, exist_ok=True)

mut_data_df = pd.read_parquet(f"{PREP_DIR}/mut_data.parquet")
drivers_pcawg = pd.read_parquet(f"{PREP_DIR}/driver_mutations.parquet")

cancer_subtype_mut_data_df = mut_data_df[mut_data_df["Project_Code"] == SUBTYPE]
print(f"Number of mutations in {SUBTYPE} Cancer: {cancer_subtype_mut_data_df.shape[0]}")

# annotate PCAWG driver mutations
drivers = drivers_pcawg[drivers_pcawg["Project_Code"] == SUBTYPE]
drivers.loc[:, "driver"] = True
cancer_subtype_mut_data_df = pd.merge(cancer_subtype_mut_data_df, drivers, on=["Tumor_Sample_Barcode", "mutation_loc"], how="left")
cancer_subtype_mut_data_df["driver"].fillna(False, inplace=True)
cancer_subtype_mut_data_df.drop(["Project_Code_x", "Project_Code_y"], axis=1, inplace=True)

print(f"Saving {SUBTYPE} Cancer data for {cancer_subtype_mut_data_df.shape[0]} mutations")
cancer_subtype_mut_data_df.to_csv(f"{MUT_DATA_BY_CANCER_SUBTYPE}/{SUBTYPE}.tsv", sep="\t", index=False)
