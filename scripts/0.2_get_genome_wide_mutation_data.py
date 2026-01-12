import os
import numpy as np
import pandas as pd

# input files/dir
MUTATION_DATA = "data/snv_mv_indels_by_cancer_subtype"
CNV_DATA = "data/datasets/PCAWG/mutations/cnv"
SV_DATA = "data/datasets/PCAWG/mutations/sv"

# output files/dir
GENOME_WIDE_MUT_DATA = "data/genome_wide_mutation_data.tsv"

CANCER_TYPES = os.listdir(MUTATION_DATA)
CANCER_TYPES = [cancer_type for cancer_type in CANCER_TYPES if cancer_type.endswith(".tsv")]
CANCER_TYPES = [cancer_type.replace(".tsv", "") for cancer_type in CANCER_TYPES]

# get genome wide SNPs, MNVs and Indels
genome_wide_mut_df = pd.DataFrame()
for cancer_type in CANCER_TYPES:
	df = pd.read_csv(os.path.join(MUTATION_DATA, cancer_type + ".tsv"), sep="\t")
	df["cancer_type"] = cancer_type
	df_grouped = df.groupby("Tumor_Sample_Barcode").agg({
		"mutation": "count",
		"cancer_type": "first"
	}).reset_index()
	df_grouped.rename(columns={"mutation": "total_mutations"}, inplace=True)
	genome_wide_mut_df = pd.concat([genome_wide_mut_df, df_grouped], axis=0, ignore_index=True)
print(f"Processed {len(genome_wide_mut_df)} samples for genome-wide mutation rate")

# add CNV burden data
cnv_data_files = os.listdir(CNV_DATA)
cnv_data_files = [cnv_data_file for cnv_data_file in cnv_data_files if cnv_data_file.endswith(".txt")]
cnv_data_df = []
for file in cnv_data_files:
	cnv_data = pd.read_csv(f"{CNV_DATA}/{file}", sep="\t")
	# get altered regions, 1/1 represents normal diploid state
	cnv_data = cnv_data[(cnv_data["major_cn"] != 1) | (cnv_data["minor_cn"] != 1)]
	cnv_data["cna_length"] = cnv_data["end"] - cnv_data["start"] + 1
	cnv_data_df.append({
		"Tumor_Sample_Barcode": file.split(".")[0],
		"alterned_region_length": cnv_data["cna_length"].sum() 
	})
cnv_data_df = pd.DataFrame(cnv_data_df)
cnv_data_df["cna_burden"] = cnv_data_df["alterned_region_length"] * 100 / 3.1e+09
cnv_data_df.drop(columns=["alterned_region_length"], inplace=True)
genome_wide_mut_df = pd.merge(genome_wide_mut_df, cnv_data_df, on="Tumor_Sample_Barcode", how="left")
print(f"Measured CNV burden data for {len(genome_wide_mut_df)} samples")
genome_wide_mut_df.to_csv(GENOME_WIDE_MUT_DATA, sep="\t", index=False)
