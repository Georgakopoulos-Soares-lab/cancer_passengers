import os
import pandas as pd

# input files
ICGC_MUT_DATA = "data/datasets/PCAWG/mutations/snv_mnv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf"
TCGA_MUT_DATA = "data/datasets/PCAWG/mutations/snv_mnv_indel/final_consensus_passonly.snv_mnv_indel.tcga.controlled.maf"
ICGC_DRIVER_MUTATIONS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"
TCGA_DRIVER_MUTATIONS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_TCGA_samples.controlled.tsv"

# output files/dirs
PREP_DIR = "data/intermediate_files"
os.makedirs(PREP_DIR, exist_ok=True)
CANCER_SUBTYPE_COUNTS = "metadata/cancer_subtypes_counts.tsv"

# mutation data from ICGC and TCGA
if not os.path.exists(f"{PREP_DIR}/mut_data.parquet"):
	icgc_mut_data_df = pd.read_csv(ICGC_MUT_DATA, sep="\t")
	tcga_mut_data_df = pd.read_csv(TCGA_MUT_DATA, sep="\t")
	mut_data_df = pd.concat([icgc_mut_data_df, tcga_mut_data_df], ignore_index=True)
	mut_data_df = mut_data_df[["Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Variant_Classification", "Tumor_Sample_Barcode", "Project_Code", "Donor_ID"]]
	mut_data_df.rename(columns={"Donor_ID": "Patient_ID"}, inplace=True)
	mut_data_df["mutation"] = mut_data_df["Chromosome"].astype(str) + ":" + mut_data_df["Start_position"].astype(str) + "-" + mut_data_df["End_position"].astype(str) + ":" + mut_data_df["Reference_Allele"] + ":" + mut_data_df["Tumor_Seq_Allele2"]
	mut_data_df["mutation_loc"] = mut_data_df["Chromosome"].astype(str) + ":" + mut_data_df["Start_position"].astype(str) + ":" + mut_data_df["Reference_Allele"] + ":" + mut_data_df["Tumor_Seq_Allele2"]
	mut_data_df.drop(["Chromosome", "Start_position", "End_position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"], axis=1, inplace=True)
	mut_data_df.rename(columns={"Hugo_Symbol": "gene"}, inplace=True)
	print("Number of mutations:", mut_data_df.shape[0])
	print("Saving all mutation data...")
	mut_data_df.to_parquet(f"{PREP_DIR}/mut_data.parquet", index=False)

# get driver mutations from PCAWG resource (ICGC samples)
if not os.path.exists(f"{PREP_DIR}/driver_mutations.parquet"):
	drivers_icgc = pd.read_csv(ICGC_DRIVER_MUTATIONS, sep="\t")
	drivers_tcga = pd.read_csv(TCGA_DRIVER_MUTATIONS, sep="\t")
	drivers_pcawg = pd.concat([drivers_icgc, drivers_tcga])
	drivers_pcawg = drivers_pcawg[["sample_id", "ttype", "chr", "pos", "ref", "alt", "top_category"]]
	drivers_pcawg = drivers_pcawg[drivers_pcawg["top_category"].isin(["mutational", "germline"])]
	drivers_pcawg["mutation_loc"] = drivers_pcawg["chr"].astype(str) + ":" + drivers_pcawg["pos"].astype(str) + ":" + drivers_pcawg["ref"] + ":" + drivers_pcawg["alt"]
	drivers_pcawg.drop(["chr", "pos", "ref", "alt", "top_category"], axis=1, inplace=True)
	drivers_pcawg.rename(columns={"sample_id": "Tumor_Sample_Barcode", "ttype": "Project_Code"}, inplace=True)
	drivers_pcawg.drop_duplicates(inplace=True)
	print("Saving PCAWG driver mutation data...")
	drivers_pcawg.to_parquet(f"{PREP_DIR}/driver_mutations.parquet", index=False)

# Summary: cancer subtype counts
if not os.path.exists(CANCER_SUBTYPE_COUNTS):
	mut_data_df = pd.read_parquet(f"{PREP_DIR}/mut_data.parquet")
	cancer_subtype_counts = mut_data_df.groupby("Project_Code")["Tumor_Sample_Barcode"].nunique().reset_index()
	cancer_subtype_counts.rename(columns={"Tumor_Sample_Barcode": "Sample_Count"}, inplace=True)
	cancer_subtype_counts.sort_values("Sample_Count", ascending=False, inplace=True)
	print("Saving cancer subtype counts")
	cancer_subtype_counts.to_csv(CANCER_SUBTYPE_COUNTS, sep="\t", index=False)
	print(f"Total number of samples: {cancer_subtype_counts['Sample_Count'].sum()}")
