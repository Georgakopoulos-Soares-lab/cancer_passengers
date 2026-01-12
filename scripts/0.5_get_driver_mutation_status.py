import os
import argparse
import numpy as np
import pandas as pd

# input files/dir
MUT_DATA_BY_CANCER_SUBTYPE = "data/snv_mv_indels_by_cancer_subtype"
ICGC_DRIVERS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"
TCGA_DRIVERS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_TCGA_samples.controlled.tsv"
CYTOBAND_GENES = "metadata/cytoband_gene_map.tsv"

# output files/dir
MUT_DATA_BY_CANCER_SUBTYPE = "data/snv_mv_indels_by_cancer_subtype"
os.makedirs(MUT_DATA_BY_CANCER_SUBTYPE, exist_ok=True)
DRIVER_COUNT_BY_CANCER_SUBTYPE = "data/driver_count_by_cancer_subtype"
os.makedirs(DRIVER_COUNT_BY_CANCER_SUBTYPE, exist_ok=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer_type", type=str, help="Cancer tissue type e.g Pancreas")
    args = parser.parse_args()
    cancer_type = args.cancer_type

    # get driver mutation status
    drivers_icgc = pd.read_csv(ICGC_DRIVERS, sep="\t")
    drivers_tcga = pd.read_csv(TCGA_DRIVERS, sep="\t")
    drivers_pcawg = pd.concat([drivers_icgc, drivers_tcga])
    drivers_pcawg = drivers_pcawg[["sample_id", "ttype", "gene", "top_category"]]

    # driver mutation count
    drivers_count = drivers_pcawg.groupby("sample_id").agg({
        "top_category": "count"
    }).reset_index()
    drivers_count.rename(columns={"top_category": "num_drivers"}, inplace=True)
    drivers_count.sort_values(by="num_drivers", ascending=False, inplace=True)
    drivers_count.to_csv(f"{DRIVER_COUNT_BY_CANCER_SUBTYPE}/{cancer_type}.tsv", sep="\t", index=False)

    drivers_pcawg.rename(columns={
        "sample_id": "Tumor_Sample_Barcode", 
        "ttype": "Project_Code",
        "top_category": "driver_mutation_type"
    }, inplace=True)
    drivers_pcawg.drop_duplicates(inplace=True)
    drivers_pcawg = drivers_pcawg[drivers_pcawg["Project_Code"] == cancer_type]
    drivers_pcawg.drop(columns=["Project_Code"], inplace=True)
    drivers = drivers_pcawg.groupby(["Tumor_Sample_Barcode", "gene"]).agg({
        "driver_mutation_type": lambda x: ", ".join(x)
    }).reset_index()

    # add genes effected by large mutations
    cytoband_genes = pd.read_csv(CYTOBAND_GENES, sep="\t")
    cytoband_genes = cytoband_genes[["region", "genes"]]
    drivers = pd.merge(drivers, cytoband_genes, left_on="gene", right_on="region", how="left")
    drivers.drop(columns=["region"], inplace=True)
    drivers.rename(columns={"genes": "cytoband_genes"}, inplace=True)
    drivers["gene"] = drivers.apply(lambda x: x["cytoband_genes"] if pd.notnull(x["cytoband_genes"]) else x["gene"], axis=1)
    drivers.drop(columns=["cytoband_genes"], inplace=True)
    drivers = drivers.assign(gene=drivers["gene"].str.split(",")).explode("gene")
    drivers.drop_duplicates(inplace=True)

    # get mutations for cancer type
    mut_df = pd.read_csv(f"{MUT_DATA_BY_CANCER_SUBTYPE}/{cancer_type}.tsv", sep="\t")
    tumor_gene = mut_df.groupby(["Tumor_Sample_Barcode", "gene"]).agg({
        "driver": "any"
    }).reset_index()
    tumor_gene["driver"] = tumor_gene["driver"].apply(lambda x: "mutational" if x else np.nan)
    tumor_gene = tumor_gene.merge(drivers, on=["Tumor_Sample_Barcode", "gene"], how="left")
    tumor_gene["driver_mutation_type"] = tumor_gene["driver"] + ", " + tumor_gene["driver_mutation_type"]
    tumor_gene["driver_mutation_type"] = tumor_gene["driver_mutation_type"].apply(lambda x: ', '.join(list(set(x.split(", ")))) if type(x) == str else x)

    tumor_gene.drop(columns=["driver"], inplace=True)
    tumor_gene["driver_mutation_type"].fillna("None", inplace=True)

    tumor_gene["has_driver"] = tumor_gene["driver_mutation_type"].apply(lambda x: True if x != "None" else False)
    mut_df = mut_df.merge(tumor_gene, on=["Tumor_Sample_Barcode", "gene"], how="left")
    mut_df.to_csv(f"{MUT_DATA_BY_CANCER_SUBTYPE}/{cancer_type}.tsv", sep="\t", index=False)
    print(mut_df.shape)
    print(f"Added driver mutation status for {cancer_type}")
