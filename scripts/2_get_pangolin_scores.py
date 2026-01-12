import os
import argparse
import pandas as pd 

# input files/dir
MUT_DATA = "data/annotated_snv_mv_indels_by_cancer_subtype"
DRIVER_GENES = "data/driver_genes"
REF_FASTA = "data/ref/hg19_uppercase.fa"
GENE_ANNOT_DB = "data/ref/gencode.v38lift37.annotation.db"

# output files/dir
SPLICE_PRED = "data/pangolin_predictions"

if not os.path.exists(SPLICE_PRED):
    os.makedirs(SPLICE_PRED)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer_type", type=str, help="Cancer tissue type e.g Pancreas")
    args = parser.parse_args()
    cancer_type = args.cancer_type
    print(f"Processing {cancer_type}...")

    # get driver genes
    driver_genes = pd.read_csv(f"{DRIVER_GENES}/{cancer_type}.tsv", sep="\t")
    driver_genes_list = driver_genes["gene"].tolist()

    # get mutations in driver genes
    mut_data = pd.read_csv(f"{MUT_DATA}/{cancer_type}.tsv", sep="\t")
    mut_data = mut_data[mut_data["gene"].isin(driver_genes_list)]
    print(mut_data.shape)
    mut_data = mut_data[["mutation_loc"]]
    mut_data[["CHROM", "POS", "REF", "ALT"]] = mut_data["mutation_loc"].str.split(":", expand=True)
    mut_data.drop(columns=["mutation_loc"], inplace=True)
    mut_data.drop_duplicates(inplace=True)
    print(mut_data.shape)
    # filter indels as they are not supported by pangolin
    mut_data = mut_data[(mut_data["REF"] != "-") & (mut_data["ALT"] != "-")] 
    print(f"Number of SNPs in driver genes: {mut_data.shape[0]}")

    # generate input file for pangolin
    mut_data.to_csv(f"{SPLICE_PRED}/{cancer_type}_input.csv", sep=",", index=False)

    # run pangolin
    input_file = f"{SPLICE_PRED}/{cancer_type}_input.csv"
    output_file = f"{SPLICE_PRED}/{cancer_type}_pred"
    cmd = f"pangolin {input_file} {REF_FASTA} {GENE_ANNOT_DB} {output_file}"
    print(cmd)
    os.system(cmd)
    print(f"Finished fetching pangoling scores for {cancer_type} cancer")
