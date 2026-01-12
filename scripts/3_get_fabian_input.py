import os
import argparse
import pyfaidx
import pysam
import pandas as pd 

# input files/dir
MUT_DATA = "data/annotated_snv_mv_indels_by_cancer_subtype"
DRIVER_GENES = "data/driver_genes"
REF_FASTA = "data/ref/hg19.fa"

# output files/dir
FABIAN_INPUT = "data/fabian_input"

if not os.path.exists(FABIAN_INPUT):
    os.makedirs(FABIAN_INPUT)

# Fasta sequence extractor
pyfaidx.Faidx(REF_FASTA)
REF_GENOME_OPEN = pysam.Fastafile(REF_FASTA) 

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
    print(cancer_type, mut_data.shape)
    mut_data = mut_data[["mutation_loc"]]
    mut_data.drop_duplicates(inplace=True)
    print(cancer_type, mut_data.shape)
    
    # create vcf input
    vcf_file = f"{FABIAN_INPUT}/{cancer_type}_input.vcf"
    with open(vcf_file, "w") as f:
        f.write(f"##fileformat=VCFv4.2\n")
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDATA\n")
        for index, row in mut_data.iterrows():
            chrom, pos, ref, alt = row["mutation_loc"].split(":")
            # using 0-based indexing for pysam fetch
            start = int(pos)-1
            # include the previous base in case of indels and MNVs
            if ref == "-" or alt == "-" or len(ref) > 1 or len(alt) > 1:
                if ref == "-":
                    ref = ""
                if alt == "-":
                    alt = ""
                ref = REF_GENOME_OPEN.fetch(f"chr{chrom}", start-1, start) + ref
                alt = REF_GENOME_OPEN.fetch(f"chr{chrom}", start-1, start) + alt
                ref = ref.upper()
                alt = alt.upper()
                pos = int(pos) - 1
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t116\t.\t.\tGT:DP\t0/1:154\n")
