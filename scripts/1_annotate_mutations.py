import os
import argparse
import pandas as pd
from pandarallel import pandarallel

# input files/dir
MUT_DATA_BY_CANCER_SUBTYPE = "data/snv_mv_indels_by_cancer_subtype"
CADD_SCORE_FILE = "cadd/whole_genome_SNVs.tsv.gz"

# intermediate (input for annovar)
MUT_ANNOT = "data/mutations_annot"
os.makedirs(MUT_ANNOT, exist_ok=True)

# output files/dir
ANNOTATED_MUTATIONS = "data/annotated_snv_mv_indels_by_cancer_subtype"

if not os.path.exists(MUT_ANNOT):
    os.makedirs(MUT_ANNOT)

if not os.path.exists(ANNOTATED_MUTATIONS):
    os.makedirs(ANNOTATED_MUTATIONS)

def get_genic_region_annotations(cancer_type):
    print(f"Processing {cancer_type}...")
    # annotate mutations with driver/passenger status
    mut_df = pd.read_csv(f"{MUT_DATA_BY_CANCER_SUBTYPE}/{cancer_type}.tsv", sep="\t")

    # generate input file for annovar
    mut_df["mutation_pos"] = mut_df["mutation"].apply(lambda x: ':'.join(x.split(":")[:2]))
    mut_bed = mut_df[["mutation_pos", "Strand"]]
    mut_bed["chrom"] = mut_bed["mutation_pos"].apply(lambda x: x.split(":")[0])
    mut_bed["start"] = mut_bed["mutation_pos"].apply(lambda x: int(x.split(":")[1].split("-")[0]))
    mut_bed["end"] = mut_bed["mutation_pos"].apply(lambda x: int(x.split(":")[1].split("-")[1]))
    mut_bed["ref"] = mut_df["mutation"].apply(lambda x: x.split(":")[2])
    mut_bed["alt"] = mut_df["mutation"].apply(lambda x: x.split(":")[3])
    mut_bed["end"] = mut_bed.apply(lambda x: x["start"] if x["ref"] == "-" else x["end"], axis=1) # handle insertions for annovar
    mut_bed = mut_bed[["chrom", "start", "end", "ref", "alt"]]
    mut_bed.drop_duplicates(inplace=True)
    mut_bed.to_csv(f"{MUT_ANNOT}/{cancer_type}.input", sep="\t", index=False, header=False)
    print(f"Generated input file for annovar for {cancer_type}")

    # annotate mutations with genic regions using annovar
    input_file = f"{MUT_ANNOT}/{cancer_type}.input"
    output_prefix = f"{MUT_ANNOT}/{cancer_type}"
    genome_db = "annovar/humandb/"
    cmd = f"annovar/annotate_variation.pl -out {output_prefix} -build hg19 {input_file} {genome_db}"
    os.system(cmd)

    # add genic region annotations from annovar to the mutation data
    annot_file = f"{MUT_ANNOT}/{cancer_type}.variant_function"
    annot_df = pd.read_csv(annot_file, sep="\t", header=None)
    annot_df = annot_df[[0, 2, 3]]
    annot_df.columns = ["genic_region", "chr", "start"]
    annot_df["mutation_pos"] = annot_df["chr"].astype(str) + ":" + annot_df["start"].astype(str)
    annot_df = annot_df[["mutation_pos", "genic_region"]]
    annot_df.drop_duplicates(inplace=True)
    mut_df["mutation_pos"] = mut_df["mutation"].apply(lambda x: x.split(':')[0] + ":" + x.split(':')[1].split('-')[0])
    mut_df = mut_df.merge(annot_df, on="mutation_pos", how="left")
    mut_df.drop(columns=["mutation_pos"], inplace=True)
    print(f"Added genic region annotations for {cancer_type}")
    return mut_df

def get_cadd_scores(mutation):
    chrom, pos, ref, alt = mutation.split(":")
    # handle indels
    if ref == "-" or alt == "-":
        return [None, None]
    # handle MNVs
    if len(ref) > 1 or len(alt) > 1:
        return [None, None]
    loc = pos.split("-")[0]
    cmd = f"tabix {CADD_SCORE_FILE} {chrom}:{pos} | grep '{chrom}\t{loc}\t{ref}\t{alt}' | cut -f 5,6"
    cadd_score = os.popen(cmd).read().strip().split("\t")
    return cadd_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cancer_type", type=str, help="Cancer tissue type e.g Pancreas")
    args = parser.parse_args()
    cancer_type = args.cancer_type

    # get genic region annotations from annovar
    mut_df = get_genic_region_annotations(cancer_type)

    # get CADD scores for mutations
    pandarallel.initialize()    
    cadd_df = mut_df[["mutation"]].copy()
    cadd_df.drop_duplicates(inplace=True)
    cadd_df["CADD_scores"]  = cadd_df["mutation"].parallel_apply(get_cadd_scores)
    cadd_df["CADD_score_raw"] = cadd_df["CADD_scores"].apply(lambda x: x[0])
    cadd_df["CADD_score_PHRED"] = cadd_df["CADD_scores"].apply(lambda x: x[1])
    cadd_df.drop(columns=["CADD_scores"], inplace=True)
    mut_df = pd.merge(mut_df, cadd_df, on="mutation", how="left")
    print(f"Added CADD scores for {cancer_type}")
    mut_df.to_csv(f"{ANNOTATED_MUTATIONS}/{cancer_type}.tsv", sep="\t", index=False)
