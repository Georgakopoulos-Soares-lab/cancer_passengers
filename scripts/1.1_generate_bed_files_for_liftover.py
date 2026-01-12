import os
import argparse
import pandas as pd 

# input files/dir
MUT_DATA = "data/annotated_snv_mv_indels_by_cancer_subtype"
DRIVER_GENES = "data/driver_genes"

# output files/dir
MUTATION_BED_HG19 = "data/mutation_bed_files_hg19"
MUTATION_BED_HG38 = "data/mutation_bed_files_hg38"

os.makedirs(MUTATION_BED_HG19, exist_ok=True)
os.makedirs(MUTATION_BED_HG38, exist_ok=True)

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
	mut_data["CHROM"] = "chr" + mut_data["CHROM"].astype(str)
	mut_data["POS"] = mut_data["POS"].astype(int) - 1 # convert to 0-based
	mut_data["END"] = mut_data["POS"]
	mut_data["mutation_loc"] = mut_data["CHROM"].str.replace("chr", "") + ":" + (mut_data["POS"] + 1).astype(str)
	mut_data = mut_data[["CHROM", "POS", "END", "mutation_loc"]]
	mut_data.drop_duplicates(inplace=True)
	print(mut_data.shape)

	# generate bed files
	mut_data.to_csv(f"{MUTATION_BED_HG19}/{cancer_type}.bed", sep="\t", index=False, header=False)
 
	# liftover to hg38 using UCSC liftover tool
	input_bed = f"{MUTATION_BED_HG19}/{cancer_type}.bed"
	output_bed = f"{MUTATION_BED_HG38}/{cancer_type}.bed"
	unmapped_bed = f"{MUTATION_BED_HG38}/{cancer_type}_unmapped.bed"
	chain_file = "metadata/hg19ToHg38.over.chain.gz"

	# get hg38 coordinates but also keep the hg19 coordinates in the name field
	os.system(f"./scripts/liftOver {input_bed} {chain_file} {output_bed} {unmapped_bed}")
 
	# remove empty unmapped files
	if os.path.getsize(unmapped_bed) == 0:
		os.remove(unmapped_bed)
		print("No unmapped regions")
	else:
		print(f"Unmapped regions present, check {unmapped_bed}")
