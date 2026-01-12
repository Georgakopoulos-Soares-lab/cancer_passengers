import os
import argparse
import pandas as pd

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--mut_file", type=str, help="Path to the mutation file")
	parser.add_argument("--results", type=str, help="Path to save the annotated mutation file")
	args = parser.parse_args()
	mut_file = args.mut_file
	results = args.results

	mut_df = pd.read_csv(mut_file, sep="\t")

	# generate input file for annovar
	mut_bed = pd.DataFrame()
	mut_bed["chrom"] = mut_df["mutation"].apply(lambda x: x.split(":")[0])
	mut_bed["start"] = mut_df["mutation"].apply(lambda x: int(x.split(":")[1]))
	mut_bed["end"] = mut_df["mutation"].apply(lambda x: int(x.split(":")[1]))
	mut_bed["ref"] = mut_df["mutation"].apply(lambda x: x.split(":")[2])
	mut_bed["alt"] = mut_df["mutation"].apply(lambda x: x.split(":")[3])
	mut_bed["end"] = mut_bed.apply(lambda x: x["start"] if x["ref"] == "-" else x["end"], axis=1) # handle insertions for annovar
	mut_bed.drop_duplicates(inplace=True)
	mut_bed.to_csv("mutations.input", sep="\t", index=False, header=False)

	# annotate mutations with genic regions using annovar
	input_file = "mutations.input"
	genome_db = "annovar/humandb/"
	cmd = f"annovar/annotate_variation.pl -out mutations -build hg19 {input_file} {genome_db}"
	os.system(cmd)

	# add genic region annotations from annovar to the mutation data
	annot_df = pd.read_csv("mutations.variant_function", sep="\t", header=None)
	annot_df = annot_df[[0, 1, 2, 3]]
	annot_df.columns = ["genic_region", "gene", "chr", "start"]
	annot_df["mutation_pos"] = annot_df["chr"].astype(str) + ":" + annot_df["start"].astype(str)
	annot_df = annot_df[["mutation_pos", "genic_region", "gene"]]
	annot_df.drop_duplicates(inplace=True)
	mut_df["mutation_pos"] = mut_df["mutation"].apply(lambda x: x.split(":")[0] + ":" + x.split(":")[1])
	mut_df = mut_df.merge(annot_df, on="mutation_pos", how="left")
	mut_df.drop(columns=["mutation_pos"], inplace=True)
	mut_df.to_csv(results, sep="\t", index=False)
