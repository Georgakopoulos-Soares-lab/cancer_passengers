import pandas as pd

CYTOBANDS = "metadata/cytobands.tsv"
CYTOBAND_GENES = "metadata/cytoband_gene_map.tsv"
GENCODE_ANNOT = "data/ref/gencode.v19.annotation.gtf"

# get genes from gencode
gencode_annot = pd.read_csv(GENCODE_ANNOT, sep="\t", comment="#", header=None)
gencode_annot = gencode_annot[gencode_annot[2] == "gene"]
gencode_annot = gencode_annot.sort_values(by=[0, 3, 4]).reset_index(drop=True)

# get gene coordinates
gencode_genes = gencode_annot[[0, 3, 4, 6, 8]].copy()
gencode_genes["gene"] = gencode_genes[8].str.extract(r'gene_name "(.*?)";')
gencode_genes.drop(8, axis=1, inplace=True)
gencode_genes.columns = ["chr", "start", "end", "strand", "gene"]
gencode_genes["chr"] = gencode_genes["chr"].str.replace("chr", "")
gencode_genes["start"] = gencode_genes["start"].astype(int)
gencode_genes["end"] = gencode_genes["end"].astype(int)
gencode_genes = gencode_genes.sort_values(by=["chr", "start", "end"]).reset_index(drop=True)
gencode_genes = gencode_genes.reset_index(drop=True)

# get cytoband coordinates
cytobands = pd.read_csv(CYTOBANDS, sep="\t")
cytobands["chr"] = cytobands["coordinates"].str.extract(r'chr(.*?):')
cytobands["start"] = cytobands["coordinates"].str.extract(r':(.*)-').astype(int)
cytobands["end"] = cytobands["coordinates"].str.extract(r'-(.*)').astype(int)
cytobands.drop("coordinates", axis=1, inplace=True)
cytobands = cytobands[["chr", "start", "end", "region"]]
cytobands.sort_values(by=["chr", "start", "end"], inplace=True)
cytobands = cytobands.reset_index(drop=True)

# get genes in each cytoband
cytoband_genes = pd.DataFrame(columns=["chr", "start", "end", "region", "gene"])
for i, row in cytobands.iterrows():
	chr_ = row["chr"]
	start = row["start"]
	end = row["end"]
	region = row["region"]
	genes = gencode_genes[(gencode_genes["chr"] == chr_) & (gencode_genes["start"] >= start) & (gencode_genes["end"] <= end)]
	genes["region"] = region
	cytoband_genes = pd.concat([cytoband_genes, genes], ignore_index=True)

# list of all genes in each region
region_genes = cytoband_genes[["region", "gene"]].groupby("region")["gene"].apply(lambda x: ",".join(x)).reset_index()
region_genes.columns = ["region", "genes"]
region_genes["num_genes"] = region_genes["genes"].apply(lambda x: len(x.split(",")))

# add genes to cytoband data and save to file
cytoband_genes = pd.read_csv(CYTOBAND_GENES, sep="\t")
cytoband_genes.drop("genes", axis=1, inplace=True)
cytoband_genes = pd.merge(cytoband_genes, region_genes, on="region", how="left")
cytoband_genes.to_csv(CYTOBAND_GENES, sep="\t", index=False)
