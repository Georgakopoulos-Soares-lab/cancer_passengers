import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# input files/dir
ICGC_DRIVER_MUTATIONS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_ICGC_samples.public.tsv"
TCGA_DRIVER_MUTATIONS = "data/datasets/PCAWG/driver_mutations/TableS3_panorama_driver_mutations_TCGA_samples.controlled.tsv"
HGNC_GENES = "metadata/HGNC_gene_names.tsv"

# output files/dir
REGION_GENE_MAP = "metadata/cytobands.tsv"

driver_df = pd.concat(
    [
		pd.read_csv(ICGC_DRIVER_MUTATIONS, sep="\t"),
		pd.read_csv(TCGA_DRIVER_MUTATIONS, sep="\t"),
	]
)
driver_df = driver_df[["sample_id", "ttype", "gene", "top_category"]]

# HGNC gene names
gene_names = pd.read_csv(HGNC_GENES, sep="\t")
gene_names = gene_names.rename(columns={"Approved symbol": "gene"})["gene"].unique().tolist()

# HGNC gene aliases
aliases = pd.read_csv(HGNC_GENES, sep="\t")
aliases = aliases.rename(columns={"Alias symbols": "aliases"})["aliases"]
aliases = aliases.dropna().unique().tolist()
aliases = [alias for alias in aliases if alias != ""]
alias_list = []
for alias in aliases:
    alias_list.extend(alias.split(", "))
alias_list = list(set(alias_list))

# HGNC previous symbols
previous_symbols = pd.read_csv(HGNC_GENES, sep="\t")
previous_symbols = previous_symbols.rename(columns={"Previous symbols": "previous"})["previous"]
previous_symbols = previous_symbols.dropna().unique().tolist()
previous_symbols = [symbol for symbol in previous_symbols if symbol != ""]
print(len(previous_symbols), previous_symbols[:10])
previous_list = []
for symbol in previous_symbols:
	previous_list.extend(symbol.split(", "))
previous_list = list(set(previous_list))

# get entries that are not HGNC gene names
driver_gene_count = driver_df.groupby(["gene", "top_category"]).size().reset_index(name="count")
driver_gene_count = driver_gene_count.sort_values(by="count", ascending=False)
driver_gene_count = driver_gene_count[~driver_gene_count["gene"].isin(gene_names)]
gene_alias = driver_gene_count[driver_gene_count["gene"].isin(alias_list)]
driver_gene_count = driver_gene_count[~driver_gene_count["gene"].isin(alias_list)]
gene_previous = driver_gene_count[driver_gene_count["gene"].isin(previous_list)]
driver_gene_count = driver_gene_count[~driver_gene_count["gene"].isin(previous_list)]

# exclude telomeres as they do not have genes
driver_gene_count = driver_gene_count[~driver_gene_count["gene"].str.contains("telomere")]
print(driver_gene_count.shape)
# save file for manual annotation
region_gene = driver_gene_count.rename(columns={"gene": "region"})[["region"]]
region_gene.to_csv(REGION_GENE_MAP, sep="\t", index=False)
