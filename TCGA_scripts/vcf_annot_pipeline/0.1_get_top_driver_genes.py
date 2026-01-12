import os
import pandas as pd

# input files/dir
DRIVER_GENES_INTOGEN = "data/datasets/driver_genes_TCGA_intogen"
SAMPLE_SHEET = "TCGA_scripts/vcf_annot_pipeline/gdc_sample_sheet.2025-11-04.tsv"
REF_DIR = "data/ref"
REF_GENOME = f"{REF_DIR}/hg19.fa"
GENCODE_ANNOT = f"{REF_DIR}/gencode.v19.annotation.gtf" # for getting start and end positions of genes
CHROM_SIZES = f"{REF_DIR}/hg19.chrom.sizes"
REFGENE_ANNOT = f"{REF_DIR}/hg19.refGene.gtf" # for getting genic region annotations
INTRON_ANNOTATION = f"{REF_DIR}/hg19_gene_annot_introns.bed"

# output files/dir
GENE_ANNOT_BED_DIR = "data/tcga_genic_region_annot_bed"
DRIVER_GENES_DIR = "data/TCGA_driver_genes"
GENE_LIST = "metadata/tcga_gene_list.tsv"

if not os.path.exists(GENE_ANNOT_BED_DIR):
    os.makedirs(GENE_ANNOT_BED_DIR)

if not os.path.exists(DRIVER_GENES_DIR):
    os.makedirs(DRIVER_GENES_DIR)

# get list of cancer types from the cancer samples
samples_df = pd.read_csv(SAMPLE_SHEET, sep="\t")
CANCER_TYPES = samples_df["Project ID"].unique().tolist()
CANCER_TYPES = [cancer.replace("TCGA-", "") for cancer in CANCER_TYPES]
CANCER_TYPES.sort()
print(f"Found {len(CANCER_TYPES)}", CANCER_TYPES)

# get dictionary of chromosome sizes
chrom_sizes = pd.read_csv(CHROM_SIZES, sep="\t", header=None, names=["chrom", "size"])
chrom_sizes["size"] = chrom_sizes["size"].astype(int)
chrom_sizes_dict = chrom_sizes.set_index("chrom")["size"].to_dict()

# gene annotation from gencode for getting gene start and end positions
gencode_annot = pd.read_csv(GENCODE_ANNOT, sep="\t", comment="#", header=None)
gencode_annot = gencode_annot[gencode_annot[2] == "gene"]
gencode_annot = gencode_annot.sort_values(by=[0, 3, 4]).reset_index(drop=True)

# refGene annotation for getting genic reigion coordinates
# get exon, UTRs, start/stop codon, transcript, CDS annotation
gene_annot = pd.read_csv(REFGENE_ANNOT, sep="\t", header=None)
# get intronic regions
intron_annot = pd.read_csv(INTRON_ANNOTATION, sep="\t", header=None)
intron_annot = intron_annot[[0, 1, 2, 5]]
intron_annot.columns = ["chrom", "start", "end", "strand"]
intron_annot["region"] = "intron"

def get_upstream_and_downstream(coords, length, prev_coord, next_coord):
    # if previous gene overlaps with the specified upstream region, start at the end of the previous gene
	if (prev_coord is not None) and (coords[1] - length <= prev_coord[2]):
		start = prev_coord[2] + 1
		# if genes are overlapping, start at original start of the gene
		if start > coords[1]:
			start = coords[1]
	else:
		start = coords[1] - length
		# if start is extending beyond the start of the chromosome, start from 0
		if start < 1:	# GTF is 1-based
			start = 1
	# if next gene overlaps with the specified downstream region, end at the start of the next gene
	if (next_coord is not None) and (coords[2] + length >= next_coord[1]):
		end = next_coord[1] - 1
		# if genes are overlapping, end at the original end of the gene
		if end < coords[2]:
			end = coords[2]
	else:
		end = coords[2] + length
		# if end is extending beyond the end of the chromosome, end at the end of the chromosome
		if end > chrom_sizes_dict[coords[0]]:
			end = chrom_sizes_dict[coords[0]]
	return (coords[0], start, end)

# merge overlapping coordinates within the same genomic region
def merge_overlapping_intervals(gene_annot):
	regions = gene_annot["region"].unique()
	gene_annot_merged = pd.DataFrame()
	for region in regions:
		gene_annot_region = gene_annot[gene_annot["region"] == region]
		gene_annot_region = gene_annot_region.sort_values(["chrom", "start", "end", "strand"])
		gene_annot_region["4"] = "."
		gene_annot_region["5"] = "."
		gene_annot_region = gene_annot_region[["chrom", "start", "end", "4", "5", "strand"]]
		annot_bed = f"data/ref/{region}.bed"
		gene_annot_region.to_csv(annot_bed, sep="\t", header=False, index=False)
		merged_bed = f"data/ref/{region}_merged.bed"
		# merge overlapping coordinates within the same strand
		cmd = f"bedtools merge -i {annot_bed} -s -c 6 -o distinct > {merged_bed}"
		os.system(cmd)
		if os.path.getsize(merged_bed) != 0:
			gene_annot_region = pd.read_csv(merged_bed, sep="\t", header=None)
			gene_annot_region.columns = ["chrom", "start", "end", "strand"]
			gene_annot_region["region"] = region
			gene_annot_merged = pd.concat([gene_annot_merged, gene_annot_region])
		os.system(f"rm {annot_bed} {merged_bed}")
	if gene_annot_merged.empty:
		return gene_annot_merged
	gene_annot_merged = gene_annot_merged.sort_values(["strand", "chrom", "start", "end"])
	return gene_annot_merged

# get genic region lengths and annotations for each driver gene
def get_genic_regions(cancer_type, driver_genes_df):
	genic_region_annotations = pd.DataFrame()
	regions = ["exon", "intron", "3UTR", "5UTR", "start_codon", "stop_codon", "transcript", "CDS"]
	for region in regions:
		driver_genes_df[f"{region}_length"] = 0
	for index, row in driver_genes_df.iterrows():
		# get exon, UTRs, start/stop codon, transcript, CDS annotation
		gene_annot_data = gene_annot[gene_annot[8].str.contains(f'gene_name "{row["gene"]}"')]
		gene_annot_data = gene_annot_data[[0, 3, 4, 6, 2]]
		gene_annot_data.columns = ["chrom", "start", "end", "strand", "region"]
		gene_annot_data["start"] = gene_annot_data["start"].astype(int) - 1 # convert to 0-based for bedtools
	
		# get intron annotation
		intron_annot_data = intron_annot[intron_annot["chrom"] == row["chr"]]
		intron_annot_data = intron_annot_data[intron_annot_data["strand"] == row["strand"]]
		intron_annot_data = intron_annot_data[intron_annot_data["start"] >= row["gene_start"]]
		intron_annot_data = intron_annot_data[intron_annot_data["end"] <= row["gene_end"]]
		gene_annot_data = pd.concat([gene_annot_data, intron_annot_data])
		gene_annot_driver_gene = merge_overlapping_intervals(gene_annot_data)
		if gene_annot_driver_gene.empty:
			continue

		# get genic region annotations for the driver gene
		gene_annot_driver_gene["gene"] = row["gene"]
		genic_region_annotations = pd.concat([genic_region_annotations, gene_annot_driver_gene])

		# get genic region lengths
		gene_annot_driver_gene["length"] = gene_annot_driver_gene["end"] - gene_annot_driver_gene["start"]
		gene_annot_driver_gene = gene_annot_driver_gene.groupby("region").agg({
			"length": "sum"
		})
		gene_annot_driver_gene = gene_annot_driver_gene.reset_index()
		for region in regions:
			length = gene_annot_driver_gene[gene_annot_driver_gene["region"] == region]["length"]
			driver_genes_df.loc[index, f"{region}_length"] = length.values[0] if length.shape[0] > 0 else 0

	# save genic region annotations to a bed file
	genic_region_annotations = genic_region_annotations.sort_values(["chrom", "start", "end", "strand"])
	genic_region_annotations["chrom"] = genic_region_annotations["chrom"].str.replace("chr", "")
	genic_region_annotations = genic_region_annotations[["chrom", "start", "end", "strand", "region", "gene"]]
	genic_region_annotations.to_csv(f"{GENE_ANNOT_BED_DIR}/{cancer_type}.bed", sep="\t", header=False, index=False)

	# drivers genes and their genic region lengths
	driver_genes_df["total_length"] = driver_genes_df["gene_end"] - driver_genes_df["gene_start"] + 1
	driver_genes_df["promoter_length"] = driver_genes_df["gene_start"] - driver_genes_df["promoter_start"]
	return driver_genes_df

for cancer_type in CANCER_TYPES:
	cancer_type = cancer_type.replace("TCGA_", "")
	files = os.listdir(DRIVER_GENES_INTOGEN)
	cancer_files = [f for f in files if f.startswith(cancer_type)]
	driver_genes_intogen = pd.DataFrame()
	for fname in cancer_files:
		driver_genes_intogen = pd.concat([driver_genes_intogen, pd.read_csv(f"{DRIVER_GENES_INTOGEN}/{fname}", sep="\t")])
	driver_genes_intogen.sort_values(by="Samples (%)", ascending=False, inplace=True)
	driver_genes_intogen.drop_duplicates(subset=["Symbol"], inplace=True, keep='first')
	driver_genes_intogen["Samples (%)"] = driver_genes_intogen["Samples (%)"].astype(float)
	driver_genes_intogen = driver_genes_intogen[driver_genes_intogen["Samples (%)"] >= 2]
	# get top 10 driver genes, include genes with the same percentage as the 10th gene
	top_percentages = driver_genes_intogen["Samples (%)"].head(10)
	driver_genes_top = driver_genes_intogen[driver_genes_intogen["Samples (%)"].isin(top_percentages)]
	driver_gene_list = driver_genes_top["Symbol"].to_list()
	# TERT gene is known to be a driver in several cancers according to PCAWG, but is not included in the intogen list
	tert_promoter_cancers = ['BLCA', 'GBM', 'HNSC', 'LIHC', 'SKCM', 'THCA']
	if "TERT" not in driver_gene_list and cancer_type in tert_promoter_cancers:
		driver_gene_list.append("TERT")

	gene_objects = []
	for driver_gene in driver_gene_list:
		idx = gencode_annot[gencode_annot[8].str.contains(f'gene_name "{driver_gene}"')]
		if idx.empty:
			print(f"Gene {driver_gene} not found in the annotation file")
			continue
		idx = idx.index[0]
		strand = gencode_annot.loc[idx][6]
		gene_of_interest = gencode_annot.loc[idx]

		# get genes on the same strand to get the previous and next gene
		gene_annot_strand = gencode_annot[gencode_annot[6] == strand].reset_index(drop=True)
		idx = gene_annot_strand[gene_annot_strand[8].str.contains(f'gene_name "{driver_gene}"')].index[0]
		gene_of_interest = gene_annot_strand.loc[idx]

		prev_gene = gene_annot_strand.loc[idx-1] if idx > 0 else None
		next_gene = gene_annot_strand.loc[idx+1] if idx < gene_annot_strand.shape[0]-1 else None
		gene_of_interest_coords = (gene_of_interest[0], gene_of_interest[3], gene_of_interest[4])
		prev_coord = (prev_gene[0], prev_gene[3], prev_gene[4]) if prev_gene is not None else None
		next_coord = (next_gene[0], next_gene[3], next_gene[4]) if next_gene is not None else None

		# coordinates of the driver genes including 10kB upstream and downstream of the gene, and promoter region
		extended_coords = get_upstream_and_downstream(gene_of_interest_coords, 10000, prev_coord, next_coord)
		promoter_coords = get_upstream_and_downstream(gene_of_interest_coords, 2500, prev_coord, next_coord)
		gene_coord_obj = {
			"gene": driver_gene,
			"chr": gene_of_interest_coords[0],
			"strand": gene_of_interest[6],
			"start": extended_coords[1], # start of 10kb upstream of the gene
			"end": extended_coords[2], # end of 10kb downstream of the gene
			"gene_start": gene_of_interest_coords[1], # start of the gene
			"gene_end": gene_of_interest_coords[2], # end of the gene
			"gene_length": gene_of_interest_coords[2] - gene_of_interest_coords[1] + 1, # length of the gene
			"promoter_start": promoter_coords[1] - 1, # start of the promoter
			"promoter_end": gene_of_interest_coords[1] - 1, # promoter ends right before the start of the gene
		}
		gene_objects.append(gene_coord_obj)
	genes_df = pd.DataFrame(gene_objects)
	print(f"{cancer_type}: {len(genes_df)} driver genes")

	driver_genes_list = genes_df["gene"].tolist()
	print(driver_genes_list)
	driver_genes_df = get_genic_regions(cancer_type, genes_df)
	driver_genes_df.head()
	driver_genes_df.to_csv(f"{DRIVER_GENES_DIR}/{cancer_type}.tsv", sep="\t", index=False)

# get list of all driver genes
all_driver_genes = []
for cancer_type in CANCER_TYPES:
	cancer_type = cancer_type.replace("TCGA_", "")
	driver_genes_df = pd.read_csv(f"{DRIVER_GENES_DIR}/{cancer_type}.tsv", sep="\t")
	all_driver_genes.extend(driver_genes_df["gene"].tolist())
all_driver_genes = list(set(all_driver_genes))
print(len(all_driver_genes))
all_driver_genes_df = pd.DataFrame(all_driver_genes, columns=["gene"])
all_driver_genes_df.to_csv(GENE_LIST, sep="\t", index=False)
