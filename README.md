# Continuum model of Cancer development

This repository hosts the codebase and analytical workflows for the manuscript “The Cumulative Impact of Passenger Mutations on Cancer Development.” The project explores the functional effects of passenger mutations across diverse cancer types, quantifies their impact on patient survival, and evaluates their value as prognostic biomarkers.

The analyses performed include:

- Quantification of passenger mutation burden in known cancer genes
- Prediction of pathogenicity for passenger mutations using computational tools
- Assessment of functional impact on transcription factor binding, splicing and gene expression
- Evaluation of survival outcomes and the prognostic potential of passenger mutations

## Repository Structure

- **data/**: Raw data used for all analyses are stored here. Due to privacy and data access restrictions, only a subset of the data is shared via Zenodo (DOI: 10.5281/zenodo.18224257). Users are responsible for independently obtaining access to controlled datasets, including TCGA, in accordance with the relevant data use agreements.
- **metadata/**: Contains text files and tabular data obtained from multiple sources, and used for analyses.
- **notebooks/**: Contains Jupyter notebooks that use the processed data in the **results/** directory for further analysis and data visualization.
- **scripts/**: Contains scripts used for downloading, parsing, and pre-processing data. Also contains meta scripts for running other scripts.
- **POG570_scripts/**: Scripts and Notebooks used for reproducing core analyses in the POG570 dataset.
- **TCGA_scripts/**: Scripts and Notebooks used for reproducing core analyses in the entire TCGA WGS dataset.
- **plot_data/**: Source data for all plots and figures.
- **conda_env_cancer_model.yml**: Conda environment specification file, used for reproducing the conda environment used here.

## Getting Started

To get started with the project, clone the repository and install the required dependencies. Create a new conda environment for the project using the provided specifications.

```bash

# clone the repository
git clone https://github.com/Georgakopoulos-Soares-lab/cancer_passengers.git
cd cancer_passengers

# create a new conda environment
conda env create -f conda_env_cancer_passengers.yml
conda activate cancer_passengers

```

Configure the following at the root level of the directory:

- Download and setup [Annovar](https://www.openbioinformatics.org/annovar/annovar_download_form.php)
- Download CADD [file](https://cadd.gs.washington.edu/download)
- Setup [Pangolin tool](https://github.com/tkzeng/Pangolin)

## Datasets Required
- The PCAWG data are accessible through the ICGC-ARGO, and can be downloaded by following the instructions provided in the ICGC-25K data access documentation (https://docs.icgc-argo.org/docs/data-access/icgc-25k-data). 
- The TCGA portion of the dataset is controlled-access and is available upon dbGaP authorization under accession phs000178.v11.p8. 
- Survival data were collected from the TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR).
- Driver mutations were obtained from the following data files from the ICGC/TCGA driver mutation pan-cancer compendium -
- - TableS3_panorama_driver_mutations_ICGC_samples.public.tsv.gz (public)
- - TableS3_panorama_driver_mutations_TCGA_samples.controlled.tsv.gz (controlled access)
- POG570 data was obtained from https://www.bcgsc.ca/downloads/POG570.
- Cancer Hotspot Mutations (v2) were downloaded from https://www.cancerhotspots.org/#/download.

## Contact
For questions or to report bugs, contact -
```
abn5461@psu.edu (Akshatha Nayak)
izg5139@psu.edu (Ilias Georgakopoulos-Soares)
```
