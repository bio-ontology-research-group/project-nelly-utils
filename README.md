# Project Nelly Utils

This repository contains utilities for generating synthetic genomic data with pathogenic variants, specifically designed for Project Nelly. It allows for the creation of simulated NGS reads (FASTQ) and Phenopackets for specific genetic diseases.

## Prerequisites

### Software
1.  **Python 3.8+**
2.  **ART Read Simulator** (`art_illumina`)
    *   Used for simulating Illumina sequencing reads.
    *   Installation (Conda): `conda install bioconda::art`
    *   Installation (Ubuntu/Debian): `sudo apt-get install art-nextgen-simulation-tools`
3.  **Minimap2**
    *   Used for mapping variants from HG38 to custom assemblies (required for `map_and_simulate.py`).
    *   Installation (Conda): `conda install bioconda::minimap2`
    *   Installation (Ubuntu/Debian): `sudo apt-get install minimap2`

### Python Dependencies
Install the required Python packages:

```bash
pip install -r requirements.txt
```

## Data Setup

Before running the simulations, you need to populate the `data/` directory with the necessary reference genomes.

1.  **HG38 Reference Genome**
    *   Download `hg38.fa` and place it in the **root directory** of this project.
    *   **Download Link:** [UCSC hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
    *   After downloading, unzip it: `gunzip hg38.fa.gz`
    *   Index it: `samtools faidx hg38.fa` (optional but recommended for speed).

2.  **Human Phenotype Ontology (HPO) Data**
    *   Required for generating Phenopackets.
    *   Download `phenotype.hpoa` from the HPO project.
    *   URL: `http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa`
    *   Place it in `data/phenotype.hpoa`.

3.  **ClinVar Data** (Optional, for auto-updating variants)
    *   Required if you want to update `config/diseases.json` with real pathogenic variants.
    *   Download `variant_summary.txt.gz` from NCBI ClinVar.
    *   URL: `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz`
    *   Place it in `data/variant_summary.txt.gz`.

4.  **Custom Haplotype Assemblies** (Required for `map_and_simulate.py`)
    *   Download the following assemblies and place them in the `data/` directory. These are from the NCBI **apr041** dataset.
    *   **Using NCBI Datasets CLI:**
        ```bash
        # Haplotype 1 (GCA_050492415.1)
        datasets download genome accession GCA_050492415.1 --include genome
        unzip ncbi_dataset.zip -d data_hap1
        mv data_hap1/ncbi_dataset/data/GCA_050492415.1/GCA_050492415.1_apr041.1_v1_genomic.fna data/

        # Haplotype 2 (GCA_050492395.1)
        datasets download genome accession GCA_050492395.1 --include genome
        unzip ncbi_dataset.zip -d data_hap2
        mv data_hap2/ncbi_dataset/data/GCA_050492395.1/GCA_050492395.1_apr041.2_v1_genomic.fna data/
        ```
    *   The files should be named:
        *   `data/GCA_050492415.1_apr041.1_v1_genomic.fna`
        *   `data/GCA_050492395.1_apr041.2_v1_genomic.fna`

## Usage

The project provides two main workflows for generating data.

### 1. Simple Batch Simulation (`scripts/simulate_batch.py`)
This script uses `hg38.fa` as the template for all simulations. It patches the specific variants directly into hg38 and simulates reads.

**Command:**
```bash
python scripts/simulate_batch.py
```
*   Reads configuration from `config/diseases.json`.
*   Generates output in `output/<disease_name>/`.

### 2. Updating Variants from ClinVar (`scripts/update_config_from_clinvar.py`)
This script scans the `config/diseases.json` file, looks up the OMIM IDs in the downloaded ClinVar summary, and automatically updates the variant information (chromosome, position, ref, alt) with a known pathogenic variant.

**Command:**
```bash
python scripts/update_config_from_clinvar.py
```
*   Requires `data/variant_summary.txt.gz`.
*   Useful for ensuring your config has valid pathogenic variants for the specified diseases.

### 3. Advanced Haplotype-Aware Simulation (`scripts/map_and_simulate.py`)
This is the **recommended** method for more realistic data. It maps the context of the variant from HG38 to the custom diploid assemblies (Hap1 and Hap2), patches the correct haplotype based on inheritance patterns (dominant vs. recessive), and simulates reads from the specific Region of Interest (ROI).

**Command:**
```bash
python scripts/map_and_simulate.py
```
*   Requires the GCA assemblies in `data/`.
*   Extracts only the relevant genomic regions (approx. 20kb around the variant) for simulation to save time and space.

## Configuration

The diseases and variants to be simulated are defined in `config/diseases.json`. You can add more diseases by following the existing format:

```json
{
    "name": "Disease Name",
    "omim_id": "123456",
    "inheritance": "recessive",
    "variant": {
        "chrom": "chr1",
        "pos": 1234567,
        "ref": "A",
        "alt": "T"
    }
}
```

## Output Structure

The output is organized by disease name in the `output/` directory:

```
output/
└── <disease_name>/
    ├── <disease_name>.phenopacket.json  # Standardized phenotype description
    ├── <disease_name>1.fq               # Simulated Read 1 (FASTQ)
    ├── <disease_name>2.fq               # Simulated Read 2 (FASTQ)
    ├── roi_diploid.fa                   # (map_and_simulate only) The specific genomic region used
    └── ...
```

## Tools Overview

*   `src/main.py`: The core driver script for single-sample simulation.
*   `src/genome_patcher.py`: Handles checking reference alleles and applying mutations to FASTA sequences.
*   `src/read_simulator.py`: A wrapper around `art_illumina`.
*   `src/phenotype.py`: Generates GA4GH Phenopackets based on OMIM IDs.