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
    *   Used for lifting over variant coordinates from hg38 to the haplotype assemblies (required for `haplotype_full_sim.py`).
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

4.  **Custom Haplotype Assemblies** (Required for `haplotype_full_sim.py`)
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

> **Algorithm details:** see [`ALGORITHM.md`](ALGORITHM.md) for a full description of
> the simulation approach, coordinate liftover, strand handling, and resource estimates.

### 1. Full Diploid Simulation — canonical workflow (`scripts/haplotype_full_sim.py`)

This is the **recommended** script. For each disease it:

1. Lifts over the hg38 variant coordinate to each haplotype assembly using minimap2.
2. Patches the **entire** Hap1 and/or Hap2 assembly (dominant → Hap1 only; recessive → both).
3. Concatenates the patched assemblies into a diploid FASTA.
4. Simulates paired-end reads with `art_illumina` from the **full** diploid genome.

Reads originate from the individual's real genome, preserving natural background variation
as a realistic context for variant calling evaluation.

**Run all diseases (default 15× per haplotype, ≈ 30× diploid):**
```bash
python scripts/haplotype_full_sim.py
```

**Run a single disease at low coverage for a quick test:**
```bash
python scripts/haplotype_full_sim.py --disease "sickle cell" --coverage 1
```

**All options:**
```
--disease NAME      Run only this disease (case-insensitive substring match)
--coverage INT      Coverage per haplotype (default: 15)
--read-length INT   Read length in bp (default: 150)
--output-dir DIR    Override output root directory
--hap1 FASTA        Override Hap1 assembly path
--hap2 FASTA        Override Hap2 assembly path
--hg38 FASTA        Override hg38 path
```

*   Requires `hg38.fa`, the GCA haplotype assemblies in `data/`, and `minimap2` in PATH.
*   Output per disease in `output/<disease_name>/`.

### 2. Simple hg38 Batch Simulation (`scripts/simulate_batch.py`)

Patches variants directly into hg38 and simulates reads from hg38. Useful for quick
checks that do not require the individual's diploid assemblies.

```bash
python scripts/simulate_batch.py
```

*   Reads `config/diseases.json`; outputs to `output/<disease_name>/`.
*   Does **not** preserve the individual's background variation.

### 3. Updating Variants from ClinVar (`scripts/update_config_from_clinvar.py`)

Looks up OMIM IDs in the local ClinVar summary and updates variant coordinates in
`config/diseases.json`.

```bash
python scripts/update_config_from_clinvar.py
```

*   Requires `data/variant_summary.txt.gz`.

### Deprecated scripts

| Script | Reason deprecated |
|--------|------------------|
| `scripts/map_and_simulate.py` | Simulates reads from a 20 kb ROI only — reads cannot seed correctly against a full-genome aligner or pangenome |
| `scripts/simulate_master_genome.py` | Complex background+locus approach, superseded by `haplotype_full_sim.py` |

Do not use the deprecated scripts for alignment benchmarks.


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

```
output/
└── <disease_name>/
    ├── hap1_patched.fa                  # Full Hap1 assembly with variant applied
    ├── hap2_patched.fa                  # Full Hap2 assembly with variant applied (recessive)
    ├── diploid_patched.fa               # Concatenated diploid FASTA (input to ART)
    ├── <disease_name>_1.fq              # Simulated R1 reads (paired-end)
    ├── <disease_name>_2.fq              # Simulated R2 reads
    └── <disease_name>.phenopacket.json  # GA4GH Phenopacket with HPO terms + variant
```

## Tools Overview

*   `scripts/haplotype_full_sim.py`: **Canonical simulation script** — full diploid workflow.
*   `src/genome_patcher.py`: Patches a variant into a full FASTA assembly (exact contig match, 60-char line wrapping).
*   `src/read_simulator.py`: Wrapper around `art_illumina`; falls back to a pure-Python error-free simulator.
*   `src/phenotype.py`: Generates GA4GH Phenopackets from OMIM IDs and local `phenotype.hpoa`.
*   `src/main.py`: CLI entry point for single-sample simulation (used by `simulate_batch.py`).