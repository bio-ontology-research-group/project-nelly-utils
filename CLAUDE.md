# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Install dependencies:**
```bash
pip install -r requirements.txt
```

**Run all tests:**
```bash
pytest
```

**Run a single test file:**
```bash
pytest tests/test_genome_patcher.py
```

**Run a single test:**
```bash
pytest tests/test_genome_patcher.py::test_function_name
```

**System tools required** (external, not pip-installable):
- `art_illumina` — Illumina read simulator (primary)
- `minimap2` — sequence mapping (for haplotype-aware workflows)
- `samtools` — optional, for BAM/FASTA indexing

## Architecture

Project Nelly generates synthetic NGS reads and GA4GH Phenopackets for genetic disease research. Given a disease config with variant coordinates, it patches a reference genome, simulates Illumina reads, and produces standardized phenotype metadata.

### Core modules (`src/`)

- **`genome_patcher.py`** — applies variants (SNPs/indels) to FASTA files using `pysam`; supports haplotype-specific patching with ref allele verification
- **`read_simulator.py`** — wraps `art_illumina` to produce paired-end FASTQs; falls back to a pure-Python simulator if `art_illumina` is unavailable
- **`phenotype.py`** — builds GA4GH Phenopackets from local `data/phenotype.hpoa`; encodes inheritance pattern (heterozygous for dominant, homozygous for recessive)
- **`main.py`** — CLI entry point for single-sample simulation; orchestrates the three modules above

### Simulation workflows (4 tiers, `scripts/`)

1. **`simulate_batch.py`** — simplest; patches hg38 directly and simulates reads for each disease in `config/diseases.json`
2. **`map_and_simulate.py`** — haplotype-aware; maps variant context from hg38 to diploid assemblies (Hap1/Hap2) using `minimap2`, respects inheritance patterns
3. **`simulate_master_genome.py`** — most efficient for large batches; simulates a single shared background genome once, then only simulates disease-specific loci (~20 kb ROI) per disease; outputs an `assemble_samples.sh` script that links background + locus FASTQs into per-patient reads
4. **`update_config_from_clinvar.py`** — utility to auto-update `config/diseases.json` variant coordinates from the local ClinVar database

### Disease configuration (`config/diseases.json`)

Each disease entry includes `name`, `omim_id`, `inheritance` (`dominant`|`recessive`), and a `variant` object with `chrom`, `pos`, `ref`, `alt`. This is the primary place to add or modify simulated diseases.

### Required data files (not in repo)

- `hg38.fa` — HG38 reference genome (place in project root)
- `data/phenotype.hpoa` — HPO annotations
- `data/variant_summary.txt.gz` — ClinVar data (optional; only needed for `update_config_from_clinvar.py`)
- `data/GCA_050492415.1_apr041.1_v1_genomic.fna` — Haplotype 1 assembly (advanced workflows)
- `data/GCA_050492395.1_apr041.2_v1_genomic.fna` — Haplotype 2 assembly (advanced workflows)

### Known limitations

- `genome_patcher.py` has partial indel support; SNPs are fully supported
- `phenotype.py` uses local `phenotype.hpoa` file lookup, not a live API
