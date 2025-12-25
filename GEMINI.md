# Gemini Customization Guide for Project Nelly

This document provides guidance on how to customize and extend Project Nelly, particularly for a Gemini agent or a developer working with one.

## Project Structure

The project is organized as follows:

-   `src/`: Contains the core Python source code.
    -   `main.py`: The main CLI entry point.
    -   `genome_patcher.py`: Handles modifying the genome FASTA file.
    -   `read_simulator.py`: A wrapper for the ART read simulator.
    -   `phenotype.py`: Creates Phenopackets.
-   `data/`: Intended for input data like reference genomes and ClinVar VCFs.
-   `output/`: The default directory for all generated files.
-   `tests/`: Contains the test suite.
-   `README.md`: General user documentation.
-   `requirements.txt`: Python dependencies.
-   `GEMINI.md`: This file.

## Customization and Extension

### 1. Supporting Indels in `genome_patcher`

The `genome_patcher.py` module currently only supports SNPs. To add support for insertions and deletions (indels):

1.  **Modify `genome_patcher.py`:**
    -   Locate the `if len(variant.ref) == 1 and len(variant.alts[0]) == 1:` block.
    -   Add `elif` or `else` conditions to handle cases where the lengths of the reference and alternate alleles are not equal.
    -   For an **insertion**, the reference allele at the position needs to be prepended with the new sequence from the alternate allele.
    -   For a **deletion**, the reference allele's sequence needs to be removed from the genome.
    -   Be mindful of VCF format conventions for indels.

2.  **Add tests:**
    -   Create a new test case in `tests/test_genome_patcher.py`.
    -   Add a variant representing an indel to `tests/data/test_clinvar.vcf`.
    -   Assert that the output FASTA file reflects the correct indel.

### 2. Implementing Real HPO Term Fetching

The `phenotype.py` module uses a placeholder function `get_hpo_terms_for_disease`. To implement a real search:

1.  **Choose a data source:** The HPO project provides its data for download in various formats (e.g., OBO). You can download these files and parse them.
2.  **Implement the search logic:**
    -   Parse the HPO data file to build a mapping between diseases and HPO terms.
    -   In `get_hpo_terms_for_disease`, use this mapping to find the HPO terms for the given disease name.
3.  **Consider APIs:** While the HPO project doesn't have a simple, official REST API for this purpose, other services might. Using an external API would require the `requests` library.

### 3. Adding a New Read Simulator

To support a different read simulator (e.g., `wgsim`):

1.  **Create a new function in `read_simulator.py`:** For example, `simulate_reads_wgsim(...)`.
2.  **Implement the command-line wrapper:** This function should construct and execute the `wgsim` command using `subprocess`.
3.  **Update `main.py`:** Add a command-line argument to allow the user to choose the simulator. Based on this choice, call the appropriate function in `read_simulator.py`.

## Running Tests

To run the test suite, first install the development dependencies:

```bash
pip install -r requirements.txt
```

Then run `pytest`:

```bash
pytest
```
