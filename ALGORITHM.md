# Project Nelly — Algorithm

## Goal

Generate synthetic paired-end Illumina reads for a **real diploid individual** carrying a known
pathogenic variant, then evaluate whether a downstream variant calling pipeline recovers that
variant correctly.

The reads must originate from the individual's own genome (not from hg38) so that their natural
background variation is preserved — this is what makes the benchmark realistic.

---

## Data model

| File | Description |
|------|-------------|
| `hg38.fa` | GRCh38 reference genome — used only as a coordinate source |
| `data/GCA_050492415.1_apr041.1_v1_genomic.fna` | Haplotype 1 (Hap1) assembly of the target individual |
| `data/GCA_050492395.1_apr041.2_v1_genomic.fna` | Haplotype 2 (Hap2) assembly of the target individual |
| `config/diseases.json` | Disease definitions: variant in hg38 coordinates, inheritance mode |

The haplotype assemblies are de-novo assemblies (e.g. from PacBio/ONT) for a specific individual
(here: a Saudi individual, accessions JBHMRZ/JBHMSA). Their contigs do **not** share chromosome
names with hg38; coordinates must be lifted over before patching.

---

## The canonical algorithm — `scripts/haplotype_full_sim.py`

### Step 1 — Coordinate liftover (hg38 → assembly)

The disease variant is stored in hg38 coordinates (`chrom`, `pos`, `ref`, `alt`).
To find the equivalent position in each haplotype assembly:

1. Extract a ±10 kb genomic context window from hg38 around the variant using `pysam`.
2. Write this context as a temporary FASTA and map it to each assembly with
   `minimap2 -x asm5 --secondary=no`.
3. Parse the PAF output to find which assembly contig the context mapped to, the target
   start/end, and the strand.
4. Compute the variant's position within the aligned region using the query offset:
   - Forward strand: `target_pos = t_start + (variant_offset − q_start) + 1`
   - Reverse strand: `target_pos = t_end − (variant_offset − q_start)`  (1-based)

> **Why hg38 context and not a cDNA sequence?**
> Using a genomic context window is equivalent to a gene-level search for SNPs and small indels:
> the ±10 kb window spans the gene locus and minimap2 finds it in the assembly regardless of
> contig naming. A cDNA/exon-only query could be used for better specificity in multi-exon genes
> or repeat-rich loci, but requires an additional exon-coordinates lookup step and is not yet
> implemented.

### Step 2 — Patch the full haplotype assemblies

Using `src/genome_patcher.py`, apply the variant to the **entire** assembly FASTA (not a sub-region):

- **Dominant (heterozygous):** patch Hap1 only; Hap2 remains wildtype.
- **Recessive (homozygous):** patch both Hap1 and Hap2.

When the minimap2 alignment is on the reverse strand the ref/alt alleles are
reverse-complemented before patching, because the assembly FASTA is always in forward
orientation.

`genome_patcher.py` streams the FASTA contig-by-contig, so memory usage is bounded by the
largest single contig, not the entire assembly size. Output is written with 60-character
line wrapping (required by `samtools faidx`).

### Step 3 — Build a diploid FASTA

Concatenate patched-Hap1 and (patched or original) Hap2 into a single diploid FASTA.
This is what ART will read.

### Step 4 — Simulate paired-end reads from the full diploid

Run `art_illumina` on the full diploid FASTA at the requested per-haplotype coverage
(default 15×; diploid total ≈ 30×).

**This is the critical step that previous scripts got wrong.** Earlier scripts
(`map_and_simulate.py`, `simulate_master_genome.py`) simulated reads from a 20 kb ROI only.
Reads from a 20 kb FASTA cannot seed correctly against a pangenome/hg38 aligner — the aligner
exhausts its search space on every read, causing indefinite runtime.

By simulating from the full diploid assembly:
- Reads cover the entire genome at uniform depth.
- The disease locus is covered at the same depth as the rest.
- The individual's natural variation (SNPs, indels relative to hg38) is present throughout
  the reads, making the benchmark realistic.

### Step 5 — Phenopacket

A GA4GH Phenopacket JSON is written per sample, encoding HPO phenotype terms, the variant
(as a `VcfRecord` in hg38 coordinates), zygosity, and ACMG pathogenicity classification.

---

## Evaluation workflow (external to this tool)

After simulation, the intended downstream evaluation is:

```
FASTQ  →  align (vg giraffe / BWA-MEM2)  →  call variants  →  compare to truth
```

**Truth set:** the variant inserted in Step 2 (stored in `config/diseases.json`).

**Comparison:** check whether the called variants include the expected `(chrom, pos, ref, alt)`
at the expected zygosity (het for dominant, hom for recessive). A simple VCF intersection
(`bcftools isec` or `rtg vcfeval`) is sufficient.

This evaluation step is not yet implemented in Project Nelly itself; it is expected to be run
as part of the calling pipeline being tested.

---

## Script inventory

| Script | Status | Purpose |
|--------|--------|---------|
| `scripts/haplotype_full_sim.py` | **Canonical** | Full-genome diploid sim (this algorithm) |
| `scripts/simulate_batch.py` | Utility | Simple hg38-patching; useful for quick checks |
| `scripts/map_and_simulate.py` | Deprecated | ROI-only sim — **do not use for alignment benchmarks** |
| `scripts/simulate_master_genome.py` | Deprecated | Background+locus approach — superseded above |
| `scripts/update_config_from_clinvar.py` | Utility | Auto-update variant coords from ClinVar |
| `scripts/regenerate_phenopackets.py` | Utility | Re-generate all phenopackets from config |
| `scripts/update_vcf.py` | Utility | Export config to VCF format |

---

## Resource requirements

Patching and simulating from a full assembly (~3 GB per haplotype) is I/O intensive.

| Step | Approximate time | Approximate disk |
|------|-----------------|-----------------|
| minimap2 liftover | < 1 min | negligible |
| Patch full Hap1 | 5–10 min | ~3 GB |
| Patch full Hap2 | 5–10 min | ~3 GB |
| Diploid concat | 2–5 min | ~6 GB |
| ART 15× on 6 GB | 30–90 min | ~20 GB FASTQ |

Per disease. Running all 10 diseases requires ~300 GB disk and ~10–20 h on a single machine.

For rapid testing, use `--coverage 1` and `--disease "sickle cell"` to run a single disease
at minimal depth.

---

## Known limitations

- Indel support in `genome_patcher.py` is experimental; SNPs are fully validated.
- Trinucleotide repeat expansions (e.g. Huntington) require the assembly to have enough
  length at the insertion site; the patcher applies a simple string replacement.
- The liftover uses genomic context (not cDNA); very long introns or segmental duplications
  near a variant could produce ambiguous minimap2 hits.
- The evaluation step (variant calling + comparison) must be run externally.
