#!/usr/bin/env python3
"""
haplotype_full_sim.py — Full diploid genome simulation for variant calling evaluation.

This is the canonical simulation script for Project Nelly.

Algorithm
---------
For each disease in config/diseases.json:

  Step 1 — Coordinate liftover
    Extract a ±10 kb genomic context from hg38 around the variant position.
    Map this context to each haplotype assembly using minimap2 (asm5 preset).
    Use the alignment offset to compute the exact variant position in each assembly.

  Step 2 — Patch full assemblies
    Apply the variant to the FULL haplotype assembly FASTA using genome_patcher.py.
    Dominant (heterozygous): patch Hap1 only; Hap2 is wildtype.
    Recessive (homozygous):  patch both Hap1 and Hap2.
    Strand is respected: if minimap2 maps on '-', ref/alt are reverse-complemented
    before patching the forward-strand assembly sequence.

  Step 3 — Build diploid FASTA
    Concatenate patched Hap1 + (patched or original) Hap2.

  Step 4 — Simulate paired-end reads
    Run art_illumina on the full diploid FASTA.
    Reads originate from the individual's own genome, preserving their
    natural variation as a realistic background for variant calling.

  Step 5 — Phenopacket
    Write a GA4GH Phenopacket with HPO terms and variant interpretation.

Outputs per disease (in output/<disease_name>/):
  hap1_patched.fa         — patched Hap1 assembly
  hap2_patched.fa         — patched Hap2 assembly (recessive only)
  diploid_patched.fa      — concatenated diploid FASTA
  <sample_id>_1.fq        — R1 reads (ART paired-end)
  <sample_id>_2.fq        — R2 reads
  <sample_id>.phenopacket.json

Usage
-----
  python scripts/haplotype_full_sim.py [options]

  --disease NAME     Run only one disease (partial name match, case-insensitive)
  --coverage INT     Coverage per haplotype (default: 15; diploid total = 2×)
  --read-length INT  Read length (default: 150)
  --output-dir DIR   Override output directory
  --hap1 PATH        Override Hap1 FASTA path
  --hap2 PATH        Override Hap2 FASTA path
  --hg38 PATH        Override hg38 FASTA path
"""

import os
import sys
import json
import argparse
import shutil
import subprocess
import tempfile

import pysam

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(base_dir, "src"))

from genome_patcher import patch_genome
from read_simulator import simulate_reads
from phenotype import create_phenopacket

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_HG38   = os.path.join(base_dir, "hg38.fa")
DEFAULT_HAP1   = os.path.join(base_dir, "data", "GCA_050492415.1_apr041.1_v1_genomic.fna")
DEFAULT_HAP2   = os.path.join(base_dir, "data", "GCA_050492395.1_apr041.2_v1_genomic.fna")
DEFAULT_CONFIG = os.path.join(base_dir, "config", "diseases.json")
DEFAULT_OUTDIR = os.path.join(base_dir, "output")
CONTEXT_WINDOW = 10_000  # bp either side of variant for minimap2 query

_RC = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def reverse_complement(seq: str) -> str:
    return seq.translate(_RC)[::-1]


# ---------------------------------------------------------------------------
# Step 1 helpers: coordinate liftover
# ---------------------------------------------------------------------------

def fetch_hg38_context(hg38_fasta: str, chrom: str, pos: int, ref: str, window: int = CONTEXT_WINDOW):
    """
    Return (context_sequence, variant_offset_in_context) from hg38.
    pos is 1-based.
    """
    with pysam.FastaFile(hg38_fasta) as fa:
        refs = set(fa.references)
        if chrom not in refs:
            alt = ("chr" + chrom) if not chrom.startswith("chr") else chrom[3:]
            if alt in refs:
                chrom = alt
            else:
                raise ValueError(f"Chromosome '{chrom}' not found in hg38")
        start = max(0, pos - 1 - window)
        end   = pos - 1 + len(ref) + window
        seq   = fa.fetch(chrom, start, end)
    return seq, (pos - 1 - start)


def map_context_to_assembly(context_seq: str, assembly_fasta: str, work_dir: str):
    """
    Map context_seq to assembly using minimap2 asm5.
    Returns the best-hit dict or None.
    """
    query_fa = os.path.join(work_dir, "minimap2_query.fa")
    with open(query_fa, "w") as f:
        f.write(">query\n" + context_seq + "\n")

    cmd = ["minimap2", "-x", "asm5", "--secondary=no", assembly_fasta, query_fa]
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 9 and parts[0] == "query":
            return {
                "contig":  parts[5],
                "t_start": int(parts[7]),
                "t_end":   int(parts[8]),
                "q_start": int(parts[2]),
                "q_end":   int(parts[3]),
                "strand":  parts[4],
            }
    return None


def liftover_position(hit: dict, variant_offset: int):
    """
    Convert variant_offset (0-based, within query context) to a 1-based
    position on the assembly contig, accounting for strand.
    Returns None if the variant offset falls outside the aligned region.
    """
    if hit is None:
        return None
    if not (hit["q_start"] <= variant_offset < hit["q_end"]):
        return None
    delta = variant_offset - hit["q_start"]
    if hit["strand"] == "+":
        return hit["t_start"] + delta + 1          # 1-based
    else:
        return hit["t_end"] - delta                # 1-based (reverse strand)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_disease(d: dict, args):
    sample_id = d["name"].lower().replace(" ", "_").replace("'", "")
    print(f"\n{'='*60}")
    print(f"Disease : {d['name']}  (OMIM {d['omim_id']}, {d['inheritance']})")
    print(f"Sample  : {sample_id}")

    var = d["variant"]
    chrom, pos, ref, alt = var["chrom"], var["pos"], var["ref"], var["alt"]

    sample_dir = os.path.join(args.output_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1 — Coordinate liftover via minimap2
    # ------------------------------------------------------------------
    print("\n[1/4] Lifting over coordinates via minimap2 ...")
    try:
        context_seq, variant_offset = fetch_hg38_context(args.hg38, chrom, pos, ref)
    except ValueError as e:
        print(f"  ERROR: {e}. Skipping.")
        return False

    hit1 = map_context_to_assembly(context_seq, args.hap1, sample_dir)
    hit2 = map_context_to_assembly(context_seq, args.hap2, sample_dir)

    pos1 = liftover_position(hit1, variant_offset)
    pos2 = liftover_position(hit2, variant_offset)

    if pos1 is None and pos2 is None:
        print("  ERROR: could not lift over to either haplotype. Skipping.")
        return False

    if pos1:
        print(f"  Hap1 → {hit1['contig']}:{pos1}  strand={hit1['strand']}")
    else:
        print("  Hap1 → no hit")
    if pos2:
        print(f"  Hap2 → {hit2['contig']}:{pos2}  strand={hit2['strand']}")
    else:
        print("  Hap2 → no hit")

    # ------------------------------------------------------------------
    # Step 2 — Patch FULL haplotype assemblies
    # ------------------------------------------------------------------
    print("\n[2/4] Patching full assemblies ...")

    # Dominant = heterozygous → patch Hap1 only
    # Recessive = homozygous  → patch both
    do_patch_hap1 = (pos1 is not None)                                    # always patch if mappable
    do_patch_hap2 = (pos2 is not None) and (d["inheritance"] == "recessive")

    hap1_out = os.path.join(sample_dir, "hap1_patched.fa")
    hap2_out = os.path.join(sample_dir, "hap2_patched.fa")

    if do_patch_hap1:
        # Reverse-complement alleles when gene is on reverse strand
        p_ref1 = reverse_complement(ref) if hit1["strand"] == "-" else ref
        p_alt1 = reverse_complement(alt) if hit1["strand"] == "-" else alt
        print(f"  Patching Hap1: {hit1['contig']}:{pos1}  {p_ref1}→{p_alt1}")
        patch_genome(args.hap1, hit1["contig"], pos1, p_ref1, p_alt1, hap1_out)
    else:
        print("  Hap1: no patch applied (no hit)")
        shutil.copy2(args.hap1, hap1_out)

    if do_patch_hap2:
        p_ref2 = reverse_complement(ref) if hit2["strand"] == "-" else ref
        p_alt2 = reverse_complement(alt) if hit2["strand"] == "-" else alt
        print(f"  Patching Hap2: {hit2['contig']}:{pos2}  {p_ref2}→{p_alt2}")
        patch_genome(args.hap2, hit2["contig"], pos2, p_ref2, p_alt2, hap2_out)
    else:
        label = "dominant — Hap2 is wildtype" if d["inheritance"] == "dominant" else "no hit"
        print(f"  Hap2: no patch ({label})")
        # For the diploid FASTA we use the original directly (avoids copying several GB)
        hap2_out = args.hap2   # point to original

    # ------------------------------------------------------------------
    # Step 3 — Build diploid FASTA
    # ------------------------------------------------------------------
    print("\n[3/4] Concatenating diploid FASTA ...")
    diploid_fa = os.path.join(sample_dir, "diploid_patched.fa")
    with open(diploid_fa, "wb") as out_f:
        for src in (hap1_out, hap2_out):
            with open(src, "rb") as in_f:
                shutil.copyfileobj(in_f, out_f)
    print(f"  → {diploid_fa}")

    # ------------------------------------------------------------------
    # Step 4 — Simulate paired-end reads from full diploid
    # ------------------------------------------------------------------
    print(f"\n[4/4] Simulating reads ({args.coverage}× per haplotype) ...")
    read_prefix = os.path.join(sample_dir, sample_id + "_")
    simulate_reads(
        fasta_path=diploid_fa,
        output_prefix=read_prefix,
        coverage=args.coverage,
        read_length=args.read_length,
    )
    print(f"  → {read_prefix}1.fq  /  {read_prefix}2.fq")

    # ------------------------------------------------------------------
    # Phenopacket
    # ------------------------------------------------------------------
    pp_path = os.path.join(sample_dir, f"{sample_id}.phenopacket.json")
    create_phenopacket(
        omim_id=d["omim_id"],
        sample_id=sample_id,
        output_phenopacket_path=pp_path,
        variant_info=var,
        inheritance=d["inheritance"],
    )

    print(f"\n  Done. Outputs in: {sample_dir}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Full diploid genome simulation for variant calling evaluation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--disease",    metavar="NAME",
                        help="Run only this disease (case-insensitive substring match).")
    parser.add_argument("--coverage",   type=int, default=15, metavar="INT",
                        help="Coverage per haplotype (default: 15; diploid total ≈ 30×).")
    parser.add_argument("--read-length",type=int, default=150, metavar="INT",
                        help="Read length in bp (default: 150).")
    parser.add_argument("--output-dir", default=DEFAULT_OUTDIR, metavar="DIR",
                        help=f"Output root directory (default: {DEFAULT_OUTDIR}).")
    parser.add_argument("--hap1",       default=DEFAULT_HAP1,  metavar="FASTA",
                        help="Haplotype 1 assembly FASTA.")
    parser.add_argument("--hap2",       default=DEFAULT_HAP2,  metavar="FASTA",
                        help="Haplotype 2 assembly FASTA.")
    parser.add_argument("--hg38",       default=DEFAULT_HG38,  metavar="FASTA",
                        help="hg38 reference FASTA (for coordinate liftover query).")
    parser.add_argument("--config",     default=DEFAULT_CONFIG, metavar="JSON",
                        help=f"Disease config JSON (default: {DEFAULT_CONFIG}).")
    args = parser.parse_args()

    # Validate required files
    missing = [p for p in (args.hg38, args.hap1, args.hap2, args.config) if not os.path.exists(p)]
    if missing:
        for p in missing:
            print(f"ERROR: required file not found: {p}")
        sys.exit(1)

    with open(args.config) as f:
        diseases = json.load(f)

    os.makedirs(args.output_dir, exist_ok=True)

    results = {}
    for d in diseases:
        name = d["name"]
        if args.disease and args.disease.lower() not in name.lower():
            continue
        ok = run_disease(d, args)
        results[name] = "OK" if ok else "FAILED"

    print(f"\n{'='*60}")
    print("Summary:")
    for name, status in results.items():
        print(f"  {status:6s}  {name}")


if __name__ == "__main__":
    main()
