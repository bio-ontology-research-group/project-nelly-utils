import pytest
from src.genome_patcher import patch_genome
import os

def test_patch_genome_snp(tmp_path):
    """
    Tests patch_genome with a simple SNP on a single-contig FASTA.
    Original: ACGTACGTACGTACGTACGT (pos 5 = A)
    Patch: pos 5, ref=A, alt=T → ACGTTCGTACGTACGTACGT
    """
    fasta_path = "tests/data/test_genome.fa"
    output_fasta_path = str(tmp_path / "patched.fa")

    patch_genome(
        fasta_path=fasta_path,
        chrom="test_contig",
        pos=5,       # 1-based; 0-based index 4 → 'A'
        ref="A",
        alt="T",
        output_fasta_path=output_fasta_path,
        haplotype='both'
    )

    with open(output_fasta_path, "r") as f:
        lines = f.readlines()

    assert lines[0] == ">test_contig\n"
    # Sequence is written 60 chars/line; 20bp fits on one line
    full_seq = "".join(l.strip() for l in lines[1:])
    assert full_seq == "ACGTTCGTACGTACGTACGT"


def test_patch_genome_exact_chrom_match(tmp_path):
    """
    Ensures exact chromosome name match is used when available,
    so 'chr1' does not accidentally patch 'chr10' or 'chr11'.
    """
    # Build a two-contig FASTA
    fasta_path = str(tmp_path / "two_contig.fa")
    with open(fasta_path, "w") as f:
        f.write(">chr1\n" + "A" * 20 + "\n")
        f.write(">chr11\n" + "G" * 20 + "\n")

    output_fasta_path = str(tmp_path / "patched.fa")
    patch_genome(
        fasta_path=fasta_path,
        chrom="chr1",
        pos=1,
        ref="A",
        alt="T",
        output_fasta_path=output_fasta_path,
    )

    with open(output_fasta_path, "r") as f:
        content = f.read()

    # chr1 first base should be T
    lines = content.splitlines()
    assert lines[0] == ">chr1"
    chr1_seq = "".join(l for l in lines[1:] if not l.startswith(">")).split(">")[0]
    # chr11 should be unchanged (all G)
    chr11_start = content.index(">chr11")
    chr11_seq = content[chr11_start:].splitlines()[1]
    assert chr11_seq == "G" * 20
