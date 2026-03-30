"""
tests/test_integration.py — End-to-end simulation integration test.

Genome design
-------------
30 kb, composed of 30 × 1000 bp homopolymer blocks, cycling:
    block  0 (    0–  999) : A × 1000
    block  1 ( 1000– 1999) : C × 1000
    block  2 ( 2000– 2999) : G × 1000
    block  3 ( 3000– 3999) : T × 1000
    block  4 ( 4000– 4999) : A × 1000  … (pattern repeats)
    …
    block 12 (12000–12999) : A × 1000
    block 13 (13000–13999) : C × 1000
    block 14 (14000–14999) : G × 1000   ← VARIANT HERE
    block 15 (15000–15999) : T × 1000
    …

Variant: G → A at 1-based position 14501 (0-based: 14500), inside the all-G block.

This creates the unique 31-mer:
    G × 15  +  A  +  G × 15   =  "GGGGGGGGGGGGGGGAGGGGGGGGGGGGGGGG"

Uniqueness guarantee
--------------------
In the original 30 kb genome, G-blocks are only flanked by C-blocks and T-blocks
(the cycle is A→C→G→T). Therefore:
  • "A" never appears inside a G-block  →  the k-mer G+A+G is absent in orignal.
  • RC of the k-mer is C×15 + T + C×15. C-blocks are flanked by A and G blocks,
    never by T-blocks, so the RC k-mer is also absent in the original.

Tests
-----
  test_unique_kmer_not_in_original_genome
      Sanity-check: k-mer and its RC are absent from the unmodified sequence.

  test_patch_creates_unique_kmer
      genome_patcher introduces the k-mer at the expected position.

  test_variant_recoverable_from_reads          ← the main integration test
      1. Simulate error-free reads from UNPATCHED diploid → k-mer must be ABSENT.
      2. Simulate error-free reads from PATCHED diploid   → k-mer must be PRESENT.
      3. Number of variant-carrying reads is consistent with 5× per-haplotype coverage.

  test_minimap2_liftover    (skipped if minimap2 not in PATH)
      Lift over a known position from a "fake hg38" (same sequence) to a fake
      assembly using the liftover helpers from haplotype_full_sim.
      Uses the entire 30 kb as the minimap2 query window so that ambiguous
      sub-block k-mers are not an issue.
"""

import os
import sys
import shutil
import random
import pytest

# Make src/ importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

from genome_patcher import patch_genome
from read_simulator import simulate_reads_python

# ---------------------------------------------------------------------------
# Genome constants
# ---------------------------------------------------------------------------
BLOCK_SIZE   = 1_000
NUM_BLOCKS   = 30
GENOME_LEN   = BLOCK_SIZE * NUM_BLOCKS    # 30,000 bp per haplotype
MOTIFS       = ['A', 'C', 'G', 'T']

# Variant lives in block 14 (14 % 4 == 2 → all-G block)
VARIANT_BLOCK = 14
VARIANT_POS_0 = BLOCK_SIZE * VARIANT_BLOCK + BLOCK_SIZE // 2   # 14500 (0-based)
VARIANT_POS_1 = VARIANT_POS_0 + 1                               # 14501 (1-based)
VARIANT_REF   = 'G'
VARIANT_ALT   = 'A'

# The unique 31-mer created by the variant (G×15 + A + G×15)
KMER_R        = 15
UNIQUE_KMER   = 'G' * KMER_R + 'A' + 'G' * KMER_R

_RC = str.maketrans('ACGTacgt', 'TGCAtgca')

def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]

RC_UNIQUE_KMER = revcomp(UNIQUE_KMER)   # C×15 + T + C×15


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_genome_seq() -> str:
    """Return the 30 kb reference sequence (no variant)."""
    return ''.join(MOTIFS[i % 4] * BLOCK_SIZE for i in range(NUM_BLOCKS))


def write_fasta(path: str, contig: str, seq: str) -> None:
    """Write sequence as a well-formed FASTA (60-char line wrap)."""
    with open(path, 'w') as f:
        f.write(f'>{contig}\n')
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + '\n')


def read_fastq_seqs(fq_path: str) -> list:
    seqs = []
    with open(fq_path) as f:
        lines = f.readlines()
    for i in range(0, len(lines), 4):
        seqs.append(lines[i + 1].strip())
    return seqs


def kmer_in_reads(kmer: str, fq1: str, fq2: str) -> bool:
    """Return True if kmer or its RC appears in any read (R1 or R2)."""
    rc = revcomp(kmer)
    for fq in (fq1, fq2):
        for seq in read_fastq_seqs(fq):
            if kmer in seq or rc in seq:
                return True
    return False


def count_reads_with_kmer(kmer: str, fq1: str, fq2: str) -> int:
    """Count reads (across R1 + R2) that contain kmer or its RC."""
    rc = revcomp(kmer)
    count = 0
    for fq in (fq1, fq2):
        for seq in read_fastq_seqs(fq):
            if kmer in seq or rc in seq:
                count += 1
    return count


# ---------------------------------------------------------------------------
# Tests — k-mer and patching
# ---------------------------------------------------------------------------

def test_unique_kmer_not_in_original_genome():
    """
    Sanity-check: the unique 31-mer and its RC must be absent from the
    unmodified genome.  If this fails, the genome design is wrong.
    """
    seq = build_genome_seq()
    assert UNIQUE_KMER not in seq, (
        f"UNIQUE_KMER already present in the original genome at "
        f"index {seq.index(UNIQUE_KMER)}. Redesign the test genome."
    )
    assert RC_UNIQUE_KMER not in seq, (
        "RC of UNIQUE_KMER already present in the original genome. "
        "Redesign the test genome."
    )


def test_patch_creates_unique_kmer(tmp_path):
    """
    After genome_patcher applies the SNP, the unique 31-mer must appear at
    exactly the expected position and nowhere else.
    """
    ref_fa     = str(tmp_path / 'ref.fa')
    patched_fa = str(tmp_path / 'patched.fa')
    write_fasta(ref_fa, 'hap', build_genome_seq())

    patch_genome(ref_fa, 'hap', VARIANT_POS_1, VARIANT_REF, VARIANT_ALT, patched_fa)

    with open(patched_fa) as f:
        seq = ''.join(l.strip() for l in f if not l.startswith('>'))

    assert seq[VARIANT_POS_0] == VARIANT_ALT, (
        f"Expected '{VARIANT_ALT}' at 0-based position {VARIANT_POS_0}, "
        f"got '{seq[VARIANT_POS_0]}'."
    )
    assert UNIQUE_KMER in seq, "Unique k-mer not found in patched sequence."
    assert seq.count(UNIQUE_KMER) == 1, "Unique k-mer appears more than once — unexpected."


# ---------------------------------------------------------------------------
# Main integration test
# ---------------------------------------------------------------------------

def test_variant_recoverable_from_reads(tmp_path):
    """
    Full pipeline test: patch → build diploid FASTA → simulate → find variant.

    Diploid layout
    --------------
    Both haplotypes carry the same variant (homozygous / recessive).
    This maximises the expected number of variant-carrying reads.

    Coverage used: 5× per haplotype → ~10× diploid.
    Diploid genome size: 2 × 30 kb = 60 kb.
    Expected read pairs: (10 × 60,000) / 150 = 4,000.

    For each haplotype, reads that span position 14500:
      fragment covers variant if start ∈ [14500 - fragment_size, 14500]
      ≈ 500 / 30,000 ≈ 1.7 % of reads per haplotype.
      At 2,000 reads per haplotype → ~34 variant-spanning reads expected.
    We assert ≥ 3, a very conservative bound.

    Negative control
    ----------------
    Reads from the UNPATCHED diploid must contain zero copies of the k-mer.
    """
    random.seed(42)   # make the Python simulator deterministic
    ref_seq = build_genome_seq()

    # ── Negative control: reads from unpatched diploid ─────────────────────
    unpatched_fa = str(tmp_path / 'unpatched_diploid.fa')
    with open(unpatched_fa, 'w') as f:
        for hap in ('hap1', 'hap2'):
            f.write(f'>{hap}\n')
            for i in range(0, len(ref_seq), 60):
                f.write(ref_seq[i:i + 60] + '\n')

    ctrl_prefix = str(tmp_path / 'ctrl_')
    simulate_reads_python(unpatched_fa, ctrl_prefix, coverage=5, read_length=150)

    ctrl_r1, ctrl_r2 = ctrl_prefix + '1.fq', ctrl_prefix + '2.fq'
    assert not kmer_in_reads(UNIQUE_KMER, ctrl_r1, ctrl_r2), (
        "Unique k-mer found in reads from UNPATCHED genome — "
        "the k-mer is not unique or the simulator is broken."
    )

    # ── Positive test: patch both haplotypes, build diploid ────────────────
    random.seed(42)

    hap1_fa  = str(tmp_path / 'hap1.fa')
    hap1_pat = str(tmp_path / 'hap1_patched.fa')
    hap2_fa  = str(tmp_path / 'hap2.fa')
    hap2_pat = str(tmp_path / 'hap2_patched.fa')

    write_fasta(hap1_fa, 'hap1', ref_seq)
    write_fasta(hap2_fa, 'hap2', ref_seq)

    patch_genome(hap1_fa, 'hap1', VARIANT_POS_1, VARIANT_REF, VARIANT_ALT, hap1_pat)
    patch_genome(hap2_fa, 'hap2', VARIANT_POS_1, VARIANT_REF, VARIANT_ALT, hap2_pat)

    diploid_fa = str(tmp_path / 'diploid_patched.fa')
    with open(diploid_fa, 'wb') as out:
        for src in (hap1_pat, hap2_pat):
            with open(src, 'rb') as inp:
                shutil.copyfileobj(inp, out)

    pat_prefix = str(tmp_path / 'patient_')
    simulate_reads_python(diploid_fa, pat_prefix, coverage=5, read_length=150)

    pat_r1, pat_r2 = pat_prefix + '1.fq', pat_prefix + '2.fq'

    assert kmer_in_reads(UNIQUE_KMER, pat_r1, pat_r2), (
        "Unique k-mer NOT found in any read from the PATCHED genome. "
        "The inserted variant was not recovered."
    )

    n = count_reads_with_kmer(UNIQUE_KMER, pat_r1, pat_r2)
    print(f"\n    Variant k-mer found in {n} reads "
          f"(5× coverage, 60 kb diploid, seed=42)")
    assert n >= 3, (
        f"Only {n} reads carry the variant k-mer — expected ≥ 3 at 5× coverage. "
        "Coverage calculation or simulator may be broken."
    )


# ---------------------------------------------------------------------------
# Minimap2 liftover test
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not shutil.which('minimap2'), reason='minimap2 not in PATH')
def test_minimap2_liftover(tmp_path):
    """
    Verify the coordinate liftover helpers from haplotype_full_sim.py.

    Setup
    -----
    Both 'fake hg38' and 'fake assembly' carry identical sequence.
    We use the ENTIRE 30 kb as the minimap2 query window, which sidesteps
    the repeated-block ambiguity: a 30 kb unique query maps to exactly one
    position in a 30 kb target.

    Expected result
    ---------------
    The lifted-over position in the assembly equals VARIANT_POS_1 = 14501.
    """
    from haplotype_full_sim import (
        fetch_hg38_context,
        map_context_to_assembly,
        liftover_position,
    )

    ref_seq = build_genome_seq()

    hg38_fa = str(tmp_path / 'hg38_fake.fa')
    hap_fa  = str(tmp_path / 'hap_fake.fa')

    # Both use contig name 'chr_test' in the fake hg38,
    # but the assembly contig is named 'hap_contig' (different name, same sequence)
    write_fasta(hg38_fa, 'chr_test', ref_seq)
    write_fasta(hap_fa,  'hap_contig', ref_seq)

    # Use a large window (15 000 bp either side) so the context = full genome
    context_seq, variant_offset = fetch_hg38_context(
        hg38_fa, 'chr_test', VARIANT_POS_1, VARIANT_REF, window=15_000
    )

    assert variant_offset == VARIANT_POS_0, (
        f"variant_offset={variant_offset}, expected {VARIANT_POS_0}"
    )
    assert context_seq[variant_offset].upper() == VARIANT_REF.upper(), (
        f"Base at variant_offset in context is "
        f"'{context_seq[variant_offset]}', expected '{VARIANT_REF}'"
    )

    hit = map_context_to_assembly(context_seq, hap_fa, str(tmp_path))
    assert hit is not None, "minimap2 returned no hit for the context query."
    assert hit['contig'] == 'hap_contig', (
        f"Expected contig 'hap_contig', got '{hit['contig']}'"
    )

    lifted_pos = liftover_position(hit, variant_offset)
    assert lifted_pos == VARIANT_POS_1, (
        f"Lifted-over position = {lifted_pos}, expected {VARIANT_POS_1}. "
        f"hit={hit}, variant_offset={variant_offset}"
    )
    print(f"\n    minimap2 liftover: chr_test:{VARIANT_POS_1} → "
          f"hap_contig:{lifted_pos}  strand={hit['strand']}")
