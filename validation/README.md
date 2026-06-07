# Simulator validation harness

End-to-end check that the haplotype-aware simulators plant disease variants at
the **exact** hg38 coordinate with the **correct zygosity**, by round-tripping
through read simulation + alignment + variant calling.

## Why this exists

`scripts/map_and_simulate.py`, `scripts/simulate_master_genome.py` and
`scripts/haplotype_full_sim.py` lift hg38 variant coordinates onto the apr041
diploid assemblies. The original liftover interpolated linearly across the
minimap2 alignment, ignoring indels, so variants drifted off-target (and by a
*different* amount per haplotype, splitting homozygous variants into two
heterozygous calls). The fix walks the minimap2 CIGAR (`-c` / `cg:Z`). These
scripts prove the bug and the fix.

## Requirements

`hg38.fa` + both `data/GCA_*apr041*.fna` in the repo root/`data`, plus
`art_illumina`, `minimap2`, `samtools`, `bcftools` on PATH. Use the **system
`python3`** (it has `pysam`); the `.arvados-venv` does not.

## Scripts

```bash
# Analytic: buggy linear pos vs indel-aware CIGAR pos, per haplotype
python3 validation/confirm_drift.py [disease-substring]

# Full round-trip for one disease: simulate -> align(chr-only hg38) -> call -> compare
#   mode = buggy (old linear) | fixed (real patched calculate_target_pos)
python3 validation/roundtrip.py tay fixed

# Whole SNP panel pass/fail table
python3 validation/validate_all.py fixed
```

(Large intermediates are written under `validation/work/` and can be deleted.)

## Result (2026-06-06)

| mode  | recovered at exact pos | correct zygosity |
|-------|------------------------|------------------|
| buggy | 1 / 9 SNP diseases     | —                |
| fixed | **9 / 9**              | **9 / 9**        |

Regression guard: `tests/test_integration.py::test_minimap2_liftover_indel_aware`.
