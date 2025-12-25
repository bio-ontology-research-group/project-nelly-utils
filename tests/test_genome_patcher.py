import pytest
from src.genome_patcher import patch_genome
import os

def test_patch_genome_snp(tmp_path):
    """
    Tests the patch_genome function with a simple SNP.
    """
    test_data_dir = "tests/data"
    fasta_path = os.path.join(test_data_dir, "test_genome.fa")
    vcf_path = os.path.join(test_data_dir, "test_clinvar.vcf")
    output_fasta_path = tmp_path / "patched.fa"
    variant_id = "TEST001"

    # Run the patcher
    patch_genome(
        fasta_path=fasta_path,
        vcf_path=vcf_path,
        variant_id=variant_id,
        output_fasta_path=str(output_fasta_path),
        haplotype='both'
    )

    # Check the output
    with open(output_fasta_path, "r") as f:
        lines = f.readlines()
        
    # Should have a header and a sequence line
    assert len(lines) == 2
    # The header should be the same
    assert lines[0] == ">test_contig\n"
    # The sequence should be patched
    expected_sequence = "ACGTGCGTACGTACGTACGT\n"
    assert lines[1] == expected_sequence
