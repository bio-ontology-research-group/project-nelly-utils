import pytest
import pysam
import json
from src.phenotype import create_phenopacket
import os

def test_create_phenopacket(tmp_path):
    """
    Tests the create_phenopacket function.
    """
    test_data_dir = "tests/data"
    vcf_path = os.path.join(test_data_dir, "test_clinvar.vcf")
    output_phenopacket_path = tmp_path / "phenopacket.json"
    sample_id = "test_sample"

    # Get the variant record from the VCF file
    vcf_file = pysam.VariantFile(vcf_path)
    variant_record = None
    for record in vcf_file.fetch():
        if record.id == "TEST001":
            variant_record = record
            break
    vcf_file.close()

    assert variant_record is not None

    # Run the function
    create_phenopacket(
        vcf_record=variant_record,
        sample_id=sample_id,
        output_phenopacket_path=str(output_phenopacket_path)
    )

    # Check the output
    with open(output_phenopacket_path, "r") as f:
        phenopacket_data = json.load(f)

    assert phenopacket_data['id'] == f"{sample_id}-phenopacket"
    assert phenopacket_data['subject']['id'] == sample_id
    assert len(phenopacket_data['phenotypicFeatures']) == 1
    assert phenopacket_data['phenotypicFeatures'][0]['type']['id'] == 'HP:0000001'
    assert phenopacket_data['phenotypicFeatures'][0]['type']['label'] == 'All'
