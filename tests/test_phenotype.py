import pytest
import json
from unittest.mock import patch
from src.phenotype import create_phenopacket
from phenopackets import OntologyClass


def test_create_phenopacket(tmp_path):
    """
    Tests create_phenopacket with the current API (omim_id, sample_id).
    Mocks HPOA lookup to avoid dependency on the data file.
    """
    output_path = str(tmp_path / "phenopacket.json")
    sample_id = "test_sample"
    omim_id = "143100"

    fake_hpo_terms = [OntologyClass(id="HP:0000001", label="All")]

    with patch("src.phenotype.get_hpo_terms_from_local_db", return_value=fake_hpo_terms):
        create_phenopacket(
            omim_id=omim_id,
            sample_id=sample_id,
            output_phenopacket_path=output_path,
        )

    with open(output_path) as f:
        data = json.load(f)

    assert data["id"] == f"{sample_id}-phenopacket"
    assert data["subject"]["id"] == sample_id
    assert len(data["phenotypicFeatures"]) == 1
    assert data["phenotypicFeatures"][0]["type"]["id"] == "HP:0000001"
