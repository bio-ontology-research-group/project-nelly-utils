from phenopackets import (
    Phenopacket, Individual, PhenotypicFeature, OntologyClass,
    Interpretation, Diagnosis, GenomicInterpretation, VariantInterpretation,
    VariationDescriptor, VcfRecord, MoleculeContext, AcmgPathogenicityClassification
)
from google.protobuf.json_format import MessageToJson
import os
import csv

def get_hpo_terms_from_local_db(omim_id: str) -> list[OntologyClass]:
    """
    Retrieves HPO terms for a given OMIM ID from the local phenotype.hpoa file.
    """
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    hpoa_path = os.path.join(base_dir, "data", "phenotype.hpoa")
    
    if not os.path.exists(hpoa_path):
        print(f"Error: HPOA file not found at {hpoa_path}")
        return []

    if not omim_id.startswith("OMIM:"):
        query_id = f"OMIM:{omim_id}"
    else:
        query_id = omim_id

    hpo_ids = set()
    
    try:
        with open(hpoa_path, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split('\t')
                if len(parts) < 4:
                    continue
                if parts[0] == query_id:
                    hpo_ids.add(parts[3])
    except Exception as e:
        print(f"Error reading HPOA file: {e}")
        return []
    
    return [OntologyClass(id=hid, label="") for hid in hpo_ids]


def create_phenopacket(omim_id: str, sample_id: str, output_phenopacket_path: str, variant_info: dict = None, inheritance: str = None):
    """
    Creates a Phenopacket for a sample.

    Args:
        omim_id: The OMIM ID (e.g., "143100").
        sample_id: The ID of the sample.
        output_phenopacket_path: Path to write JSON.
        variant_info: Dict with chrom, pos, ref, alt.
        inheritance: "dominant" or "recessive".
    """
    
    # 1. Get HPO terms
    hpo_terms = get_hpo_terms_from_local_db(omim_id)
    if not hpo_terms:
        print(f"Warning: No HPO terms found for {omim_id} in local DB.")

    # 2. Subject
    subject = Individual(id=sample_id)
    
    # 3. Phenotypic Features
    phenotypic_features = [PhenotypicFeature(type=term) for term in hpo_terms]
    
    # 4. Interpretations
    interpretations = []
    if variant_info:
        # Determine Zygosity
        allelic_state = None
        if inheritance == "dominant":
            allelic_state = OntologyClass(id="GENO:0000135", label="heterozygous")
        elif inheritance == "recessive":
            allelic_state = OntologyClass(id="GENO:0000136", label="homozygous")

        # Create VcfRecord
        vcf_record = VcfRecord(
            genome_assembly=variant_info.get('genome_assembly', 'GRCh38'),
            chrom=variant_info.get('chrom', 'unknown'),
            pos=int(variant_info.get('pos', 0)),
            ref=variant_info.get('ref', 'N'),
            alt=variant_info.get('alt', 'N')
        )
        
        # Create VariationDescriptor
        variation_descriptor = VariationDescriptor(
            id=f"var-{variant_info.get('chrom')}-{variant_info.get('pos')}",
            vcf_record=vcf_record,
            molecule_context=MoleculeContext.genomic,
            allelic_state=allelic_state
        )
        
        # Create VariantInterpretation
        variant_interpretation = VariantInterpretation(
            variation_descriptor=variation_descriptor,
            acmg_pathogenicity_classification=AcmgPathogenicityClassification.PATHOGENIC
        )
        
        # Create GenomicInterpretation
        genomic_interpretation = GenomicInterpretation(
            subject_or_biosample_id=sample_id,
            interpretation_status=GenomicInterpretation.InterpretationStatus.CAUSATIVE,
            variant_interpretation=variant_interpretation
        )
        
        # Create Diagnosis
        diagnosis = Diagnosis(
            genomic_interpretations=[genomic_interpretation],
            disease=OntologyClass(id=omim_id, label="See OMIM")
        )
        
        # Create Interpretation
        interpretation = Interpretation(
            id=f"{sample_id}-interpretation",
            progress_status=Interpretation.ProgressStatus.SOLVED,
            diagnosis=diagnosis
        )
        
        interpretations.append(interpretation)

    # 5. Assemble Phenopacket
    phenopacket = Phenopacket(
        id=f"{sample_id}-phenopacket",
        subject=subject,
        phenotypic_features=phenotypic_features,
        interpretations=interpretations
    )

    # 6. Write to JSON
    with open(output_phenopacket_path, 'w') as f:
        f.write(MessageToJson(phenopacket))
    
    print(f"Phenopacket created for {sample_id} at {output_phenopacket_path}")

if __name__ == '__main__':
    pass
