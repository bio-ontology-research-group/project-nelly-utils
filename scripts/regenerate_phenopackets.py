import os
import sys
import json

# Ensure src is in path
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(base_dir, "src"))

from phenotype import create_phenopacket

# Config
OUTPUT_DIR = os.path.join(base_dir, "output")
CONFIG_PATH = os.path.join(base_dir, "config", "diseases.json")

def main():
    if not os.path.exists(CONFIG_PATH):
        print(f"Error: Config file not found at {CONFIG_PATH}")
        return

    with open(CONFIG_PATH, 'r') as f:
        diseases = json.load(f)

    for d in diseases:
        sample_id = d['name'].lower().replace(" ", "_").replace("'", "")
        print(f"Regenerating Phenopacket for {d['name']} ({sample_id})...")
        
        sample_dir = os.path.join(OUTPUT_DIR, sample_id)
        # Ensure dir exists (it should from simulation)
        os.makedirs(sample_dir, exist_ok=True)
        
        pheno_path = os.path.join(sample_dir, f"{sample_id}.phenopacket.json")
        
        # Prepare variant info
        variant = d['variant']
        variant_info = {
            'chrom': variant['chrom'],
            'pos': variant['pos'],
            'ref': variant['ref'],
            'alt': variant['alt'],
            'genome_assembly': 'GRCh38'
        }
        
        # We need to pass "OMIM:ID" format. The JSON has just ID (e.g. "143100").
        omim_id_arg = d['omim_id']
        if not omim_id_arg.startswith("OMIM:"):
            omim_id_arg = f"OMIM:{omim_id_arg}"
            
        create_phenopacket(
            omim_id=omim_id_arg, 
            sample_id=sample_id, 
            output_phenopacket_path=pheno_path, 
            variant_info=variant_info,
            inheritance=d.get('inheritance')
        )

if __name__ == "__main__":
    main()
