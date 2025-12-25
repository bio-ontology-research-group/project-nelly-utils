import json
import gzip
import os
import sys

# Paths
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
config_path = os.path.join(base_dir, "config", "diseases.json")
clinvar_path = os.path.join(base_dir, "data", "variant_summary.txt.gz")

def main():
    if not os.path.exists(config_path):
        print("Config not found.")
        return

    with open(config_path, 'r') as f:
        diseases = json.load(f)

    # Map OMIM ID to disease index in list
    # Normalize OMIM IDs in config to just numbers for easier matching
    omim_map = {}
    for i, d in enumerate(diseases):
        clean_id = d['omim_id'].replace("OMIM:", "")
        omim_map[clean_id] = i

    print(f"Looking for variants for {len(omim_map)} diseases...")

    found_variants = {} # omim_id -> {chrom, pos, ref, alt, type}

    # Open ClinVar
    with gzip.open(clinvar_path, 'rt', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        
        # Map column names to indices
        col_idx = {name: i for i, name in enumerate(header)}
        
        required_cols = ['PhenotypeIDS', 'Assembly', 'ClinicalSignificance', 'Chromosome', 'PositionVCF', 'ReferenceAlleleVCF', 'AlternateAlleleVCF', 'Type']
        for c in required_cols:
            if c not in col_idx:
                print(f"Error: Missing column {c} in ClinVar file.")
                return

        count = 0
        for line in f:
            count += 1
            if count % 100000 == 0:
                print(f"Processed {count} lines...", end='\r')

            parts = line.strip().split('\t')
            
            # Check Assembly
            if parts[col_idx['Assembly']] != 'GRCh38':
                continue

            # Check Significance
            sig = parts[col_idx['ClinicalSignificance']].lower()
            if "pathogenic" not in sig or "conflicting" in sig:
                continue

            # Check Phenotype
            pheno_ids = parts[col_idx['PhenotypeIDS']]
            
            # Optimization: check if any of our target IDs are in string
            # Precise check: split tokens
            matched_omim = None
            
            # Iterate our targets (small list)
            for target_id in omim_map:
                # Look for "OMIM:123456"
                if f"OMIM:{target_id}" in pheno_ids:
                    matched_omim = target_id
                    break
            
            if not matched_omim:
                continue

            # We found a candidate
            chrom = parts[col_idx['Chromosome']]
            pos = parts[col_idx['PositionVCF']]
            ref = parts[col_idx['ReferenceAlleleVCF']]
            alt = parts[col_idx['AlternateAlleleVCF']]
            vtype = parts[col_idx['Type']]

            # Skip if any field is empty or weird
            if not chrom or not pos or not ref or not alt:
                continue
                
            # Prefer SNVs over Indels if we already have one, but take Indel if nothing else
            # Or just take the first Pathogenic one we find?
            # Let's prioritize SNVs.
            
            is_snv = (vtype == 'single nucleotide variant')
            
            # If we don't have this disease yet, take it.
            # If we have it, but it's not SNV and this one IS SNV, take this one.
            
            if matched_omim not in found_variants:
                found_variants[matched_omim] = {
                    'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'type': vtype
                }
            else:
                current = found_variants[matched_omim]
                if current['type'] != 'single nucleotide variant' and is_snv:
                     found_variants[matched_omim] = {
                        'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'type': vtype
                    }

    print("\nUpdating config...")
    
    updated_count = 0
    for omim_id, var_data in found_variants.items():
        idx = omim_map[omim_id]
        d = diseases[idx]
        
        print(f"  Updated {d['name']} (OMIM:{omim_id}): {var_data['chrom']}:{var_data['pos']} {var_data['ref']}>{var_data['alt']} ({var_data['type']})")
        
        # Ensure chrom has "chr" prefix if missing (ClinVar often just has "4", "X")
        c = var_data['chrom']
        if not c.startswith("chr"):
            c = f"chr{c}"
            
        d['variant']['chrom'] = c
        d['variant']['pos'] = var_data['pos']
        d['variant']['ref'] = var_data['ref']
        d['variant']['alt'] = var_data['alt']
        d['note'] = "Updated from ClinVar reference"
        
        updated_count += 1

    print(f"Updated {updated_count} diseases.")
    
    with open(config_path, 'w') as f:
        json.dump(diseases, f, indent=2)

if __name__ == "__main__":
    main()
