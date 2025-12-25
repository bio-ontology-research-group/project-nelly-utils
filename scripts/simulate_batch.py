import os
import subprocess
import sys
import json

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(base_dir, "output")
    config_path = os.path.join(base_dir, "config", "diseases.json")
    
    # We use the hg38.fa found in the root directory
    fasta_path = os.path.join(base_dir, "hg38.fa")
    
    if not os.path.exists(config_path):
        print(f"Error: Config file not found at {config_path}")
        return

    if not os.path.exists(fasta_path):
        print(f"Error: Reference genome not found at {fasta_path}")
        return

    with open(config_path, 'r') as f:
        diseases = json.load(f)

    os.makedirs(output_dir, exist_ok=True)
    
    for d in diseases:
        sample_id = d['name'].lower().replace(" ", "_").replace("'", "")
        print(f"\nProcessing {sample_id} ({d['inheritance']})...")
        
        # Determine haplotype arg
        # Dominant: 'hap1' (heterozygous)
        # Recessive: 'both' (homozygous)
        hap_arg = 'hap1' if d['inheritance'] == 'dominant' else 'both'
        
        variant = d['variant']
        
        cmd = [
            sys.executable,
            os.path.join(base_dir, "src", "main.py"),
            "--fasta", fasta_path,
            "--chrom", variant['chrom'],
            "--pos", str(variant['pos']),
            "--ref", variant['ref'],
            "--alt", variant['alt'],
            "--omim-id", d['omim_id'],
            "--output-dir", os.path.join(output_dir, sample_id),
            "--sample-id", sample_id,
            "--haplotype", hap_arg,
            "--coverage", "30",
            "--art-path", "art_illumina" # Will trigger fallback if not found
        ]
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {sample_id}: {e}")

if __name__ == "__main__":
    main()
