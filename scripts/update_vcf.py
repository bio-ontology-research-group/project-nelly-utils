import json
import os

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
config_path = os.path.join(base_dir, "config", "diseases.json")
vcf_path = os.path.join(base_dir, "output", "project_variants.vcf")

def main():
    if not os.path.exists(config_path):
        print("Config not found")
        return

    with open(config_path, 'r') as f:
        diseases = json.load(f)

    with open(vcf_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"ClinVar Disease Name\">\n")
        f.write("##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Clinical Significance\">\n")
        f.write("##INFO=<ID=OMIM,Number=1,Type=String,Description=\"OMIM ID\">\n")
        f.write("##INFO=<ID=INHERITANCE,Number=1,Type=String,Description=\"Mode of Inheritance\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for i, d in enumerate(diseases):
            var = d['variant']
            # Remove 'chr' prefix for standard VCF if desired, but hg38 uses chr usually.
            # ClinVar often uses numbers. Let's stick to what's in config (e.g. "chr15")
            
            chrom = var['chrom']
            pos = var['pos']
            ref = var['ref']
            alt = var['alt']
            omim_id = d['omim_id']
            inheritance = d.get('inheritance', 'unknown')
            name = d['name']
            
            # Simple ID
            var_id = f"v{i+1}"
            
            info = f"CLNDN=\"{name}\";CLNSIG=Pathogenic;OMIM={omim_id};INHERITANCE={inheritance}"
            
            f.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t{info}\n")

    print(f"Updated VCF at {vcf_path}")

if __name__ == "__main__":
    main()
