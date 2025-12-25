import os
import sys
import json
import subprocess
import pysam
import shutil

# Ensure src is in path
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(base_dir, "src"))

from read_simulator import simulate_reads
from phenotype import create_phenopacket

# Config
HG38_FASTA = os.path.join(base_dir, "hg38.fa")
HAP1_FASTA = os.path.join(base_dir, "data", "GCA_050492415.1_apr041.1_v1_genomic.fna")
HAP2_FASTA = os.path.join(base_dir, "data", "GCA_050492395.1_apr041.2_v1_genomic.fna")
OUTPUT_DIR = os.path.join(base_dir, "output")
CONFIG_PATH = os.path.join(base_dir, "config", "diseases.json")
CONTEXT_WINDOW = 10000  # +/- 10kb = 20kb total context

def get_context_sequence(chrom, pos, ref, window=CONTEXT_WINDOW):
    """
    Extracts sequence from HG38 around the variant.
    Returns: (sequence_string, start_offset_in_window)
    """
    with pysam.FastaFile(HG38_FASTA) as fasta:
        start = max(0, pos - 1 - window)
        end = pos - 1 + len(ref) + window
        
        # Verify chrom exists
        if chrom not in fasta.references:
            if chrom.startswith("chr"): alt_chrom = chrom[3:]
            else: alt_chrom = "chr" + chrom
            if alt_chrom in fasta.references: chrom = alt_chrom
            else: raise ValueError(f"Chromosome {chrom} not found in HG38.")

        seq = fasta.fetch(chrom, start, end)
        variant_offset = pos - 1 - start
        return seq, variant_offset

def map_sequence(sequence, target_fasta):
    """
    Maps the sequence to the target fasta using minimap2.
    Returns: (target_contig, target_pos_1based, is_reverse)
    """
    query_path = os.path.join(OUTPUT_DIR, "temp_query.fa")
    with open(query_path, "w") as f:
        f.write(">query\n")
        f.write(sequence + "\n")
    
    cmd = ["minimap2", "-x", "asm5", "--secondary=no", target_fasta, query_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    best_hit = None
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split("\t")
        if parts[0] == "query":
            best_hit = parts
            break 
            
    if not best_hit: return None
        
    return {
        "contig": best_hit[5],
        "t_start": int(best_hit[7]),
        "t_end": int(best_hit[8]),
        "q_start": int(best_hit[2]),
        "q_end": int(best_hit[3]),
        "strand": best_hit[4]
    }

def calculate_target_pos(mapping, variant_offset):
    """
    Calculates the 1-based target position for the variant.
    """
    if not mapping: return None
    if not (mapping['q_start'] <= variant_offset <= mapping['q_end']):
        return None

    offset_in_mapping = variant_offset - mapping['q_start']
    
    if mapping['strand'] == '+':
        target_pos_0 = mapping['t_start'] + offset_in_mapping
        return target_pos_0 + 1 
    else:
        # Reverse strand:
        # q_start (low index) matches t_end (high index)
        # offset increases from q_start
        # so target pos decreases from t_end
        target_pos_0 = mapping['t_end'] - offset_in_mapping - 1
        return target_pos_0 + 1

def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def extract_and_patch(fasta_path, contig, target_pos, ref, alt, do_patch, strand, window=CONTEXT_WINDOW):
    """
    Extracts a chunk around target_pos and applies patch if requested.
    Returns: (header, sequence)
    """
    with pysam.FastaFile(fasta_path) as fasta:
        contig_len = fasta.get_reference_length(contig)
        
        # Define chunk bounds (0-based)
        # target_pos is 1-based. target_pos-1 is 0-based index of variant start.
        start = max(0, target_pos - 1 - window)
        end = min(contig_len, target_pos - 1 + len(ref) + window) # Approx end
        
        seq_chunk = fasta.fetch(contig, start, end)
        
        # Calculate variant offset within chunk
        # Global pos: target_pos - 1
        # Chunk start: start
        rel_pos = (target_pos - 1) - start
        
        header = f">{contig}:{start+1}-{end}"
        if do_patch:
            header += "_mutated"
            
            # Prepare alleles
            p_ref = ref
            p_alt = alt
            
            # If the mapping was reverse strand, the genomic sequence at this pos
            # is the RC of the hg38 Ref.
            if strand == '-':
                p_ref = rc(ref)
                p_alt = rc(alt)
                
            # Verify Ref
            # chunk sequence is always Forward relative to assembly
            observed_ref = seq_chunk[rel_pos:rel_pos+len(p_ref)]
            
            if observed_ref.upper() != p_ref.upper():
                print(f"    WARNING: Ref mismatch! Expected {p_ref}, Found {observed_ref}. Patching anyway.")
            
            # Apply patch
            seq_chunk = seq_chunk[:rel_pos] + p_alt + seq_chunk[rel_pos+len(p_ref):]
            
        return header, seq_chunk

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    with open(CONFIG_PATH, 'r') as f:
        diseases = json.load(f)

    for d in diseases:
        sample_id = d['name'].lower().replace(" ", "_").replace("'", "")
        print(f"\nProcessing {d['name']} ({sample_id})...")
        
        sample_dir = os.path.join(OUTPUT_DIR, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        
        variant = d['variant']
        chrom = variant['chrom']
        pos = variant['pos']
        ref = variant['ref']
        alt = variant['alt']
        
        # 1. Get Context
        try:
            query_seq, variant_offset = get_context_sequence(chrom, pos, ref)
        except Exception as e:
            print(f"  Error fetching context: {e}")
            continue

        # 2. Map
        map1 = map_sequence(query_seq, HAP1_FASTA)
        target_pos1 = calculate_target_pos(map1, variant_offset)
        
        map2 = map_sequence(query_seq, HAP2_FASTA)
        target_pos2 = calculate_target_pos(map2, variant_offset)

        if not target_pos1 and not target_pos2:
            print("  Failed to map to either haplotype. Skipping.")
            continue

        # 3. Patch Decision
        patch_hap1 = False
        patch_hap2 = False
        if d['inheritance'] == 'dominant':
            if target_pos1: patch_hap1 = True
        else: # Recessive
            if target_pos1: patch_hap1 = True
            if target_pos2: patch_hap2 = True

        # 4. Extract & Patch ROI
        chunk_fa_path = os.path.join(sample_dir, "roi_diploid.fa")
        
        with open(chunk_fa_path, "w") as out_f:
            # Hap1
            if target_pos1:
                h, s = extract_and_patch(HAP1_FASTA, map1['contig'], target_pos1, ref, alt, patch_hap1, map1['strand'])
                out_f.write(f"{h}_hap1\n{s}\n")
            
            # Hap2
            if target_pos2:
                h, s = extract_and_patch(HAP2_FASTA, map2['contig'], target_pos2, ref, alt, patch_hap2, map2['strand'])
                out_f.write(f"{h}_hap2\n{s}\n")
                
        # 5. Simulate Reads
        # Coverage 15x on the diploid file = 15x per haplotype = 30x total
        print("  Simulating 30x ROI reads...")
        simulate_reads(chunk_fa_path, os.path.join(sample_dir, sample_id), 15, 150, simulator_path="art_illumina")

        # 6. Phenopacket
        create_phenopacket(d['omim_id'], sample_id, os.path.join(sample_dir, f"{sample_id}.phenopacket.json"))

if __name__ == "__main__":
    main()

