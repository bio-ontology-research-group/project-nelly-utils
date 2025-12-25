import os
import sys
import json
import subprocess
import pysam
import shutil
from collections import defaultdict

# Ensure src is in path
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(base_dir, "src"))

from read_simulator import simulate_reads_python 
# Note: importing simulate_reads_python directly to have finer control if needed, 
# or I'll just use the wrapper but I need to be careful about file appending.
# Actually, I'll use a custom simulation function here that accepts a list of intervals to process.

# Config
HAP1_FASTA = os.path.join(base_dir, "data", "GCA_050492415.1_apr041.1_v1_genomic.fna")
HAP2_FASTA = os.path.join(base_dir, "data", "GCA_050492395.1_apr041.2_v1_genomic.fna")
OUTPUT_DIR = os.path.join(base_dir, "output", "master_sim")
CONFIG_PATH = os.path.join(base_dir, "config", "diseases.json")
HG38_FASTA = os.path.join(base_dir, "hg38.fa")
ROI_PADDING = 10000 # +/- 10kb

def get_context_sequence(chrom, pos, ref, window=ROI_PADDING):
    # Same helper as before
    with pysam.FastaFile(HG38_FASTA) as fasta:
        start = max(0, pos - 1 - window)
        end = pos - 1 + len(ref) + window
        if chrom not in fasta.references:
            if "chr"+chrom in fasta.references: chrom = "chr"+chrom
            elif chrom.startswith("chr") and chrom[3:] in fasta.references: chrom = chrom[3:]
        
        if chrom not in fasta.references: return None, 0
        seq = fasta.fetch(chrom, start, end)
        variant_offset = pos - 1 - start
        return seq, variant_offset

def map_sequence(sequence, target_fasta):
    # Same helper as before
    query_path = os.path.join(OUTPUT_DIR, "temp_query.fa")
    with open(query_path, "w") as f:
        f.write(">query\n")
        f.write(sequence + "\n")
    
    cmd = ["minimap2", "-x", "asm5", "--secondary=no", target_fasta, query_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0: return None
    
    for line in result.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split("\t")
        if parts[0] == "query":
            return {
                "contig": parts[5],
                "t_start": int(parts[7]),
                "t_end": int(parts[8]),
                "q_start": int(parts[2]),
                "q_end": int(parts[3]),
                "strand": parts[4]
            }
    return None

def calculate_target_pos(mapping, variant_offset):
    # Same helper
    if not mapping: return None
    if not (mapping['q_start'] <= variant_offset <= mapping['q_end']): return None
    offset = variant_offset - mapping['q_start']
    if mapping['strand'] == '+':
        return mapping['t_start'] + offset + 1
    else:
        return mapping['t_end'] - offset # 1-based logic handled? t_end is exclusive?
        # minimap: start (0-based), end (0-based exclusive)
        # if q match t reversed:
        # q[0] matches t[end-1]
        # q[offset] matches t[end-1-offset]
        # +1 for 1-based = t[end-1-offset] + 1 = end - offset
        return mapping['t_end'] - offset

def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def simulate_region(fasta_obj, contig, start, end, out_prefix, coverage=15, read_len=150):
    """
    Simulates reads for a specific region [start, end) of a contig.
    Using Python fallback logic for precise control over regions.
    """
    # Just extract sequence
    seq = fasta_obj.fetch(contig, start, end)
    length = len(seq)
    if length < read_len: return
    
    # Num pairs
    n_reads = int((coverage * length) / read_len)
    n_pairs = n_reads // 2
    
    import random
    
    fq1 = open(f"{out_prefix}1.fq", "a")
    fq2 = open(f"{out_prefix}2.fq", "a")
    
    for _ in range(n_pairs):
        frag_len = int(random.gauss(500, 50))
        if frag_len < read_len: frag_len = read_len + 10
        if frag_len >= length: frag_len = length - 1
        
        f_start = random.randint(0, length - frag_len)
        fragment = seq[f_start : f_start + frag_len]
        
        r1 = fragment[:read_len]
        r2 = rc(fragment[-read_len:])
        
        qual = 'I' * read_len
        
        # Unique ID
        rid = f"@{contig}:{start+f_start}-{start+f_start+frag_len}:{random.randint(0,1000000)}"
        
        fq1.write(f"{rid}/1\n{r1}\n+\n{qual}\n")
        fq2.write(f"{rid}/2\n{r2}\n+\n{qual}\n")
        
    fq1.close()
    fq2.close()


def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        
    with open(CONFIG_PATH, 'r') as f:
        diseases = json.load(f)

    # 1. Identify ROIs on Hap1 and Hap2
    # Structure: exclusions[hap_name][contig] = list of (start, end, disease_index)
    exclusions = {
        'hap1': defaultdict(list),
        'hap2': defaultdict(list)
    }
    
    disease_metadata = []

    print("Mapping variants...")
    for i, d in enumerate(diseases):
        var = d['variant']
        # Get context
        seq, offset = get_context_sequence(var['chrom'], var['pos'], var['ref'])
        if not seq:
            print(f"Skipping {d['name']} (Context fetch failed)")
            disease_metadata.append(None)
            continue
            
        # Map Hap1
        m1 = map_sequence(seq, HAP1_FASTA)
        pos1 = calculate_target_pos(m1, offset)
        
        # Map Hap2
        m2 = map_sequence(seq, HAP2_FASTA)
        pos2 = calculate_target_pos(m2, offset)
        
        meta = {
            'name': d['name'],
            'inheritance': d.get('inheritance', 'recessive'),
            'ref': var['ref'],
            'alt': var['alt'],
            'hap1_map': None,
            'hap2_map': None
        }
        
        if pos1:
            # ROI Window (0-based)
            # pos1 is 1-based variant start.
            roi_start = max(0, pos1 - 1 - ROI_PADDING)
            roi_end = pos1 - 1 + len(var['ref']) + ROI_PADDING
            exclusions['hap1'][m1['contig']].append((roi_start, roi_end, i))
            meta['hap1_map'] = {'contig': m1['contig'], 'roi_start': roi_start, 'roi_end': roi_end, 'pos': pos1, 'strand': m1['strand']}
            
        if pos2:
            roi_start = max(0, pos2 - 1 - ROI_PADDING)
            roi_end = pos2 - 1 + len(var['ref']) + ROI_PADDING
            exclusions['hap2'][m2['contig']].append((roi_start, roi_end, i))
            meta['hap2_map'] = {'contig': m2['contig'], 'roi_start': roi_start, 'roi_end': roi_end, 'pos': pos2, 'strand': m2['strand']}
            
        disease_metadata.append(meta)

    # Sort exclusions
    for hap in exclusions:
        for contig in exclusions[hap]:
            exclusions[hap][contig].sort()

    # 2. Simulate Background (Skipping Exclusions)
    # We will simulate "Background" and "Loci" separately.
    
    print("Simulating Background (this may take time)...")
    # To save time in this demo, we ONLY simulate background for contigs that HAVE exclusions (plus maybe small others).
    # In a full run, we'd do all contigs.
    
    bg_out_prefix = os.path.join(OUTPUT_DIR, "background")
    # Clear prev
    open(f"{bg_out_prefix}1.fq", "w").close()
    open(f"{bg_out_prefix}2.fq", "w").close()

    for hap_label, fasta_path in [('hap1', HAP1_FASTA), ('hap2', HAP2_FASTA)]:
        f = pysam.FastaFile(fasta_path)
        # Filter contigs to save time? Or do all?
        # Let's do ONLY contigs involved in diseases to prove the point, 
        # otherwise 3GB simulation takes hours.
        involved_contigs = set(exclusions[hap_label].keys())
        
        for contig in involved_contigs: # Iterate only relevant contigs
            length = f.get_reference_length(contig)
            excls = exclusions[hap_label][contig]
            
            curr = 0
            for (estart, eend, d_idx) in excls:
                # Simulate [curr, estart)
                if estart > curr:
                    simulate_region(f, contig, curr, estart, bg_out_prefix)
                curr = max(curr, eend)
            
            # Simulate tail [curr, length)
            if curr < length:
                simulate_region(f, contig, curr, length, bg_out_prefix)
        f.close()

    # 3. Simulate Loci (WT and Mutated)
    print("Simulating Loci...")
    
    for i, meta in enumerate(disease_metadata):
        if not meta: continue
        
        name = meta['name'].replace(" ", "_")
        
        # For this disease locus, we generate:
        # 1. WT reads (Hap1 & Hap2) -> For everyone else
        # 2. Mutated reads (Hap1 & Hap2 based on inheritance) -> For this patient
        
        wt_prefix = os.path.join(OUTPUT_DIR, f"locus_{i}_wt")
        mut_prefix = os.path.join(OUTPUT_DIR, f"locus_{i}_mut")
        
        # Clear
        for p in [wt_prefix, mut_prefix]:
            open(f"{p}1.fq", "w").close()
            open(f"{p}2.fq", "w").close()

        # Handle Hap1
        if meta['hap1_map']:
            m = meta['hap1_map']
            # Simulate WT
            with pysam.FastaFile(HAP1_FASTA) as f:
                simulate_region(f, m['contig'], m['roi_start'], m['roi_end'], wt_prefix)
            
            # Simulate Mutated (if dominant or recessive)
            # Patch in memory
            with pysam.FastaFile(HAP1_FASTA) as f:
                seq = f.fetch(m['contig'], m['roi_start'], m['roi_end'])
                # Calc relative pos
                # m['pos'] is 1-based global. m['roi_start'] is 0-based global.
                rel = (m['pos'] - 1) - m['roi_start']
                
                # Prepare patch
                p_ref = meta['ref']
                p_alt = meta['alt']
                if m['strand'] == '-':
                    p_ref = rc(p_ref)
                    p_alt = rc(p_alt)
                
                # Check match
                chunk_ref = seq[rel:rel+len(p_ref)]
                # Apply
                seq_mut = seq[:rel] + p_alt + seq[rel+len(p_ref):]
                
                # Write to temp fasta for simulation (since our helper takes file/contig)
                # Or refactor helper. Actually helper takes fasta_obj. 
                # We can't use helper easily for string.
                # Let's just write to temp file
                tmp_mut = os.path.join(OUTPUT_DIR, "tmp_mut.fa")
                with open(tmp_mut, "w") as tf:
                    tf.write(f">mut\n{seq_mut}\n")
                
                with pysam.FastaFile(tmp_mut) as tf:
                    simulate_region(tf, "mut", 0, len(seq_mut), mut_prefix)

        # Handle Hap2
        if meta['hap2_map']:
            m = meta['hap2_map']
            # Simulate WT
            with pysam.FastaFile(HAP2_FASTA) as f:
                simulate_region(f, m['contig'], m['roi_start'], m['roi_end'], wt_prefix)
                
            # Simulate Mutated (ONLY if Recessive)
            is_recessive = (meta['inheritance'] == 'recessive')
            
            if is_recessive:
                with pysam.FastaFile(HAP2_FASTA) as f:
                    seq = f.fetch(m['contig'], m['roi_start'], m['roi_end'])
                    rel = (m['pos'] - 1) - m['roi_start']
                    p_ref = meta['ref']
                    p_alt = meta['alt']
                    if m['strand'] == '-':
                        p_ref = rc(p_ref)
                        p_alt = rc(p_alt)
                    seq_mut = seq[:rel] + p_alt + seq[rel+len(p_ref):]
                    
                    tmp_mut = os.path.join(OUTPUT_DIR, "tmp_mut2.fa")
                    with open(tmp_mut, "w") as tf:
                        tf.write(f">mut\n{seq_mut}\n")
                    with pysam.FastaFile(tmp_mut) as tf:
                        simulate_region(tf, "mut", 0, len(seq_mut), mut_prefix)
            else:
                # If dominant, Hap2 is WT. 
                # The 'mut_prefix' file should contain WT reads for Hap2?
                # Yes, because for the "Patient", they need Mutated Hap1 AND WT Hap2.
                # The 'mut_prefix' represents "The reads for this locus for the patient".
                with pysam.FastaFile(HAP2_FASTA) as f:
                     simulate_region(f, m['contig'], m['roi_start'], m['roi_end'], mut_prefix)

    # 4. Generate Construction Recipes
    print("Generating assembly scripts...")
    script_path = os.path.join(OUTPUT_DIR, "assemble_samples.sh")
    with open(script_path, "w") as sh:
        sh.write("#!/bin/bash\n")
        
        for i, meta in enumerate(disease_metadata):
            if not meta: continue
            name = meta['name'].replace(" ", "_")
            
            # Construct File Lists
            # Start with background
            files_1 = ["background1.fq"]
            files_2 = ["background2.fq"]
            
            for j, other_meta in enumerate(disease_metadata):
                if not other_meta: continue
                
                if i == j:
                    # Use MUTATED for this locus
                    files_1.append(f"locus_{j}_mut1.fq")
                    files_2.append(f"locus_{j}_mut2.fq")
                else:
                    # Use WT for all other loci
                    files_1.append(f"locus_{j}_wt1.fq")
                    files_2.append(f"locus_{j}_wt2.fq")
            
            # Command
            cmd1 = f"cat {' '.join(files_1)} > {name}_final_1.fq"
            cmd2 = f"cat {' '.join(files_2)} > {name}_final_2.fq"
            
            sh.write(f"echo 'Assembling {name}...'\n")
            sh.write(cmd1 + "\n")
            sh.write(cmd2 + "\n")
            
    print(f"Master simulation complete. Run {script_path} to build samples.")

if __name__ == "__main__":
    main()
