import os
import sys
import json
import subprocess
import pysam
import shutil
import tempfile
from collections import defaultdict

# Ensure src is in path
base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(base_dir, "src"))

# Config
HAP1_FASTA = os.path.join(base_dir, "data", "GCA_050492415.1_apr041.1_v1_genomic.fna")
HAP2_FASTA = os.path.join(base_dir, "data", "GCA_050492395.1_apr041.2_v1_genomic.fna")
OUTPUT_DIR = os.path.join(base_dir, "output", "master_sim")
CONFIG_PATH = os.path.join(base_dir, "config", "diseases.json")
HG38_FASTA = os.path.join(base_dir, "hg38.fa")
ROI_PADDING = 10000 # +/- 10kb
ART_PATH = "art_illumina" # Assumes it is in PATH

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
        return mapping['t_end'] - offset 

def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def simulate_region(fasta_obj, contig, start, end, final_out_prefix, coverage=15, read_len=150):
    """
    Simulates reads for a specific region [start, end) of a contig using ART.
    Writes to a temporary file and then appends to final_out_prefix.
    """
    # Just extract sequence
    seq = fasta_obj.fetch(contig, start, end)
    length = len(seq)
    if length < read_len + 50: return # Too short for ART
    
    # Create temp fasta
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp_fa:
        tmp_fa.write(f">{contig}:{start}-{end}\n{seq}\n")
        tmp_fa_path = tmp_fa.name
        
    # Temp output prefix
    with tempfile.NamedTemporaryFile(delete=False) as tmp_out:
        tmp_out_prefix = tmp_out.name
    
    # Run ART
    cmd = [
        ART_PATH,
        '-i', tmp_fa_path,
        '-o', tmp_out_prefix,
        '-l', str(read_len),
        '-f', str(coverage),
        '-p', '-m', '500', '-s', '10',
        '-q' # Quiet
    ]
    
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        # ART output files
        art_fq1 = f"{tmp_out_prefix}1.fq"
        art_fq2 = f"{tmp_out_prefix}2.fq"
        
        # Append to final
        # Note: final_out_prefix implies we write to {final_out_prefix}1.fq and {final_out_prefix}2.fq
        if os.path.exists(art_fq1) and os.path.exists(art_fq2):
            with open(f"{final_out_prefix}1.fq", "ab") as f_out, open(art_fq1, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)
            with open(f"{final_out_prefix}2.fq", "ab") as f_out, open(art_fq2, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)
            
        # Cleanup ART outputs
        if os.path.exists(art_fq1): os.remove(art_fq1)
        if os.path.exists(art_fq2): os.remove(art_fq2)
            
    except subprocess.CalledProcessError:
        print(f"Warning: ART simulation failed for {contig}:{start}-{end}")
    except Exception as e:
        print(f"Error simulating {contig}: {e}")
    finally:
        if os.path.exists(tmp_fa_path): os.remove(tmp_fa_path)
        if os.path.exists(tmp_out_prefix): os.remove(tmp_out_prefix)


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
    print("Simulating Background (this may take time)...")
    
    bg_out_prefix = os.path.join(OUTPUT_DIR, "background")
    # Clear prev
    open(f"{bg_out_prefix}1.fq", "w").close()
    open(f"{bg_out_prefix}2.fq", "w").close()

    for hap_label, fasta_path in [('hap1', HAP1_FASTA), ('hap2', HAP2_FASTA)]:
        print(f"Processing {hap_label} from {fasta_path}...")
        f = pysam.FastaFile(fasta_path)
        
        # Iterate ALL contigs
        for contig in f.references:
            length = f.get_reference_length(contig)
            excls = exclusions[hap_label].get(contig, [])
            
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
            
            # Simulate Mutated
            with pysam.FastaFile(HAP1_FASTA) as f:
                seq = f.fetch(m['contig'], m['roi_start'], m['roi_end'])
                rel = (m['pos'] - 1) - m['roi_start']
                p_ref = meta['ref']
                p_alt = meta['alt']
                if m['strand'] == '-':
                    p_ref = rc(p_ref)
                    p_alt = rc(p_alt)
                
                seq_mut = seq[:rel] + p_alt + seq[rel+len(p_ref):]
                
                # Write to temp fasta for simulation
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tf:
                    tf.write(f">mut\n{seq_mut}\n")
                    tf_path = tf.name
                
                with pysam.FastaFile(tf_path) as tf:
                    simulate_region(tf, "mut", 0, len(seq_mut), mut_prefix)
                
                os.remove(tf_path)

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
                    
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tf:
                        tf.write(f">mut\n{seq_mut}\n")
                        tf_path = tf.name
                    with pysam.FastaFile(tf_path) as tf:
                        simulate_region(tf, "mut", 0, len(seq_mut), mut_prefix)
                    os.remove(tf_path)
            else:
                # If dominant, Hap2 is WT.
                with pysam.FastaFile(HAP2_FASTA) as f:
                     simulate_region(f, m['contig'], m['roi_start'], m['roi_end'], mut_prefix)

    # 4. Generate Construction Recipes
    print("Generating assembly scripts...")
    script_path = os.path.join(OUTPUT_DIR, "assemble_samples.sh")
    with open(script_path, "w") as sh:
        sh.write("#!/bin/bash\n")
        sh.write("set -e\n") # Exit on error
        
        # Check for gzip/bgzip
        sh.write("ZIP=gzip\n")
        sh.write("if command -v bgzip &> /dev/null; then ZIP=bgzip; fi\n\n")

        for i, meta in enumerate(disease_metadata):
            if not meta: continue
            name = meta['name'].replace(" ", "_")
            
            # Construct File Lists
            files_1 = ["background1.fq"]
            files_2 = ["background2.fq"]
            
            for j, other_meta in enumerate(disease_metadata):
                if not other_meta: continue
                
                if i == j:
                    files_1.append(f"locus_{j}_mut1.fq")
                    files_2.append(f"locus_{j}_mut2.fq")
                else:
                    files_1.append(f"locus_{j}_wt1.fq")
                    files_2.append(f"locus_{j}_wt2.fq")
            
            # Command
            cmd1 = f"cat {' '.join(files_1)} | $ZIP > {name}_final_1.fq.gz"
            cmd2 = f"cat {' '.join(files_2)} | $ZIP > {name}_final_2.fq.gz"
            
            sh.write(f"echo 'Assembling {name}...'\n")
            sh.write(cmd1 + "\n")
            sh.write(cmd2 + "\n")
            
    print(f"Master simulation complete. Run 'bash {script_path}' to build samples.")

if __name__ == "__main__":
    main()
