import subprocess
import shutil
import random
import pysam
import os

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def simulate_reads_python(fasta_path: str, output_prefix: str, coverage: int, read_length: int, paired_end: bool = True):
    print(f"Using Python fallback simulator for {fasta_path}...")
    fasta = pysam.FastaFile(fasta_path)
    total_length = sum(fasta.lengths)
    
    # Calculate number of reads
    # Coverage = (N_reads * L_read) / Genome_Size
    # N_reads = (C * G) / L
    total_reads_needed = int((coverage * total_length) / read_length)
    if paired_end:
        num_pairs = total_reads_needed // 2
    else:
        num_pairs = total_reads_needed # variable name reuse for simplicity

    print(f"Genome length: {total_length}. Target coverage: {coverage}x. Generating {num_pairs} pairs.")

    # Open output files
    suffix1 = "1.fq"
    suffix2 = "2.fq"
    
    # ART uses .fq extension usually
    # If output_prefix ends in a dot or similar, handle it? 
    # Usually output_prefix is "dir/sample". 
    # ART outputs "sample1.fq" and "sample2.fq"
    
    fq1_path = f"{output_prefix}1.fq"
    fq2_path = f"{output_prefix}2.fq"

    with open(fq1_path, 'w') as fq1, open(fq2_path, 'w') as fq2:
        for i in range(num_pairs):
            # Pick random contig weighted by length? 
            # Simple approach: Pick random contig from list, verify length > fragment size. 
            # Better: Cumulative distribution. But efficient enough for now:
            
            # Use random.choices with weights once?
            # Doing it every loop is slow. 
            # Let's just iterate contigs or pick random index?
            # For 3GB genome, maybe just pick contig index based on weights.
            contig_name = random.choices(fasta.references, weights=fasta.lengths, k=1)[0]
            contig_len = fasta.get_reference_length(contig_name)
            
            fragment_size = int(random.gauss(500, 50)) # Mean 500, SD 50
            if fragment_size < read_length: fragment_size = read_length + 10
            if fragment_size >= contig_len: 
                 # Skip short contigs or retry
                 continue

            start_pos = random.randint(0, contig_len - fragment_size)
            
            # Fetch sequence
            # fasta.fetch is 0-based, end exclusive
            fragment = fasta.fetch(contig_name, start_pos, start_pos + fragment_size)
            
            if len(fragment) < fragment_size: continue # Should not happen

            # R1: start of fragment
            seq1 = fragment[:read_length]
            # R2: end of fragment, reverse complemented
            seq2_raw = fragment[-read_length:]
            seq2 = reverse_complement(seq2_raw)

            qual = 'I' * read_length # High quality score

            # Write FASTQ
            read_name = f"@{contig_name}-{i}:{start_pos}"
            
            fq1.write(f"{read_name}/1\n{seq1}\n+\n{qual}\n")
            if paired_end:
                fq2.write(f"{read_name}/2\n{seq2}\n+\n{qual}\n")
            
            if i % 10000 == 0:
                print(f"Generated {i}/{num_pairs} pairs...", end='\r')

    print(f"\nSimulation complete. Output: {fq1_path}, {fq2_path}")


def simulate_reads(fasta_path: str, output_prefix: str, coverage: int, read_length: int, paired_end: bool = True, simulator_path: str = 'art_illumina'):
    """
    Simulates Illumina reads from a FASTA file using ART.

    Args:
        fasta_path: Path to the genome FASTA file.
        output_prefix: Prefix for the output FASTQ files.
        coverage: The desired read coverage.
        read_length: The length of the reads.
        paired_end: Whether to generate paired-end reads.
        simulator_path: Path to the 'art_illumina' executable.
    """
    if shutil.which(simulator_path):
        cmd = [
            simulator_path,
            '-i', fasta_path,
            '-o', output_prefix,
            '-l', str(read_length),
            '-f', str(coverage),
        ]

        if paired_end:
            cmd.extend(['-p', '-m', '500', '-s', '10']) # Default mean fragment size and std dev

        print(f"Running ART command: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("ART simulation completed successfully.")
            print(result.stdout)
            if result.stderr:
                print("ART stderr:")
                print(result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error running ART: {e}")
            print(f"ART stdout: {e.stdout}")
            print(f"ART stderr: {e.stderr}")
            raise
    else:
        print(f"'{simulator_path}' not found. Falling back to Python simulator.")
        simulate_reads_python(fasta_path, output_prefix, coverage, read_length, paired_end)

if __name__ == '__main__':
    # This is for testing purposes.
    # We will need to create a dummy fasta file for this.
    pass
