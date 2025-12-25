import argparse
import os
import pysam
from genome_patcher import patch_genome
from read_simulator import simulate_reads
from phenotype import create_phenopacket

def main():
    parser = argparse.ArgumentParser(description="A tool to generate synthetic genomic data with pathogenic variants.")
    parser.add_argument("--fasta", required=True, help="Path to the input FASTA file.")
    parser.add_argument("--fasta2", help="Path to the second input FASTA file (for the second haplotype).")
    
    # Replaced VCF args with direct variant args
    parser.add_argument("--chrom", required=True, help="Chromosome name (e.g., chr4).")
    parser.add_argument("--pos", type=int, required=True, help="1-based position of the variant.")
    parser.add_argument("--ref", required=True, help="Reference allele.")
    parser.add_argument("--alt", required=True, help="Alternate allele.")
    parser.add_argument("--omim-id", required=True, help="OMIM ID for phenotype generation (e.g., 143100).")

    parser.add_argument("--output-dir", required=True, help="The directory to write the output files to.")
    parser.add_argument("--sample-id", required=True, help="The ID for the sample.")
    parser.add_argument("--coverage", type=int, default=30, help="Read coverage for simulation.")
    parser.add_argument("--read-length", type=int, default=150, help="Read length for simulation.")
    parser.add_argument("--haplotype", default='both', choices=['hap1', 'hap2', 'both'], help="Which haplotype to patch.")
    parser.add_argument("--art-path", default='art_illumina', help="Path to the art_illumina executable.")

    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # --- 1. Handle FASTA inputs ---
    input_fasta = args.fasta
    if args.fasta2:
        print("Two FASTA files provided. Concatenating them.")
        concatenated_fasta_path = os.path.join(args.output_dir, "concatenated_genome.fa")
        with open(concatenated_fasta_path, 'wb') as outfile:
            with open(args.fasta, 'rb') as infile:
                outfile.write(infile.read())
            with open(args.fasta2, 'rb') as infile:
                outfile.write(infile.read())
        input_fasta = concatenated_fasta_path

    # --- 2. Define output paths ---
    patched_fasta_path = os.path.join(args.output_dir, f"{args.sample_id}.patched.fa")
    read_output_prefix = os.path.join(args.output_dir, f"{args.sample_id}")
    phenopacket_path = os.path.join(args.output_dir, f"{args.sample_id}.phenopacket.json")

    # --- 3. Patch genome ---
    print("\n--- Patching Genome ---")
    patch_genome(
        fasta_path=input_fasta,
        chrom=args.chrom,
        pos=args.pos,
        ref=args.ref,
        alt=args.alt,
        output_fasta_path=patched_fasta_path,
        haplotype=args.haplotype
    )

    # --- 4. Simulate reads ---
    print("\n--- Simulating Reads ---")
    simulate_reads(
        fasta_path=patched_fasta_path,
        output_prefix=read_output_prefix,
        coverage=args.coverage,
        read_length=args.read_length,
        simulator_path=args.art_path
    )

    # --- 5. Create Phenopacket ---
    print("\n--- Creating Phenopacket ---")
    create_phenopacket(
        omim_id=args.omim_id,
        sample_id=args.sample_id,
        output_phenopacket_path=phenopacket_path
    )

    print("\nWorkflow finished successfully!")

if __name__ == "__main__":
    main()
