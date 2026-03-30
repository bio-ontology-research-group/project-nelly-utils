import pysam

def patch_genome(fasta_path: str, chrom: str, pos: int, ref: str, alt: str, output_fasta_path: str, haplotype: str = 'both'):
    """
    Patches a reference genome with a specific variant provided directly.

    Args:
        fasta_path: Path to the reference genome FASTA file.
        chrom: The chromosome name (e.g., 'chr4', '4').
        pos: The 1-based position of the variant.
        ref: The reference allele.
        alt: The alternate allele.
        output_fasta_path: Path to write the patched FASTA file.
        haplotype: Which haplotype to patch ('hap1', 'hap2', or 'both').
                   Assumes that for a given variant chrom, there are two matching contigs in the FASTA.
                   They are sorted alphabetically: first is hap1, second is hap2.
    """
    # Open the reference genome
    fasta_file = pysam.FastaFile(fasta_path)

    # Validate Variant (Simple SNP check for now, can be expanded)
    if not (len(ref) == 1 and len(alt) == 1):
        print(f"Warning: Variant {chrom}:{pos} {ref}>{alt} is not a SNP. Indel support is experimental.")
        # Proceeding with experimental support (the logic below handles string replacement, so it might just work for simple indels)

    # Adjust to 0-based
    pos_0 = pos - 1

    # Identify matching contigs — prefer exact match, fall back to prefix match
    if chrom in fasta_file.references:
        matching_contigs = [chrom]
    else:
        # Fall back to prefix match (e.g., "chr4" matches "chr4_ctg9_hap1")
        matching_contigs = sorted([c for c in fasta_file.references if c.startswith(chrom + "_") or c == chrom])

    if not matching_contigs:
        print(f"Warning: No contigs matching '{chrom}' found in FASTA. Skipping variant.")
    
    # Map contigs to patch decision
    contigs_to_patch = set()
    if len(matching_contigs) == 1:
        if haplotype != 'both': 
             print(f"Info: Only one matching contig found for {chrom}. Patching it.")
        contigs_to_patch.add(matching_contigs[0])
        
    elif len(matching_contigs) >= 2:
        # Assume sorted: hap1, hap2
        hap1_contig = matching_contigs[0]
        hap2_contig = matching_contigs[1]
        
        if haplotype in ('hap1', 'both'):
            contigs_to_patch.add(hap1_contig)
        if haplotype in ('hap2', 'both'):
            contigs_to_patch.add(hap2_contig)

    # Process and write output
    # We stream read/write to avoid high memory usage
    with open(output_fasta_path, "w") as f_out:
        for contig in fasta_file.references:
            seq = fasta_file.fetch(contig)
            
            if contig in contigs_to_patch:
                # Apply patch
                # Check bounds
                if pos_0 >= len(seq):
                     print(f"Error: Position {pos} out of bounds for {contig} (len {len(seq)}). Skipping.")
                elif seq[pos_0:pos_0+len(ref)].upper() != ref.upper():
                    print(f"Warning: Ref mismatch at {contig}:{pos}. Expected {ref}, found {seq[pos_0:pos_0+len(ref)]}. Patching anyway.")
                    # Apply anyway
                    seq = seq[:pos_0] + alt + seq[pos_0+len(ref):]
                    print(f"Applied patch to {contig}: {pos} {ref}>{alt}")
                else:
                    seq = seq[:pos_0] + alt + seq[pos_0+len(ref):]
                    print(f"Applied patch to {contig}: {pos} {ref}>{alt}")

            f_out.write(f">{contig}\n")
            # Write sequence with 60-char line wrapping (required for samtools faidx)
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + "\n")

    fasta_file.close()

if __name__ == '__main__':
    # This is for testing purposes.
    # We will need to create dummy data for this.
    pass
