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

    # Identify matching contigs
    # We look for contigs that contain the chrom name
    matching_contigs = sorted([c for c in fasta_file.references if chrom in c])
    
    # Special handling for exact matches if available
    if chrom in fasta_file.references:
        # If the exact chrom name exists, prioritize it (or maybe it's the only one we want?)
        # For now, let's stick to the substring logic as it handles the 'hap1' suffix case well,
        # but if we have 'chr4' and 'chr4_hap1', we might get both. 
        # If the user provides 'chr4', and we have 'chr4' in fasta, we probably want that.
        if chrom in matching_contigs:
             # Refine: if we have exact match, maybe we don't want others? 
             # But the requirement is to handle the hap1/hap2 contigs.
             pass

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
            f_out.write(seq)
            f_out.write("\n")

    fasta_file.close()

if __name__ == '__main__':
    # This is for testing purposes.
    # We will need to create dummy data for this.
    pass
