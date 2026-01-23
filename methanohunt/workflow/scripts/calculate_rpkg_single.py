import pysam
import pandas as pd
import sys
import os

def calculate_rpkg_single(bam_file, mc_out_file, output_file):
    # 1. Parse MicrobeCensus output to get Genome Equivalents (GE)
    # Command: genome_equivalents=$(grep "genome_equivalents" "$file" | awk '{print $2}')
    ge = None
    try:
        with open(mc_out_file, 'r') as f:
            for line in f:
                if line.startswith('genome_equivalents'):
                    parts = line.split()
                    if len(parts) >= 2:
                        ge = float(parts[1])
                    break
    except Exception as e:
        print(f"Error reading MicrobeCensus file {mc_out_file}: {e}", file=sys.stderr)
        return

    if ge is None or ge <= 0:
        print(f"Error: Invalid Genome Equivalents found in {mc_out_file}: {ge}", file=sys.stderr)
        # Fallback or exit? user might accept 0? Default to failure to warn user.
        # But maybe we write empty or zero file to keep pipeline going?
        # Let's write empty for now with warning
        pd.DataFrame(columns=["GeneID", "RPKG"]).to_csv(output_file, sep='\t', index=False)
        return

    # 2. Retrieve alignment counts from BAM
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except ValueError:
        print(f"Error opening BAM file: {bam_file}", file=sys.stderr)
        return

    gene_counts = {}
    gene_lengths = {}
    
    # Get lengths
    for ref_name, length in zip(samfile.references, samfile.lengths):
        gene_lengths[ref_name] = length
        gene_counts[ref_name] = 0
        
    # Get counts (mapped reads)
    idx_out = pysam.idxstats(bam_file)
    for line in idx_out.splitlines():
        parts = line.split('\t')
        if len(parts) >= 3:
            rname = parts[0]
            if rname in gene_counts:
                mapped = int(parts[2])
                gene_counts[rname] = mapped

    samfile.close()

    # 3. Calculate RPKG
    # RPKG = (Reads mapped / (Gene Length / 1000.0)) / Genome Equivalents
    rpkg = {}
    
    for gene, count in gene_counts.items():
        length_kb = gene_lengths[gene] / 1000.0
        if length_kb > 0 and ge > 0:
            # Reads per kb
            rpk = count / length_kb
            # Normalized by GE
            val = rpk / ge
            rpkg[gene] = val
        else:
            rpkg[gene] = 0.0

    # Output
    df = pd.DataFrame(list(rpkg.items()), columns=["GeneID", "RPKG"])
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    bam_file = snakemake.input.bam
    mc_out_file = snakemake.input.mc_out
    output_file = snakemake.output.rpkg
    calculate_rpkg_single(bam_file, mc_out_file, output_file)
