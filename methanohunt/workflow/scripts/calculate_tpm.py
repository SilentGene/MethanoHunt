import sys
import argparse
import pandas as pd
import pysam

def calculate_tpm(bam_file, classification_file, output_file):
    # 1. Load classification to know which genes to track
    try:
        class_df = pd.read_csv(classification_file, sep="\t")
    except pd.errors.EmptyDataError:
        print("Classification file is empty. Creating empty abundance file.")
        with open(output_file, 'w') as f:
            f.write("Sample\tClassification\tDetails\n")
        return

    # Map gene_id to classification
    # class_df columns: gene_id, marker, classification
    gene_map = dict(zip(class_df['gene_id'], class_df['classification']))
    
    # 2. Parse BAM for counts and lengths
    # We use pysam to get reference lengths and map counts
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Initialize counts
    gene_counts = {}
    gene_lengths = {}
    
    # Get lengths from header
    for ref_name, length in zip(samfile.references, samfile.lengths):
        gene_lengths[ref_name] = length
        gene_counts[ref_name] = 0
        
    # Get counts (using idxstats for speed if indexed, else iterate)
    # idxstats returns: ref_name, seq_len, mapped, unmapped
    idxstats = samfile.get_index_statistics()
    # Pysam idxstats extraction might differ by version, but get_index_statistics returns records
    # Actually samfile.idxstats is a list of strings usually.
    # Let's use idxstats text parsing
    idx_out = pysam.idxstats(bam_file)
    for line in idx_out.splitlines():
        parts = line.split('\t')
        if len(parts) >= 3:
            rname = parts[0]
            # rlen = int(parts[1])
            mapped = int(parts[2])
            gene_counts[rname] = mapped

    samfile.close()
    
    # 3. Calculate TPM
    # RPK = count / (length/1000)
    # Scaling factor = sum(RPK) / 1e6
    # TPM = RPK / Scaling factor
    
    rpk = {}
    total_rpk = 0
    
    for gene, count in gene_counts.items():
        length = gene_lengths.get(gene, 0)
        if length > 0:
            r = count / (length / 1000.0)
            rpk[gene] = r
            total_rpk += r
            
    scaling_factor = total_rpk / 1000000.0
    
    tpm = {}
    if scaling_factor > 0:
        for gene, r in rpk.items():
            tpm[gene] = r / scaling_factor
    else:
        for gene in rpk:
            tpm[gene] = 0

    # 4. Aggregate by classification
    # We only care about genes in gene_map
    
    # Columns: classification, tpm
    # We aggregate by group
    
    group_tpm = {}
    
    for gene, val in tpm.items():
        if gene in gene_map:
            group = gene_map[gene]
            group_tpm[group] = group_tpm.get(group, 0) + val
            
    # 5. Output
    # Format: Sample | Group1 | Group2 ... 
    # Or Tidy format?
    # Taxonomy module uses wide format: Subgroup, Classification, Sample1_Abundance...
    # Here we are gene-centric.
    # Let's output a tidy format: Classification, TPM
    
    # Also include Sample name (derived from BAM filename or generic)
    sample_name = "Sample_1" # Placeholder
    
    results = []
    for group, val in group_tpm.items():
        results.append({"Classification": group, "TPM": val, "Sample": sample_name})
        
    df = pd.DataFrame(results)
    if not df.empty:
        df.to_csv(output_file, sep="\t", index=False)
    else:
        # Write header only
        with open(output_file, 'w') as f:
            f.write("Classification\tTPM\tSample\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    
    calculate_tpm(args.bam, args.classification, args.output)
