from Bio import SeqIO
import pandas as pd
import sys

def extract_dna(classification_file, nucl_file, output_file):
    # Load classification to get list of gene names/IDs
    try:
        df = pd.read_csv(classification_file, sep='\t')
        if 'gene_id' not in df.columns:
            # Fallback or empty
            valid_ids = set()
        else:
            valid_ids = set(df['gene_id'].astype(str))
    except:
        valid_ids = set()
    
    if not valid_ids:
        # Create empty file
        open(output_file, 'w').close()
        return

    # Filter nucl file
    # Assuming IDs in nucl file match 'gene_id' in classification (or first part)
    # Using SeqIO index might be slow for huge files, straightforward parse is usually fine.
    
    with open(output_file, 'w') as out_handle:
        for record in SeqIO.parse(nucl_file, "fasta"):
            if record.id in valid_ids:
                SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":
    classification_file = snakemake.input.classification
    nucl_file = snakemake.input.nucl
    output_file = snakemake.output.ref
    
    extract_dna(classification_file, nucl_file, output_file)
