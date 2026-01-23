import pandas as pd
from Bio import SeqIO
import os

def extract_sequences(classification_file, protein_file, output_files):
    # Read classification table
    df = pd.read_csv(classification_file, sep='\t')
    
    # Ensure classification column exists (it should be renamed to 'classification' by classify.smk)
    # The user request clarified that column 'classification' is used for filtering.
    # And 'name' column for ID matching.
    
    if 'subtype' not in df.columns or 'name' not in df.columns:
        # Fallback or error? Assuming structure is correct per previous steps.
        # If empty or not fully formed, we might just return empty files.
        pass

    # Create distinct sets of IDs for each category
    ids_methanogen = {}
    ids_anme = {}
    ids_pmoa = {}
    ids_mmox = {}

    for _, row in df.iterrows():
        # User requested to rename original 'classification' to 'subtype' and put original values there.
        # New 'classification' column has broader categories like 'McrA (methanogen)'.
        # Our filtering logic relies on "Methanogen...", "ANME..." which are in 'subtype' now.
        subtype = str(row.get('subtype', ''))
        seq_id = str(row.get('name', ''))
        
        if not subtype or not seq_id:
            continue

        if subtype.startswith('Methanogen'):
            ids_methanogen[seq_id] = subtype
        elif subtype.startswith('ANME'):
            ids_anme[seq_id] = subtype
        elif subtype.startswith('PmoA'):
            ids_pmoa[seq_id] = subtype
        elif subtype.startswith('MmoX'):
            ids_mmox[seq_id] = subtype

    # Open output files
    handle_methanogen = open(output_files['mcrA_methanogen'], 'w')
    handle_anme = open(output_files['mcrA_anme'], 'w')
    handle_pmoa = open(output_files['pmoA'], 'w')
    handle_mmox = open(output_files['mmoX'], 'w')

    try:
        # Parse protein FASTA
        for record in SeqIO.parse(protein_file, "fasta"):
            # logic to match ID. Usually headers are like ">ID description". record.id gets ID.
            # gappa/pplacer usually preserves the ID.
            
            # Check if likely a hit we want
            # Only exact match?
            
            # Modify header/description
            # User wants: "header保持原来的ID（first part），然后添加一个空格，空格之后要放上classification中的信息"
            # record.id is usually the first part. record.description is the full line.
            # We construct new description: "{id} {classification}"
            
            original_id = record.id
            
            if original_id in ids_methanogen:
                record.description = f"{original_id} {ids_methanogen[original_id]}"
                SeqIO.write(record, handle_methanogen, "fasta")
            elif original_id in ids_anme:
                record.description = f"{original_id} {ids_anme[original_id]}"
                SeqIO.write(record, handle_anme, "fasta")
            elif original_id in ids_pmoa:
                record.description = f"{original_id} {ids_pmoa[original_id]}"
                SeqIO.write(record, handle_pmoa, "fasta")
            elif original_id in ids_mmox:
                record.description = f"{original_id} {ids_mmox[original_id]}"
                SeqIO.write(record, handle_mmox, "fasta")
                
    finally:
        handle_methanogen.close()
        handle_anme.close()
        handle_pmoa.close()
        handle_mmox.close()

if __name__ == "__main__":
    # Snakemake input/output
    classification_file = snakemake.input.classification
    protein_file = snakemake.input.protein
    
    output_files = {
        'mcrA_methanogen': snakemake.output.mcrA_methanogen,
        'mcrA_anme': snakemake.output.mcrA_anme,
        'pmoA': snakemake.output.pmoA,
        'mmoX': snakemake.output.mmoX
    }

    extract_sequences(classification_file, protein_file, output_files)
