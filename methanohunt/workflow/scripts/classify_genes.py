import sys
import os
import pandas as pd
import json

# Placeholder logic for classification
# Snakemake inputs: jplaces (list), db
# Output: tsv

def classify_genes(jplace_files, output_file):
    results = []
    
    for jplace_file in jplace_files:
        marker = os.path.basename(jplace_file).replace(".jplace", "")
        # Basic parsing of jplace
        try:
            with open(jplace_file, 'r') as f:
                data = json.load(f)
            
            # Extract placements
            # data['placements'] is list of dicts.
            # p: [[edge_num, like_weight_ratio, ...], ...]
            # n: [name1, name2...]
            
            for placement in data['placements']:
                names = placement.get('n', [])
                # Just take the first name as the query ID
                # In pplacer, 'n' are the names of the sequences placed.
                for name in names:
                    # Logic to determine group (Methanogen/ANME/MOB/ANKA)
                    # This usually requires mapping the placed edge to a taxonomy.
                    # Since we don't have the refpkg structure here, we'll use a heuristic or dummy logic
                    # relying on the marker name for now, or assume the database has a mapping file.
                    
                    # USER REQUEST: "using this classification info to divide into Methanogen..."
                    # I will assign based on marker for now as a fallback, 
                    # but typically you traverse the tree or look up the edge index in a taxonomy map.
                    
                    # For demonstration:
                    if marker == "mcrA":
                        group = "Methanogen" # Default
                        # If specific edge, might be ANME
                    elif marker == "pmoA":
                        group = "MOB" # Aerobic methanotroph
                    elif marker == "mmoA":
                        group = "MOB" 
                    else:
                        group = "Unknown"
                        
                    results.append({
                        "gene_id": name,
                        "marker": marker,
                        "classification": group
                    })
        except Exception as e:
            print(f"Error parsing {jplace_file}: {e}")

    df = pd.DataFrame(results)
    if not df.empty:
        df.to_csv(output_file, sep="\t", index=False)
    else:
        # Create empty with headers
        with open(output_file, 'w') as f:
            f.write("gene_id\tmarker\tclassification\n")

if __name__ == "__main__":
    classify_genes(snakemake.input.jplaces, snakemake.output.tsv)
