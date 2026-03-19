import pandas as pd
import os

def run_annotate(input_file, column_name, output_file, database_path=None):
    if not database_path:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        database_path = os.path.join(current_dir, "database", "methanohunt_taxonomy_db.tsv")
    
    if not os.path.exists(database_path):
        raise FileNotFoundError(f"Database file not found at {database_path}")
        
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found at {input_file}")

    db_df = pd.read_csv(database_path, sep="\t")
    
    # Pre-calculate lengths of GTDB_taxonomy and sort descending
    # to find the longest matching substring first
    valid_db = db_df.dropna(subset=['GTDB_taxonomy']).copy()
    valid_db['GTDB_taxonomy_clean'] = valid_db['GTDB_taxonomy'].apply(lambda x: "".join(str(x).split()))
    valid_db['tax_len'] = valid_db['GTDB_taxonomy_clean'].str.len()
    db_sorted = valid_db.sort_values(by='tax_len', ascending=False)
    
    # Load user data
    user_df = pd.read_csv(input_file, sep="\t")
    
    if column_name not in user_df.columns:
        raise ValueError(f"Column '{column_name}' not found in the input file '{input_file}'. Available columns are: {', '.join(user_df.columns)}")
    
    classifications = []
    subgroups = []
    match_count = 0
    
    for val in user_df[column_name]:
        matched = False
        if pd.notna(val):
            val_str_clean = "".join(str(val).split())
            for _, row in db_sorted.iterrows():
                if row['GTDB_taxonomy_clean'] in val_str_clean:
                    classifications.append(row.get('Classification', None))
                    subgroups.append(row.get('Subgroup', None))
                    match_count += 1
                    matched = True
                    break
        if not matched:
            classifications.append(None)
            subgroups.append(None)
            
    # Add new columns to the end
    user_df['MethanoHunt_classification'] = classifications
    user_df['MethanoHunt_subgroup'] = subgroups
    
    # Save to output file
    user_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"Annotation completed: matched {match_count} rows out of {len(user_df)} total rows.")
