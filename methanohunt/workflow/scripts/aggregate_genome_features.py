import os
import sys
import argparse
import pandas as pd
import glob
from natsort import natsort_keygen

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--taxonomy", default=None)
    parser.add_argument("--name", default=None)
    parser.add_argument("--col", default=None)
    parser.add_argument("--output", required=True)
    return parser.parse_args()

def main():
    args = get_args()
    
    gene_db_path = os.path.join(args.db, "methanohunt_gene_db.tsv")
    tax_db_path = os.path.join(args.db, "methanohunt_taxonomy_db.tsv")
    
    if not os.path.exists(gene_db_path):
        print(f"Error: Database file not found at {gene_db_path}", file=sys.stderr)
        sys.exit(1)
        
    db_df = pd.read_csv(gene_db_path, sep="\t")
    # Dereplicate KO values as requested due to the new Methanohunt_classification column causing duplicate KOs
    db_df = db_df.drop_duplicates(subset=["KO"])
    
    ko_files = [f for f in glob.glob(os.path.join(args.outdir, "KofamScan_anno", "*_ko.tsv")) if not f.endswith("_raw_ko.tsv")]
    genomes = []
    
    classification_map = {
        "McrA (methanogen)": "Methanogen",
        "McrA (ANME)": "Anaerobic methanotroph (ANME)",
        "McrA (ANKA)": "Anaerobic alkanotroph (ANKA)",
        "PmoA": "MOB (w/ pMMO)",
        "AmoA": "Ammonia oxidizer",
        "HMO": "Hydrocarbon oxidizer",
        "MmoX": "MOB (w/ sMMO)",
        "PrmA": "Propane oxidizer",
        "MimA": "Alkane oxidizer"
    }

    results = []

    for kofam_f in ko_files:
        basename = os.path.basename(kofam_f)
        genome = basename.replace("_ko.tsv", "")
        
        # Load KO
        try:
            ko_df = pd.read_csv(kofam_f, sep="\t", header=0)
            if len(ko_df.columns) == 4 and "kofam_KO" in ko_df.columns:
                 # In case kofamscan_filter output has specific names, but user said "genome, gene_id, ko"
                 pass
            # Just to be safe, rename irst 3 cols if they lack "ko" verbatim
            # Often it's `genome`, `gene_id`, `ko`, `evalue`
            if "ko" not in ko_df.columns and "KO" in ko_df.columns:
                ko_df.rename(columns={"KO": "ko"}, inplace=True)
            elif "ko" not in ko_df.columns and len(ko_df.columns) >= 3:
                 ko_df.rename(columns={ko_df.columns[2]: "ko"}, inplace=True)
                 ko_df.rename(columns={ko_df.columns[1]: "gene_id"}, inplace=True)
        except Exception:
            ko_df = pd.DataFrame(columns=["genome", "gene_id", "ko"])
            
        gene_f = os.path.join(args.outdir, "cust_db_anno", f"{genome}_methanohunt_gene_classification.tsv")
        try:
            gene_df = pd.read_csv(gene_f, sep="\t", header=0)
        except Exception:
            gene_df = pd.DataFrame(columns=["gene_id", "classification"])

        # Collect direct classifications from gene module (HMM)
        hmm_classes = set()
        if "classification" in gene_df.columns:
            for c in gene_df["classification"].dropna().unique():
                 if c in classification_map:
                     hmm_classes.add(classification_map[c])

        # Priority 2: Filter KO vs Database
        # If strict, filter db_df first
        # But wait, strict applies to evaluating if the genome has the enzyme.
        # So we can merge first.
        
        merged_db = db_df.copy()
        if args.strict:
            # Drop rows where If_key_in_enzyme != TRUE
            # Handling float/str TRUES
            merged_db = merged_db[merged_db["If_key_in_enzyme"].astype(str).str.upper() == "TRUE"]

        ko_hits = []
        if "ko" in ko_df.columns:
             ko_list = ko_df["ko"].dropna().unique()
             matching_rows = merged_db[merged_db["KO"].isin(ko_list)]
        else:
             matching_rows = pd.DataFrame(columns=merged_db.columns)
             ko_list = []
             
        ko_classes = set()
        
        # Priority 2 checks if hmm_classes is empty
        if not hmm_classes:
             if "Methanohunt_classification" in matching_rows.columns:
                 for c in matching_rows["Methanohunt_classification"].dropna().unique():
                     if c in classification_map:
                         ko_classes.add(classification_map[c])

        final_classes = list(hmm_classes) if hmm_classes else list(ko_classes)
        final_class_str = "; ".join(sorted(final_classes)) if final_classes else ""

        # Subgroup definitions
        subgroups = set()
        is_methanogen = "Methanogen" in final_classes or "Methanogen (putative, mcrA absent)" in final_classes
        if is_methanogen:
            c_found = matching_rows["Classification"].dropna()
            e_found = set([str(x).lower() for x in matching_rows["Enzyme"].dropna()])
            
            # Note: "Acetoclastic" or "Hydrogenotrophic" can be part of strings like "Hydrogenotrophic & Acetotrophic"
            # It's safer to use str.contains
            if c_found.str.contains("Hydrogenotrophic").any():
                subgroups.add("Hydrogenotrophic methanogen")
            if c_found.str.contains("Acetoclastic").any():
                subgroups.add("Acetotrophic methanogen")
            if "mta" in e_found:
                subgroups.add("Methylotrophic methanogen (Methanol)")
                
            has_mtbA = "mtba" in e_found
            if has_mtbA:
                if "mtm" in e_found:
                    subgroups.add("Methylotrophic methanogen (Methylamine)")
                if "mtb" in e_found:
                    subgroups.add("Methylotrophic methanogen (Dimethylamine)")
                if "mtt" in e_found:
                    subgroups.add("Methylotrophic methanogen (Trimethylamine)")
                if not any(e in e_found for e in ["mtm", "mtb", "mtt"]):
                    subgroups.add("Methylotrophic methanogen (Unknown substrate)")
                    
        final_subgroups_str = "; ".join(sorted(list(subgroups)))
        
        results.append({
            "Genome": genome,
            "Classification by gene (methanohunt)": final_class_str,
            "Subgroup by gene (methanohunt)": final_subgroups_str
        })

    out_df = pd.DataFrame(results)

    # Taxonomy Mapping
    if args.taxonomy and args.col:
        try:
            tax_user_df = pd.read_csv(args.taxonomy, sep="\t")
            tax_user_df.rename(columns=lambda x: str(x).strip('\ufeff').strip(), inplace=True)  # remove BOM
            if args.col not in tax_user_df.columns:
                 print(f"Warning: Taxonomy column '{args.col}' not found in {args.taxonomy}. Skipping taxonomy classification.")
            else:
                 if not args.name:
                      print(f"Warning: Taxonomy matching requires --name argument. Skipping taxonomy classification.")
                 elif args.name not in tax_user_df.columns:
                      print(f"Warning: Name column '{args.name}' not found in {args.taxonomy}. Available columns: {list(tax_user_df.columns)}. Skipping taxonomy classification.")
                 else:
                      merge_col = args.name
                      
                      out_df = out_df.merge(tax_user_df[[merge_col, args.col]], left_on="Genome", right_on=merge_col, how="left")
                 if merge_col != "Genome":
                      out_df.drop(columns=[merge_col], inplace=True)
                      
                 # Process taxonomy against DB
                 tax_db = pd.read_csv(tax_db_path, sep="\t")
                 valid_tax_db = tax_db.dropna(subset=['GTDB_taxonomy']).copy()
                 valid_tax_db['GTDB_taxonomy_clean'] = valid_tax_db['GTDB_taxonomy'].apply(lambda x: "".join(str(x).split()))
                 valid_tax_db['tax_len'] = valid_tax_db['GTDB_taxonomy_clean'].str.len()
                 tax_db_sorted = valid_tax_db.sort_values(by='tax_len', ascending=False)
                 
                 tax_classes = []
                 tax_subs = []
                 for val in out_df[args.col]:
                     matched = False
                     if pd.notna(val):
                         val_str_clean = "".join(str(val).split())
                         for _, row in tax_db_sorted.iterrows():
                             if row['GTDB_taxonomy_clean'] in val_str_clean:
                                 tax_classes.append(row.get('Classification', None))
                                 tax_subs.append(row.get('Subgroup', None))
                                 matched = True
                                 break
                     if not matched:
                         tax_classes.append("")
                         tax_subs.append("")
                         
                 out_df["Classification by taxonomy (methanohunt)"] = tax_classes
                 out_df["Subgroup by taxonomy (methanohunt)"] = tax_subs
                 # Remove the intermediate taxonomy string column 
                 out_df.drop(columns=[args.col], inplace=True)
        except Exception as e:
             print(f"Warning: Failed to process taxonomy classification: {e}")
    else:
         # Empty columns if not requested (as per example)
         out_df["Classification by taxonomy (methanohunt)"] = ""
         out_df["Subgroup by taxonomy (methanohunt)"] = ""

    cols_to_check = [c for c in out_df.columns if c != "Genome"]
    out_df[cols_to_check] = out_df[cols_to_check].replace(r'^\s*$', pd.NA, regex=True)
    out_df.dropna(subset=cols_to_check, how='all', inplace=True)
    out_df.fillna("", inplace=True)

    out_df.sort_values(by="Genome", key=natsort_keygen(), inplace=True)

    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"Summary saved to {args.output}")

if __name__ == "__main__":
    main()
