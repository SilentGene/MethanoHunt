import os
import sys
import argparse
import pandas as pd
import glob
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--pa", required=True, help="Output presence absence matrix")
    parser.add_argument("--count", required=True, help="Output count matrix")
    parser.add_argument("--count_tr", required=True, help="Output transposed count matrix")
    return parser.parse_args()

def main():
    args = get_args()
    
    db_path = os.path.join(args.db, "methanohunt_gene_db.tsv")
    if not os.path.exists(db_path):
        print(f"Error: Database file not found at {db_path}", file=sys.stderr)
        sys.exit(1)
        
    db_df = pd.read_csv(db_path, sep="\t")
    
    # Grab first 9 columns
    if len(db_df.columns) >= 9:
        cols_9 = list(db_df.columns[:9])
    else:
        cols_9 = list(db_df.columns)
        
    db_9cols_df = db_df[cols_9].copy()
    
    pa_df = db_9cols_df.copy()
    count_df = db_9cols_df.copy()

    ko_files = [f for f in glob.glob(os.path.join(args.outdir, "KofamScan_anno", "*_ko.tsv")) if not f.endswith("_raw_ko.tsv")]
    
    annotation_outdir = os.path.join(args.outdir, "gene_annotation_final")
    os.makedirs(annotation_outdir, exist_ok=True)
    
    genomes_list = []

    for kofam_f in ko_files:
        basename = os.path.basename(kofam_f)
        genome = basename.replace("_ko.tsv", "")
        genomes_list.append(genome)
        
        # 1. Load KO file
        try:
            ko_df = pd.read_csv(kofam_f, sep="\t", header=0)
            if "ko" not in ko_df.columns and "KO" in ko_df.columns:
                ko_df.rename(columns={"KO": "ko"}, inplace=True)
            elif "ko" not in ko_df.columns and len(ko_df.columns) >= 3:
                ko_df.rename(columns={ko_df.columns[2]: "ko"}, inplace=True)
                ko_df.rename(columns={ko_df.columns[1]: "gene_id"}, inplace=True)
        except Exception:
            ko_df = pd.DataFrame(columns=["genome", "gene_id", "ko"])
            
        if "ko" in ko_df.columns:
            ko_df.rename(columns={"ko": "KO"}, inplace=True)
        else:
            ko_df["KO"] = np.nan
            
        if "gene_id" not in ko_df.columns:
            ko_df["gene_id"] = np.nan
            
        if "genome" not in ko_df.columns:
            ko_df["genome"] = genome

        # Subset to basic 3 cols for safety
        available_cols = [c for c in ["genome", "gene_id", "KO"] if c in ko_df.columns]
        ko_df = ko_df[available_cols]

        # 2. Load Gene Classification (HMM) file
        gene_f = os.path.join(args.outdir, "cust_db_anno", f"{genome}_methanohunt_gene_classification.tsv")
        try:
            hmm_df = pd.read_csv(gene_f, sep="\t", header=0)
        except Exception:
            hmm_df = pd.DataFrame(columns=["gene_id", "classification"])
            
        if "gene_id" not in hmm_df.columns:
            hmm_df["gene_id"] = np.nan
        if "classification" not in hmm_df.columns:
            hmm_df["classification"] = np.nan

        hmm_df.rename(columns={"classification": "methanohunt-classification"}, inplace=True)
        
        # Subset to basic 2 cols
        available_hmm_cols = [c for c in ["gene_id", "methanohunt-classification"] if c in hmm_df.columns]
        hmm_df = hmm_df[available_hmm_cols]

        # 3. Merge
        merged_df = pd.merge(ko_df, hmm_df, on="gene_id", how="outer")
        merged_df["genome"] = genome
        
        # Ensure all 4 columns exist
        for c in ["genome", "gene_id", "KO", "methanohunt-classification"]:
            if c not in merged_df.columns:
                merged_df[c] = np.nan
                
        merged_df = merged_df[["genome", "gene_id", "KO", "methanohunt-classification"]]
        
        # Write annotation tsv
        out_anno = os.path.join(annotation_outdir, f"{genome}_methanohunt_annotation.tsv")
        merged_df.to_csv(out_anno, sep="\t", index=False)
        
        # 4. Compute presence/absence counts
        genome_pa_list = []
        genome_count_list = []
        
        for idx, row in db_df.iterrows():
            ko_db = row.get("KO", np.nan)
            class_db = row.get("Methanohunt_classification", np.nan)
            
            # Mask
            if pd.notna(class_db) and str(class_db).strip() != "":
                # match both KO and methanohunt-classification
                matched = merged_df[
                    (merged_df["KO"] == ko_db) &
                    (merged_df["methanohunt-classification"] == class_db)
                ]
            else:
                # match only KO
                matched = merged_df[merged_df["KO"] == ko_db]
                
            gene_ids = [str(gid) for gid in matched["gene_id"].dropna().unique() if str(gid) != "nan" and str(gid) != ""]
            
            if len(gene_ids) > 0:
                genome_pa_list.append(";".join(gene_ids))
                genome_count_list.append(len(gene_ids))
            else:
                genome_pa_list.append("")
                genome_count_list.append(0)
                
        pa_df[genome] = genome_pa_list
        count_df[genome] = genome_count_list

    # Save PA and Count metrics
    pa_df.to_csv(args.pa, sep="\t", index=False)
    count_df.to_csv(args.count, sep="\t", index=False)
    
    # Generate transposed matrix
    tr_df = count_df.copy()
    
    if "Methanohunt_classification" in tr_df.columns and "Subunit" in tr_df.columns:
        mask = tr_df["Methanohunt_classification"].isna() | (tr_df["Methanohunt_classification"] == "")
        tr_df.loc[mask, "Methanohunt_classification"] = tr_df.loc[mask, "Subunit"]
        
        # Subset
        cols_to_keep = ["Methanohunt_classification"] + genomes_list
        tr_df = tr_df[cols_to_keep]
        
        # Transpose
        tr_df = tr_df.set_index("Methanohunt_classification").T
        tr_df.index.name = "Genome"
        tr_df.reset_index(inplace=True)
    else:
        # Fallback if columns are missing
        tr_df = tr_df.T
        tr_df.index.name = "Genome"
        tr_df.reset_index(inplace=True)
        
    tr_df.to_csv(args.count_tr, sep="\t", index=False)
    
    print(f"Generated annotation pieces in {annotation_outdir}")
    print(f"Saved matrices to {args.pa}, {args.count}, and {args.count_tr}")

if __name__ == "__main__":
    main()
