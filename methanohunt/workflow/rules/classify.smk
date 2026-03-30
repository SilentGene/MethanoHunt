# Classification Rules using gappa

rule gappa_assign:
    input:
        jplace = f"{config['output_dir']}/placement/{{genome}}/{{marker}}/epa_result.jplace",
        taxon = f"{config['database']}/{{marker}}/{{marker}}_refs_classification.tsv"
    output:
        results = f"{config['output_dir']}/classification/{{genome}}_{{marker}}-per_query.tsv"
    params:
        outdir = f"{config['output_dir']}/classification",
        prefix = "{genome}_{marker}"
    shell:
        """
        if [ -s {input.jplace} ]; then
            gappa examine assign \
              --jplace-path {input.jplace} \
              --taxon-file {input.taxon} \
              --per-query-results \
              --best-hit \
              --out-dir {params.outdir} \
              --file-prefix {params.prefix}"-"
        else
            echo "Empty jplace file, skipping gappa assign"
            touch {output.results}
        fi
        """

rule classify_genes:
    input:
        results = lambda wildcards: expand(f"{config['output_dir']}/classification/{{genome}}_{{marker}}-per_query.tsv", genome=wildcards.genome, marker=MARKERS)
    output:
        tsv = f"{config['output_dir']}/classification/{{genome}}_methanohunt_gene_classification.tsv"
    run:
        import pandas as pd
        all_dfs = []
        for f in input.results:
            try:
                df = pd.read_csv(f, sep='\t')
            except pd.errors.EmptyDataError:
                continue
            
            if df.empty:
                continue
                
            # Extract marker from filename (which is now {genome}_{marker}-per_query.tsv)
            filename = f.split('/')[-1].replace('-per_query.tsv', '')
            marker = 'Unknown'
            for m in MARKERS:
                if filename.endswith(m):
                    marker = m
                    break
            df['marker'] = marker
            all_dfs.append(df)
        
        if not all_dfs:
            # Create empty df with expected columns if no inputs
            final_df = pd.DataFrame(columns=["gene_id", "classification", "marker_hmm_hit"])
        else:
            final_df = pd.concat(all_dfs)
            
            # Rename columns: output of gappa usually has 'taxopath' if using --per-query-results --best fit?
            # actually gappa output likely has 'taxopath' or 'classification'?
            # Rule says: --taxon-file calls it "classification" often? Wait, gappa output header depends on version.
            # Assuming 'taxopath' is what we had before.
            
            if "taxopath" in final_df.columns:
                final_df.rename(columns={"taxopath": "subtype"}, inplace=True)
            elif "classification" in final_df.columns:
                 # If previously renamed or output differently
                 final_df.rename(columns={"classification": "subtype"}, inplace=True)
                 
            if "name" in final_df.columns:
                final_df.rename(columns={"name": "gene_id"}, inplace=True)
            
            # Derive new 'classification' column
            def get_classification(val):
                val = str(val)
                if val.startswith("Methanogen"): return "McrA (methanogen)"
                if val.startswith("ANME"): return "McrA (ANME)"
                if val.startswith("ANKA"): return "McrA (ANKA)"
                if val.startswith("PmoA"): return "PmoA"
                if val.startswith("AmoA"): return "AmoA"
                if val.startswith("HMO"): return "HMO"
                if val.startswith("MmoX"): return "MmoX"
                if val.startswith("PrmA"): return "PrmA"
                if val.startswith("MimA"): return "MimA"
                return "Unknown" # Or raise error as per request

            # User says: "如果是其他值，则报错说找不到正确的分类信息"
            # We can implement this inside the apply or iterate
            
            new_classifications = []
            for idx, row in final_df.iterrows():
                sub = str(row.get("subtype", ""))
                cls = get_classification(sub)
                if cls == "Unknown":
                     # For now, print warning or raise Error?
                     # Snakemake might catch error.
                     raise ValueError(f"Unknown classification subtype: {sub}")
                new_classifications.append(cls)
            
            final_df["classification"] = new_classifications

            # Note: 'marker' is what we added, rename it to 'marker_hmm_hit'
            # Map marker values to descriptive HMM names
            marker_desc_map = {
                "McrA": "McrA (PF02249.hmm/PF02745.hmm)",
                "PmoA": "PmoA/AmoA (PF02461.hmm)",
                "MmoX": "MmoX (PF02332.hmm)"
            }
            final_df["marker_hmm_hit"] = final_df["marker"].map(lambda x: marker_desc_map.get(str(x), str(x)))

            # Reorder columns: 
            # Request: "classification" in 2nd to last position.
            # Request Step 173: "marker_hmm_hit" in 2nd to last.
            # Current request: "classification" in 2nd to last.
            # Let's assume order: ..., marker_hmm_hit, classification, subtype? 
            # Or ..., classification, subtype?
            # Step 276: 
            # 1. Modify existing 'classification' -> 'subtype'.
            # 2. Add NEW 'classification' col at 2nd to last position?
            # Wait, if I put it 2nd to last, what is last? 
            # Maybe user means: [..., classification, last_col]
            # "modify 'classification' (now subtype) ... add new 'classification' at 2nd to last..."
            # Usually users mean relative to ends.
            # Let's target: name, [metadata], marker_hmm_hit, classification, subtype
            # Or name, [metadata], classification, subtype?
            # Let's stick to a logical order or exactly as implied?
            # "increase a new column 'classification' at 2nd to last position".
            # If we preserve old order logic: name, ..., marker_hmm_hit, classification (old).
            # Now: name, ..., marker_hmm_hit, classification (new), subtype (old classification output).
            
            # First, get all columns except our manipulated ones
            cols = [c for c in final_df.columns if c not in ["marker", "subtype", "classification", "marker_hmm_hit", "taxopath"]]
            
            # Construct new order
            # Let's place: name, ..., marker_hmm_hit, classification, subtype
            # Be careful about column existence
            
            new_order = cols + ["marker_hmm_hit", "classification", "subtype"]
            # remove duplicates if any
            
            # Filter valid columns
            new_order = [c for c in new_order if c in final_df.columns]
            
            final_df = final_df[new_order]

        final_df.to_csv(output.tsv, sep='\t', index=False)

