# Classification Rules using gappa

rule gappa_assign:
    input:
        jplace = f"{config['output_dir']}/placement/{{marker}}/epa_result.jplace",
        taxon = f"{config['database']}/{{marker}}/{{marker}}_refs_classification.tsv"
    output:
        results = f"{config['output_dir']}/classification/{{marker}}-per_query.tsv"
    params:
        outdir = f"{config['output_dir']}/classification",
        prefix = "{marker}"
    shell:
        """
        gappa examine assign \
          --jplace-path {input.jplace} \
          --taxon-file {input.taxon} \
          --per-query-results \
          --best-hit \
          --out-dir {params.outdir} \
          --file-prefix {params.prefix}"-"
        """

rule classify_genes:
    input:
        results = expand(f"{config['output_dir']}/classification/{{marker}}-per_query.tsv", marker=MARKERS)
    output:
        tsv = f"{config['output_dir']}/methanohunt_gene_classification.tsv"
    run:
        import pandas as pd
        all_dfs = []
        for f in input.results:
            df = pd.read_csv(f, sep='\t')
            # Extract marker from filename
            marker = f.split('/')[-1].replace('-per_query.tsv', '')
            df['marker'] = marker
            all_dfs.append(df)
        
        if not all_dfs:
            # Create empty df with expected columns if no inputs
            final_df = pd.DataFrame(columns=["name", "classification", "marker_hmm_hit"])
        else:
            final_df = pd.concat(all_dfs)
            
            # Rename columns
            if "taxopath" in final_df.columns:
                final_df.rename(columns={"taxopath": "classification"}, inplace=True)
            
            # Reorder columns: name, ..., marker_hmm_hit, classification
            # First, get all columns except marker and classification
            cols = [c for c in final_df.columns if c not in ["marker", "classification"]]
            
            # Construct new order
            # Note: 'marker' is what we added, rename it to 'marker_hmm_hit'
            final_df["marker_hmm_hit"] = final_df["marker"]
            
            new_order = cols + ["marker_hmm_hit", "classification"]
            
            # Filter valid columns in case some are missing
            new_order = [c for c in new_order if c in final_df.columns]
            
            final_df = final_df[new_order]

        final_df.to_csv(output.tsv, sep='\t', index=False)

