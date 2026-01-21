# Classification Rules using gappa

rule gappa_assign:
    input:
        jplace = f"{config['output_dir']}/epa/{{marker}}/epa_result.jplace",
        taxon = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}_refs_classification.tsv"
    output:
        results = f"{config['output_dir']}/classification/{{marker}}.per_query.results"
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
          --file-prefix {params.prefix}
        """

rule classify_genes:
    input:
        results = expand(f"{config['output_dir']}/classification/{{marker}}.per_query.results", marker=MARKERS)
    output:
        tsv = f"{config['output_dir']}/classification.tsv"
    run:
        import pandas as pd
        all_dfs = []
        for f in input.results:
            df = pd.read_csv(f, sep='\t')
            # Extract marker from filename or logic
            marker = f.split('/')[-1].replace('.per_query.results', '')
            df['marker'] = marker
            all_dfs.append(df)
        
        final_df = pd.concat(all_dfs)
        final_df.to_csv(output.tsv, sep='\t', index=False)

