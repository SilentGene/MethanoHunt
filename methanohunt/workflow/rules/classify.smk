# Classification Rule

rule classify_genes:
    input:
        jplaces = expand(f"{config['output_dir']}/pplacer/{{marker}}.jplace", marker=["mcrA", "pmoA", "mmoA"]),
        db = f"{config['database']}/methanohunt_db.tsv" 
        # Note: The db file here is the one used for taxonomy but user implied similar logic.
        # However, for gene classification, we probably need mapped classifications for the tree leaves or branches.
        # The user said: "classify protein... similar outcome as taxonomy module".
        # User provided 'reference tree'. Usually pplacer output needs 'guppy' or similar to classify?
        # Or we parse the jplace.
        # User said: "classify... based on three reference tree... then use this info to divide into Methanogen, ANME..."
        # I will assume the python script handles the logic of parsing jplace and assigning groups.
    output:
        tsv = f"{config['output_dir']}/classification.tsv"
    script:
        "../scripts/classify_genes.py"
