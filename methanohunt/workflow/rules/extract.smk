# Rules for extracting classified sequences

rule extract_sequences:
    input:
        classification = f"{config['output_dir']}/methanohunt_gene_classification.tsv",
        protein = config["protein"]
    output:
        mcrA_methanogen = f"{config['output_dir']}/classified_sequences/McrA_methanogen.faa",
        mcrA_anme = f"{config['output_dir']}/classified_sequences/McrA_ANME.faa",
        pmoA = f"{config['output_dir']}/classified_sequences/PmoA.faa",
        mmoX = f"{config['output_dir']}/classified_sequences/MmoX.faa"
    script:
        "../scripts/extract_sequences.py"
