# Rules for extracting classified sequences

rule extract_sequences:
    input:
        classification = f"{config['output_dir']}/classification/{{genome}}_methanohunt_gene_classification.tsv",
        protein = f"{config['output_dir']}/faa/{{genome}}.faa"
    output:
        mcrA_methanogen = f"{config['output_dir']}/classified_sequences/{{genome}}_McrA_methanogen.faa",
        mcrA_anme = f"{config['output_dir']}/classified_sequences/{{genome}}_McrA_ANME.faa",
        pmoA = f"{config['output_dir']}/classified_sequences/{{genome}}_PmoA.faa",
        mmoX = f"{config['output_dir']}/classified_sequences/{{genome}}_MmoX.faa"
    script:
        "../scripts/extract_sequences.py"
