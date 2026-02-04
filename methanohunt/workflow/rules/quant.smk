# Quantification Rules

# 1. Extract DNA reference for mapping
rule extract_dna_ref:
    input:
        classification = f"{config['output_dir']}/methanohunt_gene_classification.tsv",
        nucl = config["nucl"] if config.get("nucl") else []
    output:
        ref = f"{config['output_dir']}/classified_sequences/all_classified_sequences.ffn"
    script:
        "../scripts/extract_dna.py"

# Only proceed if samples/nucl are present
if config.get("samples") and config.get("nucl"):
    # 2. Select Mapper
    if config.get("mapper") == "bwa":
        include: "bwa.smk"
    else:
        include: "minimap.smk"

    # 3. MicrobeCensus Rule
    rule microbe_census:
        input:
            r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
            r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
        output:
            out = f"{config['output_dir']}/microbecensus/{{sample}}_MicrobeCensus.out"
        threads: 4
        shell:
            """
            run_microbe_census.py -t {threads} {input.r1},{input.r2} {output.out}
            """

    # 4. Calculate RPKG per sample
    rule calculate_rpkg:
        input:
            bam = f"{config['output_dir']}/bam/{{sample}}.bam",
            mc_out = f"{config['output_dir']}/microbecensus/{{sample}}_MicrobeCensus.out"
        output:
            rpkg = f"{config['output_dir']}/quant/{{sample}}_rpkg.tsv"
        script:
            "../scripts/calculate_rpkg_single.py"

    # 5. Aggregate RPKG
    rule aggregate_rpkg:
        input:
            rpkgs = expand(f"{config['output_dir']}/quant/{{sample}}_rpkg.tsv", sample=config["samples"].keys()),
            classification = f"{config['output_dir']}/methanohunt_gene_classification.tsv"
        output:
            combined = f"{config['output_dir']}/RPKG/RPKG_combined.tsv",
            class_tsv = f"{config['output_dir']}/RPKG/RPKG_classification.tsv",
            subtype_tsv = f"{config['output_dir']}/RPKG/RPKG_subtype.tsv"
        script:
            "../scripts/aggregate_rpkg.py"
