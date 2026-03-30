import os

outdir = config["output_dir"]
GENOMES = list(config["samples"].keys())

rule all:
    input:
        tsv = f"{outdir}/methanohunt_genome_classification.tsv"

if config["is_dna"]:
    rule prodigal:
        input:
            dna = lambda wildcards: config["samples"][wildcards.genome]
        output:
            faa = f"{outdir}/faa/{{genome}}.faa"
        shell:
            "prodigal -i {input.dna} -a {output.faa} -m -p meta -q > /dev/null 2>&1"
else:
    rule link_faa:
        input:
            faa = lambda wildcards: config["samples"][wildcards.genome]
        output:
            faa = f"{outdir}/faa/{{genome}}.faa"
        shell:
            "ln -s {input.faa} {output.faa}"

rule kofamscan_annotation:
    input:
        faa = f"{outdir}/faa/{{genome}}.faa",
        ko_list = config["ko_list"]
    output:
        raw = f"{outdir}/kofam/{{genome}}_raw_ko.tsv"
    params:
        profiles = config["profiles_dir"]
    threads: config["threads"]
    log:
        f"{outdir}/logs/kofam_{{genome}}.log"
    shell:
        "exec_annotation -p {params.profiles} -k {input.ko_list} --cpu {threads} -f detail-tsv --e-value 1e-5 --tmp-dir {outdir}/kofam/tmp_{wildcards.genome} -o {output.raw} {input.faa} > {log} 2>&1 && rm -rf {outdir}/kofam/tmp_{wildcards.genome}"

rule kofamscan_filter:
    input:
        raw = f"{outdir}/kofam/{{genome}}_raw_ko.tsv"
    output:
        filtered = f"{outdir}/kofam/{{genome}}_ko.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts", "kofamscan_filter.py")
    log:
        f"{outdir}/logs/kofam_filter_{{genome}}.log"
    shell:
        "python {params.script} -i {input.raw} -o {output.filtered} -s {wildcards.genome} -E 1e-5 >> {log} 2>&1"

MARKERS = ["McrA", "PmoA", "MmoX"]

include: "rules/hmm_mcrA.smk"
include: "rules/hmm_pmoA.smk"
include: "rules/hmm_mmoX.smk"
include: "rules/placement.smk"
include: "rules/classify.smk"
include: "rules/extract.smk"
include: "rules/tree.smk"

rule aggregate_genome_classification:
    input:
        ko_files = expand(f"{outdir}/kofam/{{genome}}_ko.tsv", genome=GENOMES),
        gene_files = expand(f"{outdir}/classification/{{genome}}_methanohunt_gene_classification.tsv", genome=GENOMES)
    output:
        summary = f"{outdir}/methanohunt_genome_classification.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts", "aggregate_genome_features.py"),
        db = config["database"],
        strict = "--strict" if config["strict"] else "",
        # Using string representation of python arg
        taxonomy = f"--taxonomy {config['taxonomy_file']}" if config.get("taxonomy_file") else "",
        tax_col = f"--col {config['taxonomy_col']}" if config.get("taxonomy_col") else "",
        outdir = outdir
    log:
        f"{outdir}/logs/aggregate_genome_features.log"
    shell:
        """
        python {params.script} \
            --outdir {params.outdir} \
            --db {params.db} \
            {params.strict} \
            {params.taxonomy} \
            {params.tax_col} \
            --output {output.summary} \
            > {log} 2>&1
        """
