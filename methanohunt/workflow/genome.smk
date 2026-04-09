import os

outdir = config["output_dir"]
GENOMES = list(config["samples"].keys())

rule all:
    input:
        tsv = f"{outdir}/methanohunt_genome_classification.tsv",
        pa = f"{outdir}/methanohunt_gene_presence_absence.tsv",
        pa_count = f"{outdir}/methanohunt_gene_presence_absence_count.tsv",
        pa_count_tr = f"{outdir}/methanohunt_gene_presence_absence_count_tr.tsv"

if config["is_dna"]:
    rule prodigal:
        input:
            dna = lambda wildcards: config["samples"][wildcards.genome]
        output:
            faa = f"{outdir}/ORF/{{genome}}.faa",
            ffn = f"{outdir}/ORF/{{genome}}.ffn",
            gff = f"{outdir}/ORF/{{genome}}.gff"
        shell:
            "prodigal -i {input.dna} -a {output.faa} -d {output.ffn} -f gff -o {output.gff} -p meta -q > /dev/null 2>&1"
else:
    rule link_faa:
        input:
            faa = lambda wildcards: config["samples"][wildcards.genome]
        output:
            faa = f"{outdir}/ORF/{{genome}}.faa"
        shell:
            "ln -s {input.faa} {output.faa}"

rule kofamscan_annotation:
    input:
        faa = f"{outdir}/ORF/{{genome}}.faa",
        ko_list = config["ko_list"]
    output:
        raw = f"{outdir}/KofamScan_anno/kofam_raw/{{genome}}_raw_ko.tsv"
    params:
        profiles = config["profiles_dir"]
    threads: max(min(4, int(config["threads"])), min(8, int(config["threads"]) // max(1, len(GENOMES))))
    log:
        f"{outdir}/logs/kofam_{{genome}}.log"
    shell:
        "exec_annotation -p {params.profiles} -k {input.ko_list} --cpu {threads} -f detail-tsv --e-value 1e-5 --tmp-dir {outdir}/KofamScan_anno/tmp_{wildcards.genome} -o {output.raw} {input.faa} > {log} 2>&1 && rm -rf {outdir}/KofamScan_anno/tmp_{wildcards.genome}"

rule kofamscan_filter:
    input:
        raw = f"{outdir}/KofamScan_anno/kofam_raw/{{genome}}_raw_ko.tsv"
    output:
        filtered = f"{outdir}/KofamScan_anno/{{genome}}_ko.tsv"
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
        ko_files = expand(f"{outdir}/KofamScan_anno/{{genome}}_ko.tsv", genome=GENOMES),
        gene_files = expand(f"{outdir}/cust_db_anno/{{genome}}_methanohunt_gene_classification.tsv", genome=GENOMES)
    output:
        summary = f"{outdir}/methanohunt_genome_classification.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts", "aggregate_genome_features.py"),
        db = config["database"],
        strict = "--strict" if config["strict"] else "",
        # Using string representation of python arg
         taxonomy = f"--taxonomy {config['taxonomy_file']}" if config.get("taxonomy_file") else "",
        name_col = f"--name {config['name_col']}" if config.get("name_col") else "",
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
            {params.name_col} \
            {params.tax_col} \
            --output {output.summary} \
            > {log} 2>&1
        """

rule gene_presence_absence:
    input:
        ko_files = expand(f"{outdir}/KofamScan_anno/{{genome}}_ko.tsv", genome=GENOMES),
        gene_files = expand(f"{outdir}/cust_db_anno/{{genome}}_methanohunt_gene_classification.tsv", genome=GENOMES)
    output:
        pa = f"{outdir}/methanohunt_gene_presence_absence.tsv",
        pa_count = f"{outdir}/methanohunt_gene_presence_absence_count.tsv",
        pa_count_tr = f"{outdir}/methanohunt_gene_presence_absence_count_tr.tsv"
    params:
        script = os.path.join(workflow.basedir, "scripts", "gene_presence_absence.py"),
        db = config["database"],
        outdir = outdir
    log:
        f"{outdir}/logs/gene_presence_absence.log"
    shell:
        """
        python {params.script} \
            --outdir {params.outdir} \
            --db {params.db} \
            --pa {output.pa} \
            --count {output.pa_count} \
            --count_tr {output.pa_count_tr} \
            > {log} 2>&1
        """
