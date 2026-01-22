# Abundance Rule

# Only define rule if inputs are present
if config.get("nucl") and config.get("reads_1"):
    rule calculate_abundance:
        input:
            classification = f"{config['output_dir']}/methanohunt_gene_classification.tsv",
            gene_nucl = config["nucl"],
            reads_1 = config["reads_1"],
            reads_2 = config.get("reads_2", "") # Handle single end if needed? Logic below assumes PE
        output:
            abundance = f"{config['output_dir']}/abundance.tsv",
            bam = f"{config['output_dir']}/mapping/reads.bam"
        threads: 8
        params:
            reads_2_arg = lambda wildcards, input: input.reads_2 if input.reads_2 else ""
        shell:
            """
            minimap2 -ax sr -t {threads} {input.gene_nucl} {input.reads_1} {params.reads_2_arg} | \
            samtools view -bS - | samtools sort -o {output.bam} -
            samtools index {output.bam}
            python methanohunt/workflow/scripts/calculate_tpm.py \
                --bam {output.bam} \
                --classification {input.classification} \
                --output {output.abundance}
            """
