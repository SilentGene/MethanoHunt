rule bwa_index:
    input:
        ref = f"{config['output_dir']}/classified_sequences/all_classified_sequences.ffn"
    output:
        idx = f"{config['output_dir']}/classified_sequences/all_classified_sequences.ffn.sa"
    threads: 1
    shell:
        "bwa index {input.ref}"

rule bwa_map:
    input:
        ref = f"{config['output_dir']}/classified_sequences/all_classified_sequences.ffn",
        idx = f"{config['output_dir']}/classified_sequences/all_classified_sequences.ffn.sa",
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        bam = f"{config['output_dir']}/bam/{{sample}}.bam"
    threads: 4
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools view -bS - | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """
