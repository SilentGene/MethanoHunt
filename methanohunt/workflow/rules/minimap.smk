rule minimap_index:
    input:
        ref = f"{config['output_dir']}/mapping_reference.fna"
    output:
        idx = f"{config['output_dir']}/mapping_reference.fna.mmi"
    threads: 1
    shell:
        "minimap2 -d {output.idx} {input.ref}"

rule minimap_map:
    input:
        idx = f"{config['output_dir']}/mapping_reference.fna.mmi",
        r1 = lambda wildcards: config["samples"][wildcards.sample]["r1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["r2"]
    output:
        bam = f"{config['output_dir']}/bam/{{sample}}.bam"
    threads: 4
    shell:
        """
        minimap2 -ax sr -t {threads} {input.idx} {input.r1} {input.r2} | \
        samtools view -bS - | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """
