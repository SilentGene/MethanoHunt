# EPA-ng and PaPaRa rules for phylogenetic placement

rule papara_align:
    input:
        query = f"{config['output_dir']}/hits/{{marker}}_hit.faa",
        tree = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.treefile",
        phylip = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.phylip"
    output:
        aln = f"{config['output_dir']}/placement/papara_alignment.{{marker}}"
    params:
        prefix = "{marker}",
        outdir = f"{config['output_dir']}/placement/"
    threads: 1
    shell:
        """
        papara -t {input.tree} -s {input.phylip} -q {input.query} -r -a -n {params.prefix}
        mv papara_alignment.{params.prefix} {output.aln}
        mv papara_log.{params.prefix} {params.outdir}
        mv papara_quality.{params.prefix} {params.outdir}
        """

rule epa_split:
    input:
        msa = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.linsi",
        papara_aln = f"{config['output_dir']}/placement/papara_alignment.{{marker}}"
    output:
        ref = f"{config['output_dir']}/placement/{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/placement/{{marker}}_query.fasta"
    shadow: "minimal"
    shell:
        """
        epa-ng --split {input.msa} {input.papara_aln}
        mv reference.fasta {output.ref}
        mv query.fasta {output.qry}
        """

rule epa_place:
    input:
        ref = f"{config['output_dir']}/placement/{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/placement/{{marker}}_query.fasta",
        tree = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.treefile",
        model = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}.raxml.bestModel"
    output:
        jplace = f"{config['output_dir']}/placement/{{marker}}/epa_result.jplace"
    params:
        outdir = f"{config['output_dir']}/placement/{{marker}}"
    threads: config.get("threads", 8)
    shell:
        """
        mkdir -p {params.outdir}
        epa-ng \
        --ref-msa {input.ref} \
        --tree    {input.tree} \
        --query   {input.qry} \
        --model   {input.model} \
        --outdir  {params.outdir} \
        --threads {threads}
        """
