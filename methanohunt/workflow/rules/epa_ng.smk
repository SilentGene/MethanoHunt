# EPA-ng and PaPaRa rules for phylogenetic placement

rule papara_align:
    input:
        query = f"{config['output_dir']}/hits/{{marker}}_hit.faa",
        tree = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}_refs_methanohunt.faa_LG+F+I+R4.treefile",
        phylip = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}_refs_methanohunt.faa.phylip"
    output:
        aln = f"{config['output_dir']}/epa/papara_alignment.{{marker}}"
    params:
        prefix = "{marker}"
    threads: 1
    shell:
        """
        papara_static_x86_64 -t {input.tree} -s {input.phylip} -q {input.query} -r -a -n {params.prefix}
        mv papara_alignment.{params.prefix} {output.aln}
        """

rule epa_split:
    input:
        msa = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}_refs_methanohunt.faa.linsi",
        papara_aln = f"{config['output_dir']}/epa/papara_alignment.{{marker}}"
    output:
        ref = f"{config['output_dir']}/epa/{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/epa/{{marker}}_query.fasta"
    shadow: "minimal"
    shell:
        """
        epa-ng --split {input.msa} {input.papara_aln}
        mv reference.fasta {output.ref}
        mv query.fasta {output.qry}
        """

rule epa_place:
    input:
        ref = f"{config['output_dir']}/epa/{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/epa/{{marker}}_query.fasta",
        tree = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}_refs_methanohunt.faa_LG+F+I+R4.treefile",
        model = lambda wildcards: f"methanohunt/database/{{marker.capitalize()}}/{{marker}}.raxml.bestModel"
    output:
        jplace = f"{config['output_dir']}/epa/{{marker}}/epa_result.jplace"
    params:
        outdir = f"{config['output_dir']}/epa/{{marker}}"
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
