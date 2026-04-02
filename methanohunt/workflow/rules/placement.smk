# EPA-ng and PaPaRa rules for phylogenetic placement

# EPA-ng and PaPaRa rules for phylogenetic placement

rule papara_align:
    input:
        query = f"{config['output_dir']}/cust_db_anno/hits/{{genome}}_{{marker}}_hit.faa",
        tree = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.treefile",
        phylip = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.phylip"
    output:
        aln = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_papara_alignment.{{marker}}",
        log = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_papara_log.{{marker}}",
        quality = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_papara_quality.{{marker}}"
    params:
        prefix = "{marker}"
    threads: 1
    shadow: "minimal"
    shell:
        """
        if [ -s {input.query} ]; then
            papara -t {input.tree} -s {input.phylip} -q {input.query} -r -a -n {params.prefix}
            mv papara_alignment.{params.prefix} {output.aln}
            mv papara_log.{params.prefix} {output.log}
            mv papara_quality.{params.prefix} {output.quality}
        else
            echo "Empty query file, skipping papara"
            touch {output.aln} {output.log} {output.quality}
        fi
        """

rule epa_split:
    input:
        msa = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.linsi",
        papara_aln = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_papara_alignment.{{marker}}"
    output:
        ref = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_query.fasta"
    shadow: "minimal"
    shell:
        """
        if [ -s {input.papara_aln} ]; then
            epa-ng --split {input.msa} {input.papara_aln}
            mv reference.fasta {output.ref}
            mv query.fasta {output.qry}
        else
            echo "Empty papara alignment, skipping epa-ng split"
            touch {output.ref} {output.qry}
        fi
        """

rule epa_place:
    input:
        ref = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_query.fasta",
        tree = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}_refs_methanohunt.faa.treefile",
        model = lambda wildcards: f"{config['database']}/{wildcards.marker}/{wildcards.marker}.raxml.bestModel"
    output:
        jplace = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}/{{marker}}/epa_result.jplace"
    params:
        outdir = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}/{{marker}}"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        if [ -s {input.qry} ]; then
            epa-ng \
            --ref-msa {input.ref} \
            --tree    {input.tree} \
            --query   {input.qry} \
            --model   {input.model} \
            --outdir  {params.outdir} \
            --threads {threads}
        else
            echo "Empty query fasta, skipping epa-ng place"
            touch {output.jplace}
        fi
        """
