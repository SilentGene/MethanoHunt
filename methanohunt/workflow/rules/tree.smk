# Rules for phylogenetic tree construction using FastTree

rule merge_sequences:
    input:
        ref = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_reference.fasta",
        qry = f"{config['output_dir']}/cust_db_anno/placement/{{genome}}_{{marker}}_query.fasta"
    output:
        combined = f"{config['output_dir']}/tree/{{genome}}_{{marker}}_combined.fasta"
    shell:
        """
        cat {input.ref} {input.qry} > {output.combined}
        """

rule fasttree_build:
    input:
        combined = f"{config['output_dir']}/tree/{{genome}}_{{marker}}_combined.fasta"
    output:
        tree = f"{config['output_dir']}/tree/{{genome}}_{{marker}}.fasttree.nwk"
    threads: 4
    shell:
        """
        FastTree -lg < {input.combined} > {output.tree}
        """
