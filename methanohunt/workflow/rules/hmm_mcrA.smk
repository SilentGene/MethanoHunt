# Rules for McrA identification
# Search McrA_N and McrA_C independently and merge

rule hmmsearch_mcrA_N:
    input:
        faa = config["protein"],
        hmm = f"{config['database']}/McrA/McrA_N-PF02745.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/McrA_N.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule hmmsearch_mcrA_C:
    input:
        faa = config["protein"],
        hmm = f"{config['database']}/McrA/McrA_C-PF02249.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/McrA_C.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule merge_mcrA_hits:
    input:
        n = f"{config['output_dir']}/hmm/McrA_N.tbl",
        c = f"{config['output_dir']}/hmm/McrA_C.tbl"
    output:
        tbl = f"{config['output_dir']}/hmm/McrA.tbl"
    shell:
        """
        cat {input.n} {input.c} > {output.tbl}
        """

rule extract_hits_mcrA:
    input:
        faa = config["protein"],
        tbl = f"{config['output_dir']}/hmm/McrA.tbl"
    output:
        faa = f"{config['output_dir']}/hits/McrA_hit.faa"
    script:
        "../scripts/extract_hits.py"
