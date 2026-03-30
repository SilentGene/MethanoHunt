# Rules for McrA identification
# Search McrA_N and McrA_C independently and merge

rule hmmsearch_mcrA_N:
    input:
        faa = f"{config['output_dir']}/faa/{{genome}}.faa",
        hmm = f"{config['database']}/McrA/McrA_N-PF02745.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/{{genome}}_McrA_N.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule hmmsearch_mcrA_C:
    input:
        faa = f"{config['output_dir']}/faa/{{genome}}.faa",
        hmm = f"{config['database']}/McrA/McrA_C-PF02249.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/{{genome}}_McrA_C.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule merge_mcrA_hits:
    input:
        n = f"{config['output_dir']}/hmm/{{genome}}_McrA_N.tbl",
        c = f"{config['output_dir']}/hmm/{{genome}}_McrA_C.tbl"
    output:
        tbl = f"{config['output_dir']}/hmm/{{genome}}_McrA.tbl"
    shell:
        """
        cat {input.n} {input.c} > {output.tbl}
        """

rule extract_hits_mcrA:
    input:
        faa = f"{config['output_dir']}/faa/{{genome}}.faa",
        tbl = f"{config['output_dir']}/hmm/{{genome}}_McrA.tbl"
    output:
        faa = f"{config['output_dir']}/hits/{{genome}}_McrA_hit.faa"
    script:
        "../scripts/extract_hits.py"
