# Rules for MmoX identification

rule hmmsearch_mmoX:
    input:
        faa = f"{config['output_dir']}/faa/{{genome}}.faa",
        hmm = f"{config['database']}/MmoX/MmoX-PF02332.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/{{genome}}_MmoX.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule extract_hits_mmoX:
    input:
        faa = f"{config['output_dir']}/faa/{{genome}}.faa",
        tbl = f"{config['output_dir']}/hmm/{{genome}}_MmoX.tbl"
    output:
        faa = f"{config['output_dir']}/hits/{{genome}}_MmoX_hit.faa"
    script:
        "../scripts/extract_hits.py"
