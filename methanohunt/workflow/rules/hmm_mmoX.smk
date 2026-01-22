# Rules for MmoX identification

rule hmmsearch_mmoX:
    input:
        faa = config["protein"],
        hmm = f"{config['database']}/MmoX/MmoX-PF02332.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/MmoX.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule extract_hits_mmoX:
    input:
        faa = config["protein"],
        tbl = f"{config['output_dir']}/hmm/MmoX.tbl"
    output:
        faa = f"{config['output_dir']}/hits/MmoX_hit.faa"
    script:
        "../scripts/extract_hits.py"
