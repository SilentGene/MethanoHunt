# Rules for PmoA identification

rule hmmsearch_pmoA:
    input:
        faa = config["protein"],
        hmm = f"{config['database']}/PmoA/PmoAAmoA-PF02461.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/PmoA.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule extract_hits_pmoA:
    input:
        faa = config["protein"],
        tbl = f"{config['output_dir']}/hmm/PmoA.tbl"
    output:
        faa = f"{config['output_dir']}/hits/PmoA_hit.faa"
    script:
        "../scripts/extract_hits.py"
