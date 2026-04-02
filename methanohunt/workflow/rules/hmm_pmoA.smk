# Rules for PmoA identification

rule hmmsearch_pmoA:
    input:
        faa = f"{config['output_dir']}/ORF/{{genome}}.faa",
        hmm = f"{config['database']}/PmoA/PmoAAmoA-PF02461.hmm"
    output:
        tbl = f"{config['output_dir']}/cust_db_anno/hmm/{{genome}}_PmoA.tbl"
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule extract_hits_pmoA:
    input:
        faa = f"{config['output_dir']}/ORF/{{genome}}.faa",
        tbl = f"{config['output_dir']}/cust_db_anno/hmm/{{genome}}_PmoA.tbl"
    output:
        faa = f"{config['output_dir']}/cust_db_anno/hits/{{genome}}_PmoA_hit.faa"
    script:
        "../scripts/extract_hits.py"
