# HMM Search, Extraction, and Alignment
# Markers: mcrA, pmoA, mmoA (defined by user)
MARKERS = ["mcrA", "pmoA", "mmoA"]

ruleorder: hmmsearch_mcrA > hmmsearch

rule hmmsearch_mcrA:
    input:
        faa = config["protein"],
        hmm1 = f"{config['database']}/McrA/McrA_N-PF02745.hmm",
        hmm2 = f"{config['database']}/McrA/McrA_C-PF02249.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/mcrA.tbl"
    threads: 6
    shell:
        """
        hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl}.tmp1 {input.hmm1} {input.faa} > /dev/null
        hmmsearch --cpu {threads} --cut_tc --tblout {output.tbl}.tmp2 {input.hmm2} {input.faa} > /dev/null
        cat {output.tbl}.tmp1 {output.tbl}.tmp2 > {output.tbl}
        rm {output.tbl}.tmp1 {output.tbl}.tmp2
        """

rule hmmsearch:
    input:
        faa = config["protein"],
        hmm = lambda wildcards: f"{config['database']}/{wildcards.marker}.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/{wildcards.marker}.tbl",
        # We don't strictly need the domtbl for simple detection but tbl contains per-sequence scores.
    threads: 4
    shell:
        "hmmsearch --cpu {threads} --tblout {output.tbl} {input.hmm} {input.faa} > /dev/null"

rule extract_hits:
    input:
        faa = config["protein"],
        tbl = f"{config['output_dir']}/hmm/{wildcards.marker}.tbl"
    output:
        faa = f"{config['output_dir']}/hits/{wildcards.marker}_hit.faa"
    script:
        "../scripts/extract_hits.py"

