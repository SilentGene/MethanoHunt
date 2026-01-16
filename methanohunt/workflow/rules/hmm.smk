# HMM Search, Extraction, and Alignment
# Markers: mcrA, pmoA, mmoA (defined by user)
MARKERS = ["mcrA", "pmoA", "mmoA"]

rule hmmsearch:
    input:
        faa = config["protein"],
        hmm = lambda wildcards: f"{config['database']}/{wildcards.marker}.hmm"
    output:
        tbl = f"{config['output_dir']}/hmm/{wildcards.marker}.tbl",
        # We don't strictly need the domtbl for simple detection but tbl contains per-sequence scores.
    threads: 4
    conda: "../envs/hmmer.yaml" # Optional, if using conda integration
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

rule align_hits:
    input:
        faa = f"{config['output_dir']}/hits/{wildcards.marker}_hit.faa",
        hmm = lambda wildcards: f"{config['database']}/{wildcards.marker}.hmm"
    output:
        sto = f"{config['output_dir']}/align/{wildcards.marker}.sto"
    threads: 2
    shell:
        # hmmalign aligns sequences to the profile
        "hmmalign --trim --outformat Stockholm -o {output.sto} {input.hmm} {input.faa}"
