# Pplacer rules

rule pplacer:
    input:
        sto = f"{config['output_dir']}/align/{wildcards.marker}.sto",
        refpkg = lambda wildcards: f"{config['database']}/{wildcards.marker}.refpkg"
    output:
        jplace = f"{config['output_dir']}/pplacer/{wildcards.marker}.jplace"
    threads: 4
    shell:
        # pplacer requires the refpkg
        "pplacer -c {input.refpkg} --out-dir {config['output_dir']}/pplacer {input.sto}"
