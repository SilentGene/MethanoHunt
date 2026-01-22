# Visualization Rule

rule visualize:
    input:
        data = f"{config['output_dir']}/abundance.tsv"
    output:
        html = f"{config['output_dir']}/MethanoHunt_report.html"
    script:
        "../scripts/visualize.py"
