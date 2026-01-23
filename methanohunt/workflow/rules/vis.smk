# Visualization Rule

rule visualize:
    input:
        class_tsv = f"{config['output_dir']}/RPKG/RPKG_classification.tsv",
        subtype_tsv = f"{config['output_dir']}/RPKG/RPKG_subtype.tsv"
    output:
        html = f"{config['output_dir']}/MethanoHunt_report.html"
    script:
        "../scripts/visualize.py"
