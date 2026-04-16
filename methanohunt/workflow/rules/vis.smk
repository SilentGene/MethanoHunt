# Visualization Rule

rule visualize:
    input:
        class_tsv = f"{config['output_dir']}/RPKG/RPKG_marker_gene.tsv",
        subtype_tsv = f"{config['output_dir']}/RPKG/RPKG_marker_gene_subtype.tsv"
    output:
        html = f"{config['output_dir']}/MethanoHunt_report.html"
    script:
        "../scripts/visualize.py"
