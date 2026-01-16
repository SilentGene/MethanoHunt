# Visualization Rule

def get_vis_input(wildcards):
    if config.get("nucl") and config.get("reads_1"):
        return f"{config['output_dir']}/abundance.tsv"
    else:
        # If no abundance, we visualise counts from classification
        return f"{config['output_dir']}/classification.tsv"

rule visualize:
    input:
        data = get_vis_input
    output:
        html = f"{config['output_dir']}/MethanoHunt_report.html"
    script:
        "../scripts/visualize.py"
