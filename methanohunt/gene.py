import os
import sys
import subprocess
import yaml
import tempfile
import click

def run_gene_pipeline(prot, nucl, reads_1, reads_2, database, output, marker, tree):
    """
    Wrapper to execute the Snakemake pipeline.
    """
    # 1. Resolve output directory
    if not output:
        base_name = os.path.splitext(os.path.basename(prot))[0]
        output = f"{base_name}_methanohunt"
    
    if not os.path.exists(output):
        os.makedirs(output)
        print(f"Created output directory: {output}")

    # 2. Resolve database path
    if not database:
        # Default to methanohunt/database
        # Assuming current structure: methanohunt/gene.py -> .. -> methanohunt/database
        current_dir = os.path.dirname(os.path.abspath(__file__))
        database = os.path.join(current_dir, "database")
    
    if not os.path.exists(database):
        click.echo(f"Error: Database directory not found at {database}", err=True)
        sys.exit(1)

    # 3. Locate Snakefile
    # methanohunt/gene.py -> methanohunt/workflow/Snakefile
    current_dir = os.path.dirname(os.path.abspath(__file__))
    snakefile_path = os.path.join(current_dir, "workflow", "Snakefile")

    if not os.path.exists(snakefile_path):
         click.echo(f"Error: Snakefile not found at {snakefile_path}", err=True)
         sys.exit(1)

    # 4. Construct configuration
    if marker:
        # User requested markers: "Method 1: preserve user case; Method 2: force capitalize"
        # User requested: "CLI command let user use --marker PmoA,McrA,MmoX... uppercase way"
        # It implies user should type it that way or we enforce expected internal names.
        # Since file system is case sensitive (or specific names used in rules), we should match what the rules expect.
        # The rules expect: McrA, PmoA, MmoX.
        # So we should probably capitalize them or trust user. 
        # Safest is to Map reasonable inputs to expected outputs, or just pass as is if user follows instructions.
        # User said "let user use... uppercase way", implies user input WILL BE uppercase.
        # But we can be helpful and title case it or verify.
        # Let's trust user but maybe title() if we want strictly standard behavior, but rules rely on exact match.
        markers_list = [m.strip() for m in marker.split(",") if m.strip()]
    else:
        # Defaults
        markers_list = ["McrA", "PmoA", "MmoX"]
    
    config = {
        "protein": os.path.abspath(prot),
        "database": os.path.abspath(database),
        "output_dir": os.path.abspath(output),
        "markers": markers_list,
        "run_tree": tree,
    }

    if nucl:
        config["nucl"] = os.path.abspath(nucl)
    if reads_1:
         config["reads_1"] = os.path.abspath(reads_1)
    if reads_2:
         config["reads_2"] = os.path.abspath(reads_2)


    print(f"Starting MethanoHunt Gene Pipeline...")
    print(f"  Snakemake: {snakefile_path}")
    print(f"  Output: {output}")
    print(f"  Database: {database}")

    # 5. Run Snakemake using subprocess (more robust across Snakemake versions)
    # Create temporary config file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as tmp_config:
        yaml.dump(config, tmp_config)
        tmp_config_path = tmp_config.name
    
    print(f"Temporary config written to: {tmp_config_path}")

    cmd = [
        "snakemake",
        "-s", snakefile_path,
        "--configfile", tmp_config_path,
        "--cores", "4",
        "--printshellcmds",
        "--keep-going"
    ]
    
    # Run in the output directory if desired, or simpler: run from here but output is absolute
    # Since config paths are absolute, running from CWD is fine.
    
    try:
        print("Running Snakemake command:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        click.echo("Pipeline completed successfully.")
    except subprocess.CalledProcessError:
        click.echo("Error: Snakemake pipeline failed.", err=True)
        sys.exit(1)
    except FileNotFoundError:
        click.echo("Error: 'snakemake' command not found. Please ensure Snakemake is installed and in your PATH.", err=True)
        sys.exit(1)
    finally:
        # Cleanup temp config
        if os.path.exists(tmp_config_path):
             os.remove(tmp_config_path)
