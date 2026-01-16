import os
import sys
import snakemake
from snakemake import snakemake as run_snakemake
import click

def run_gene_pipeline(protein, nucl, reads_1, reads_2, database, output):
    """
    Wrapper to execute the Snakemake pipeline.
    """
    # 1. Resolve output directory
    if not output:
        base_name = os.path.splitext(os.path.basename(protein))[0]
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
    config = {
        "protein": os.path.abspath(protein),
        "database": os.path.abspath(database),
        "output_dir": os.path.abspath(output),
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

    # 5. Run Snakemake
    # We use the python API
    success = run_snakemake(
        snakefile=snakefile_path,
        config=config,
        cores=4, # Default cores, maybe configurable?
        printshellcmds=True,
        keepgoing=True,
        # workdir=output # Running in output dir makes relative paths simpler, but we used absolute paths above so it's fine.
        # Actually, running in CWD is usually safer for finding package files? 
        # But we want outputs in output_dir. 
        # Let's keep workdir as CWD and direct all rule outputs to config['output_dir'].
    )

    if not success:
        click.echo("Error: Snakemake pipeline failed.", err=True)
        sys.exit(1)
    else:
        click.echo("Pipeline completed successfully.")
