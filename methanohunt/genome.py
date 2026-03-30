import os
import sys
import subprocess
import yaml
import click
import glob
import shlex

def run_genome_pipeline(genome_dir, faa_dir, suffix, profiles_dir, ko_list, taxonomy_file, taxonomy_col, threads, strict, output, snake_args=None):
    """
    Wrapper to execute the Snakemake genome pipeline.
    """
    # 1. Resolve output directory
    if not os.path.exists(output):
        os.makedirs(output)
        click.echo(f"Created output directory: {output}")

    # 2. Locate Database Path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    database = os.path.join(current_dir, "database")
    if not os.path.exists(database):
        click.echo(f"Error: Database directory not found at {database}", err=True)
        sys.exit(1)

    # 3. Locate Snakefile
    snakefile_path = os.path.join(current_dir, "workflow", "genome.smk")
    
    # Optional taxonomy parsing arguments check
    if taxonomy_file and not taxonomy_col:
        click.echo("Error: --col must be provided when --taxonomy is specified.", err=True)
        sys.exit(1)

    if taxonomy_col and not taxonomy_file:
         click.echo("Error: --taxonomy must be provided when --col is specified.", err=True)
         sys.exit(1)

    # 4. Find valid suffixes and build sample map
    samples = {}
    is_dna = genome_dir is not None
    search_dir = genome_dir if is_dna else faa_dir
    
    if not os.path.isdir(search_dir):
         click.echo(f"Error: Directory '{search_dir}' not found.", err=True)
         sys.exit(1)

    search_pattern = os.path.join(search_dir, f"*.{suffix}")
    found_files = glob.glob(search_pattern)

    if not found_files:
         click.echo(f"Error: No files matching '*.{suffix}' found in '{search_dir}'.", err=True)
         sys.exit(1)

    for f in found_files:
        basename = os.path.basename(f)
        sample_name = basename[:-len(suffix)-1]
        samples[sample_name] = os.path.abspath(f)

    # 5. Construct configuration
    config = {
        "is_dna": is_dna,
        "input_dir": os.path.abspath(search_dir),
        "suffix": suffix,
        "samples": samples,
        "profiles_dir": os.path.abspath(profiles_dir),
        "ko_list": os.path.abspath(ko_list),
        "database": os.path.abspath(database),
        "output_dir": os.path.abspath(output),
        "threads": threads,
        "strict": strict,
        "taxonomy_file": os.path.abspath(taxonomy_file) if taxonomy_file else None,
        "taxonomy_col": taxonomy_col
    }

    click.echo(f"Starting MethanoHunt Genome Pipeline...")
    click.echo(f"  Snakemake: {snakefile_path}")
    click.echo(f"  Output: {output}")
    click.echo(f"  Mode: {'DNA -> protein mapping' if is_dna else 'Protein direct mapping'}")
    click.echo(f"  Samples: {len(samples)} found.")
    click.echo(f"  Strict Mode: {strict}")

    # 6. Run Snakemake using subprocess
    config_path = os.path.join(config["output_dir"], "methanohunt_genome_config.yaml")
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    click.echo(f"Configuration written to: {config_path}")

    cmd = [
        "snakemake",
        "-s", snakefile_path,
        "--configfile", config_path,
        "--cores", str(threads),
        "--printshellcmds",
        "--keep-going"
    ]
    
    if snake_args:
        cmd.extend(shlex.split(snake_args))
    
    try:
        click.echo("Running Snakemake command: " + " ".join(cmd))
        subprocess.run(cmd, check=True, cwd=config["output_dir"])
        click.echo("Genome Pipeline completed successfully.")
    except subprocess.CalledProcessError:
        click.echo("Error: Snakemake pipeline failed.", err=True)
        sys.exit(1)
    except FileNotFoundError:
        click.echo("Error: 'snakemake' command not found. Please ensure Snakemake is installed and in your PATH.", err=True)
        sys.exit(1)
