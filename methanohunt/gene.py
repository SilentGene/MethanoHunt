import os
import sys
import subprocess
import yaml
import tempfile
import click

import glob

def run_gene_pipeline(prot, nucl, reads_1, reads_2, mapper, database, output, marker, tree, threads=4):
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
        "threads": threads,
        "mapper": mapper
    }

    if nucl:
        config["nucl"] = os.path.abspath(nucl)
    
    # Process reads if provided
    samples = {}
    if reads_1 and reads_2:
        # Helper to process read inputs (string or list -> list of files)
        def process_read_arg(arg):
            files = []
            
            # Normalize to list
            if isinstance(arg, str):
                parts = arg.split()
            elif isinstance(arg, list):
                # If list, flatten it (in case some elements are space-separated strings?)
                # Usually list of filenames unless user passed quoted string in list
                parts = []
                for item in arg:
                    parts.extend(item.split())
            else:
                return []

            for part in parts:
                # Expand glob (handle both glob pattern and literal filename)
                expanded = glob.glob(part)
                if expanded:
                    files.extend(sorted(expanded))
                else:
                    files.append(part)
            return files

        r1_files = process_read_arg(reads_1)
        r2_files = process_read_arg(reads_2)

        if len(r1_files) != len(r2_files):
             click.echo(f"Error: Number of forward reads ({len(r1_files)}) does not match reverse reads ({len(r2_files)})", err=True)
             sys.exit(1)
        
        if not r1_files:
             click.echo(f"Error: No read files found matching input.", err=True)
             sys.exit(1)

        # Pair them up and infer sample names
        # Logic: -1后的第一输入和-2后的第一个输入是一对
        for r1, r2 in zip(r1_files, r2_files):
            # Infer sample name from r1
            basename = os.path.basename(r1)
            # Try removing suffixes
            suffixes = [".1.fq.gz", "_1.fq.gz", ".1.fastq.gz", "_1.fastq.gz", 
                        ".1.fq", "_1.fq", ".1.fastq", "_1.fastq"]
            sample_name = basename
            for suf in suffixes:
                if basename.endswith(suf):
                    sample_name = basename[:-len(suf)]
                    break
            
            samples[sample_name] = {
                "r1": os.path.abspath(r1),
                "r2": os.path.abspath(r2)
            }
        
    config["samples"] = samples

    print(f"Starting MethanoHunt Gene Pipeline...")
    print(f"  Snakemake: {snakefile_path}")
    print(f"  Output: {output}")
    print(f"  Database: {database}")
    print(f"  Threads: {threads}")
    if samples:
        print(f"  Samples: {len(samples)} pairs detected.")

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
        "--cores", str(threads),
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
