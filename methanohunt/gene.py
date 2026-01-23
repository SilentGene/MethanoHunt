import os
import sys
import subprocess
import yaml
import click

import glob
import re

def strip_non_alnum_ends(s: str) -> str:
    return re.sub(r'^[^A-Za-z0-9]+|[^A-Za-z0-9]+$', '', s)

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
    if reads_1:
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

        if not r1_files:
             click.echo(f"Error: No forward read files found matching input.", err=True)
             sys.exit(1)
        
        # Enforce suffix rule
        valid_suffixes = ["1.fq.gz", "1.fastq.gz", "1.fq", "1.fastq"]
        detected_suffix = None
        
        for r1 in r1_files:
            basename = os.path.basename(r1)
            found = False
            for suf in valid_suffixes:
                if basename.endswith(suf):
                    if detected_suffix is None:
                        detected_suffix = suf
                    elif detected_suffix != suf:
                        click.echo(f"Error: All forward reads must have the same suffix. Found both '{detected_suffix}' and '{suf}'.", err=True)
                        sys.exit(1)
                    found = True
                    break
            if not found:
                click.echo(f"Error: Forward read '{r1}' does not end with any of {valid_suffixes}. Please rename your files to end with one of these suffixes.", err=True)
                sys.exit(1)

        # Pair them up and infer sample names
        for r1 in r1_files:
            basename = os.path.basename(r1)
            sample_prefix = basename[:-len(detected_suffix)]
            sample_name = strip_non_alnum_ends(sample_prefix)
            
            # Derive r2
            # e.g. sample1.1.fq.gz -> sample_prefix="sample1.", suffix="1.fq.gz"
            # r2 = "sample1." + "2" + ".fq.gz" = "sample1.2.fq.gz"
            r2_basename = sample_prefix + "2" + detected_suffix[1:]
            r2 = os.path.join(os.path.dirname(r1), r2_basename)
            
            if not os.path.exists(r2):
                click.echo(f"Error: Reverse read file '{r2}' not found (derived from '{r1}').", err=True)
                sys.exit(1)
            
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
        for name, data in samples.items():
            print(f"    - {name}:")
            print(f"      R1: {data['r1']}")
            print(f"      R2: {data['r2']}")

    # 5. Run Snakemake using subprocess (more robust across Snakemake versions)
    # Create config file in output directory
    config_path = os.path.join(output, "methanohunt_config.yaml")
    with open(config_path, 'w') as f:
        yaml.dump(config, f)
    
    print(f"Configuration written to: {config_path}")

    cmd = [
        "snakemake",
        "-s", snakefile_path,
        "--configfile", config_path,
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
        pass
