import click
import os
import sys
import shutil
import tarfile
import urllib.request
from .taxonomy import run_taxonomy

# from .gene import run_gene_pipeline  # To be implemented

@click.group()
def cli():
    """MethanoHunt: A pipeline for profiling methane cyclers."""
    pass

@cli.command()
@click.option("-i", "--input", "input_files", multiple=True, required=False, help="Input tax.tsv files (supports glob patterns)")
@click.option("-db", "--database", default=None, help="Path to MethanoHunt database file")
@click.option("-o", "--output", required=True, help="Output file prefix (will generate .tsv and .html)")
@click.argument("extra_files", nargs=-1)
def taxonomy(input_files, database, output, extra_files):
    """Run taxonomy-based profiling.
    
    Accepts input files via -i/--input or as positional arguments (e.g. *.tsv).
    """
    try:
        # Combine inputs from flag and positional args
        all_inputs = list(input_files) + list(extra_files)
        
        if not all_inputs:
             raise click.UsageError("No input files provided. Use -i or positional arguments.")

        # run_taxonomy expects list of patterns or files.
        run_taxonomy(all_inputs, database, output)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--prot", required=True, help="Input protein FASTA file (.faa)")
@click.option("--nucl", default=None, help="Input gene nucleotide FASTA file (.ffn, .fna)")
@click.option("-1", "reads_1", default=None, help="Forward reads (FASTQ)")
@click.option("-2", "reads_2", default=None, help="Reverse reads (FASTQ)")
@click.option("-db", "--database", default=None, help="Path to database folder")
@click.option("-o", "--output", default=None, help="Output directory path")
@click.option("--marker", default=None, help="Comma-separated list of markers (e.g. McrA,PmoA). Default: all")
@click.option("--tree", is_flag=True, default=False, help="Generate phylogenetic tree using FastTree")
def gene(prot, nucl, reads_1, reads_2, database, output, marker, tree):
    """Run functional gene-based profiling pipeline."""
    # Logic to be implemented.
    # We will likely import and call a function from methanohunt.gene here.
    click.echo(f"Running gene pipeline for {prot}...")
    
    # Placeholder for logic
    from .gene import run_gene_pipeline
    run_gene_pipeline(prot, nucl, reads_1, reads_2, database, output, marker, tree)


@cli.command()
@click.option("--target-dir", default=None, help="Directory to install papara (default: python binary dir or ~/.methanohunt/bin)")
def setup(target_dir):
    """Download and setup external dependencies (PaPaRa)."""
    # 1. Determine install directory
    if target_dir:
        install_dir = os.path.abspath(target_dir)
    else:
        # Try to use the directory where python executable is (often within conda env)
        python_bin_dir = os.path.dirname(sys.executable)
        if os.access(python_bin_dir, os.W_OK):
            install_dir = python_bin_dir
        else:
            # Fallback to user home
            install_dir = os.path.join(os.path.expanduser("~"), ".methanohunt", "bin")
            click.echo(f"System bin directory not writable, installing to {install_dir}")

    if not os.path.exists(install_dir):
        os.makedirs(install_dir, exist_ok=True)

    click.echo(f"Installing dependencies to: {install_dir}")

    # 2. Download Papara
    url = "https://cme.h-its.org/exelixis/resource/download/software/papara_nt-2.5-static_x86_64.tar.gz"
    tar_path = os.path.join(install_dir, "papara.tar.gz")
    
    click.echo(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, tar_path)
    except Exception as e:
        click.echo(f"Error downloading: {e}", err=True)
        sys.exit(1)

    # 3. Extract and Setup
    click.echo("Extracting...")
    try:
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=install_dir)
        
        # Binary inside tar is 'papara_static_x86_64'
        extracted_bin = os.path.join(install_dir, "papara_static_x86_64")
        final_bin = os.path.join(install_dir, "papara")
        
        if os.path.exists(extracted_bin):
            if os.path.exists(final_bin):
                os.remove(final_bin)
            os.rename(extracted_bin, final_bin)
            # Make executable
            os.chmod(final_bin, 0o755)
            click.echo(f"Successfully installed 'papara' to {final_bin}")
        else:
             # Check if it was extracted as 'papara' directly (unlikely based on request but good safety)
             if not os.path.exists(final_bin):
                 click.echo("Error: Could not find extracted binary 'papara_static_x86_64'", err=True)
                 sys.exit(1)

    except Exception as e:
        click.echo(f"Error extracting: {e}", err=True)
        sys.exit(1)
    finally:
        if os.path.exists(tar_path):
            os.remove(tar_path)

    # 4. PATH Check
    if install_dir not in os.environ["PATH"]:
        click.echo("\nNOTE: The installation directory is NOT in your PATH.")
        click.echo(f"Please add the following to your shell config file (e.g. .bashrc, .zshrc):")
        click.echo(f"  export PATH=\"{install_dir}:$PATH\"")
    else:
        click.echo("\nInstallation directory is in PATH. verify with 'papara -h'")


if __name__ == "__main__":
    cli()
