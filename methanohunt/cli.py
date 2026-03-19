import click
import os
import sys
import shutil
import tarfile
import urllib.request
from .taxonomy import run_taxonomy
from .annotate import run_annotate

# from .gene import run_gene_pipeline  # To be implemented

@click.group()
def cli():
    """MethanoHunt: A pipeline for profiling methane cyclers."""
    pass

@cli.command()
@click.option("-i", "--input", "input_files", multiple=True, required=False, help="Input tax.tsv files (supports glob patterns)")
@click.option("-db", "--database", default=None, help="Path to MethanoHunt database file")
@click.option("-o", "--output", required=True, help="Output directory path (will be created if not exists)")
@click.option("-g", "--group", "group_file", default=None, help="Tab-separated file for grouping samples (Sample\\tGroup)")
@click.argument("extra_files", nargs=-1)
def taxonomy(input_files, database, output, group_file, extra_files):
    """Run taxonomy-based profiling.
    
    Accepts input files via -i/--input or as positional arguments (e.g. *.tsv).
    """
    try:
        # Combine inputs from flag and positional args
        all_inputs = list(input_files) + list(extra_files)
        
        if not all_inputs:
             raise click.UsageError("No input files provided. Use -i or positional arguments.")

        # run_taxonomy expects list of patterns or files.
        run_taxonomy(all_inputs, database, output, group_file)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


def parse_greedy_args(args, flags):
    """
    Greedily collect arguments following specific flags until the next flag.
    Handles multiple occurrences by concatenating lists? 
    Or assuming single occurrence for now as per simple glob expansion.
    """
    values = []
    parsing = False
    for arg in args:
        if arg in flags:
            parsing = True
            continue
        if parsing:
            if arg.startswith('-') and len(arg) > 1:
                # Stop if we hit another flag (heuristic: starts with -, not just -)
                # But allow negative numbers? For filenames, usually safe.
                # Assuming typical flag usage.
                parsing = False
                # Don't break, might be another flag we care about later?
                # Actually, if we see another flag, we stop parsing FOR THIS flag.
                # If we see the SAME flag again?
                # Let's simple scan:
                # Find flag index, take next args until next flag.
            else:
                values.append(arg)
    
    # Re-implemented to handle "find flag, take validation"
    collected = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in flags:
            i += 1
            while i < len(args):
                val = args[i]
                if val.startswith('-') and len(val) > 1:
                    # Found next flag
                    break
                collected.append(val)
                i += 1
        else:
            i += 1
    return collected

@cli.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option("--prot", required=True, help="Input protein FASTA file (.faa)")
@click.option("--nucl", default=None, help="Input gene nucleotide FASTA file (.ffn, .fna)")
@click.option("-1", "reads_1", default=None, help="Forward reads (FASTQ). List or glob pattern.")
@click.option("-2", "reads_2", default=None, help="Reverse reads (FASTQ). List or glob pattern.")
@click.option("--mapper", type=click.Choice(['minimap2', 'bwa']), default='minimap2', help="Mapper to use (minimap2 or bwa)")
@click.option("-db", "--database", default=None, help="Path to database folder")
@click.option("-o", "--output", default=None, help="Output directory path")
@click.option("--marker", default=None, help="Comma-separated list of markers (e.g. McrA,PmoA). Default: all")
@click.option("--tree", is_flag=True, default=False, help="Generate phylogenetic tree using FastTree")
@click.option("--keep-bam", is_flag=True, default=False, help="Keep BAM files after run (default: delete)")
@click.option("-t", "--threads", default=4, help="Number of threads to use (default: 4)")
def gene(prot, nucl, reads_1, reads_2, mapper, database, output, marker, tree, keep_bam, threads):
    """Run functional gene-based profiling pipeline."""
    
    # Parse greedy arguments overrides
    # We inspect sys.argv to find the actual lists for -1 and -2
    # Be careful to exclude the command itself 'gene'
    r1_list = parse_greedy_args(sys.argv, ['-1', '--reads_1'])
    r2_list = parse_greedy_args(sys.argv, ['-2', '--reads_2'])

    # If greedy parsing found nothing, fall back to Click's value (which might be just the first file)
    # But wait, if Click parsed "-1 a b", it assigned "a" to reads_1 and "b" to extra_args (if allowed).
    # If we use greedy parsing, we get "a", "b".
    # So if list is not empty, use it.
    
    final_reads_1 = r1_list if r1_list else ([reads_1] if reads_1 else None)
    final_reads_2 = r2_list if r2_list else ([reads_2] if reads_2 else None)
    
    click.echo(f"Running gene pipeline for {prot}...")
    
    from .gene import run_gene_pipeline
    run_gene_pipeline(prot, nucl, final_reads_1, final_reads_2, mapper, database, output, marker, tree, keep_bam, threads)


@cli.command()
@click.option("-i", "--input", "input_file", required=True, help="Input TSV file to annotate")
@click.option("-c", "--column", "column_name", required=True, help="Column name containing the GTDB taxonomy strings")
@click.option("-db", "--database", default=None, help="Path to MethanoHunt taxonomy database file")
@click.option("-o", "--output", "output_file", required=True, help="Output TSV file path")
def annotate(input_file, column_name, database, output_file):
    """Annotate user TSV file with MethanoHunt functional classifications."""
    click.echo(f"Annotating {input_file} using column '{column_name}'...")
    try:
        run_annotate(input_file, column_name, output_file, database)
        click.echo(f"Annotation complete. Results saved to {output_file}")
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


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
