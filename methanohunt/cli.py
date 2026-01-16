import click
import os
import sys
from .taxonomy import run_taxonomy
# from .gene import run_gene_pipeline  # To be implemented

@click.group()
def cli():
    """MethanoHunt: A pipeline for profiling methane cyclers."""
    pass

@cli.command()
@click.option("-i", "--input", "input_files", multiple=True, required=True, help="Input tax.tsv files (supports glob patterns)")
@click.option("-db", "--database", default=None, help="Path to MethanoHunt database file")
@click.option("-o", "--output", required=True, help="Output TSV file path")
def taxonomy(input_files, database, output):
    """Run taxonomy-based profiling."""
    try:
        # Click handles multiple inputs as tuple of strings.
        # run_taxonomy expects list of patterns or files.
        run_taxonomy(list(input_files), database, output)
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


@cli.command()
@click.option("--protein", required=True, help="Input protein FASTA file (.faa)")
@click.option("--nucl", default=None, help="Input gene nucleotide FASTA file (.ffn, .fna)")
@click.option("-1", "reads_1", default=None, help="Forward reads (FASTQ)")
@click.option("-2", "reads_2", default=None, help="Reverse reads (FASTQ)")
@click.option("-db", "--database", default=None, help="Path to database folder")
@click.option("-o", "--output", default=None, help="Output directory path")
def gene(protein, nucl, reads_1, reads_2, database, output):
    """Run functional gene-based profiling pipeline."""
    # Logic to be implemented.
    # We will likely import and call a function from methanohunt.gene here.
    click.echo(f"Running gene pipeline for {protein}...")
    
    # Placeholder for logic
    from .gene import run_gene_pipeline
    run_gene_pipeline(protein, nucl, reads_1, reads_2, database, output)


if __name__ == "__main__":
    cli()
