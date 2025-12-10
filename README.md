# MethanoHunt

A Python script to summarize the relative abundance of methane cyclers from taxonomic profiling data.

## Overview

MethanoHunt processes taxonomy tables from metagenomic profilers (like singleM) and calculates the relative abundance of different methane metabolism groups. It generates both tabular results and interactive visualizations.

## Installation

Requires Python 3.6+ with dependencies:
```bash
pip install pandas natsort plotly
```

## Usage

```bash
python methanohunt.py -i <input_files> [-db <database>] -o <output_file>
```

### Arguments

- `-i, --input` (required): Input tax.tsv files. Supports glob patterns (e.g., `*.tax.tsv`)
- `-db, --database` (optional): Path to MethanoHunt database file. Default: `methane_cycler_db.tsv` in script directory
- `-o, --output` (required): Output TSV file path. Also generates an interactive HTML chart

### Example

```bash
python methanohunt.py -i singleM_results/*.tax.tsv -o methanohunt_results.tsv
```

## Output

- **TSV file**: Relative abundance (%) for each methane cycler group per sample
- **HTML file**: Interactive stacked bar chart grouped by classification

Here is an example of the output [methanohunt_output.html](example_output.html).

Screenshot of the interactive chart:
![MethanoHunt Chart](example.jpg)

## Database Format

The database TSV should contain:
- `GTDB_taxonomy`: Taxonomic identifier
- `Subgroup`: Metabolic subgroup
- `Classification`: Functional classification
- `Exception_taxonomy` (optional): Comma-separated taxa to exclude

## Notes

Functional classifications are taxonomy-based and may have false positives. Verification with functional gene analysis is recommended.
