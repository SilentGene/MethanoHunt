# MethanoHunt

A pipeline to profile methane cyclers from taxonomic profiling data or functional marker genes.

## Overview

MethanoHunt provides two main modules:
1.  **Taxonomy**: Summarizes relative abundance of methane cyclers from taxonomic profiles (e.g. singleM).
2.  **Gene**: A pipeline to detect, classify, and quantify methane cycling marker genes (McrA, PmoA, MmoX) from protein and nucleotide sequences.

## Installation

### Dependencies
MethanoHunt requires:
*   Python 3.8+
*   Snakemake > 9.0
*   HMMER3
*   raxml-ng, epa-ng, gappa
*   Minimap2, Samtools
*   Python packages: `pandas`, `natsort`, `plotly`, `click`, `pysam`
*   PaPaRa (install independently)

Use conda to install all depencies except PaPaRa which needs another step.

```bash
git clone https://github.com/SilentGene/MethanoHunt.git
conda env create -f methanohunt.yaml
conda activate methanohunt
methanohunt setup  # to install PaPaRa
```

## Usage

### Taxonomy Module
Analyze taxonomic abundance tables.

```bash
methanohunt taxonomy -i *.tax.tsv -o methanohunt_results
```

*   `-i`: Input taxonomy tables (supports glob patterns).
*   `-o`: Prefix for Output files. It will generate a TSV result and an HTML report.
*   `-db`: (Optional) Custom database path. If not provided, it will use the default database installed along with the pipeline.

### Gene Module
Analyze functional genes from protein sequences.

```bash
methanohunt gene --prot assembly.faa [-o output_dir] [-db database_dir]
```

**With abundance calculation (requires nucleotide genes and reads):**

```bash
methanohunt gene \
    --prot assembly_genes.faa \
    --nucl assembly_genes.ffn \
    -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
    -o my_results \
    --marker McrA,PmoA,MmoX \
    --tree
```

*   `--prot`: Input protein FASTA (.faa).
*   `--nucl`: (Optional) Corresponding gene nucleotide FASTA (.ffn), required for abundance.
*   `-1`, `-2`: (Optional) Reads for abundance quantification.
*   `-o`: Output directory.
*   `--marker`: (Optional) Marker genes to analyze. Default is McrA,PmoA,MmoX.
*   `--tree`: (Optional) Build phylogenetic trees together with reference markers using fasttree. Default is False.
*   `-db`: Database folder. Default is the database folder included in the package.

## Output

### Gene Module Output
*   `methanohunt_gene_classification.tsv`: Detected genes and their functional classification.
*   `abundance.tsv`: TPM abundance of classified genes (if reads provided).
*   `MethanoHunt_report.html`: Abundance visualization (if reads provided).
*   `hmm/`, `hits/`, `placement/`, `classification/`, `tree/`: Intermediate results.

## Database

The required database files for both taxonomy and gene modules are included in the package.


...🧙‍♂️🧬
