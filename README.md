# MethanoHunt

A pipeline to profile methane cyclers from taxonomic profiling data or functional marker genes.

## Overview

MethanoHunt provides two main modules:
1.  **Taxonomy**: Summarizes relative abundance of methane cyclers from taxonomic profiles (e.g. singleM).
2.  **Gene**: A pipeline to detect, classify, and quantify methane cycling marker genes (mcrA, pmoA, mmoA) from protein and nucleotide sequences.

## Installation

### Dependencies
MethanoHunt requires:
*   Python 3.8+
*   Snakemake > 9.0
*   HMMER
*   pplacer
*   Minimap2, Samtools
*   Python packages: `pandas`, `natsort`, `plotly`, `click`, `pysam`

### Conda Environment
We recommend using Conda/Mamba:

```bash
conda env create -f environment.yaml
conda activate methanohunt
pip install .
```

## Usage

### Taxonomy Module
Analyze taxonomic abundance tables.

```bash
methanohunt taxonomy -i sample1.tax.tsv sample2.tax.tsv -o methanohunt_results.tsv
```

*   `-i`: Input taxonomy tables (supports glob patterns).
*   `-o`: Output TSV file. Also generates an HTML report.
*   `-db`: (Optional) Custom database path.

### Gene Module
Analyze functional genes from protein sequences.

```bash
methanohunt gene --protein assembly.faa [-o output_dir] [-db database_dir]
```

**With abundance calculation (requires nucleotide genes and reads):**

```bash
methanohunt gene \
    --protein assembly.faa \
    --nucl assembly.genes.ffn \
    -1 hits_R1.fq.gz -2 hits_R2.fq.gz \
    -o my_results \
    -db /path/to/database_folder
```

*   `--protein`: Input protein FASTA (.faa).
*   `--nucl`: (Optional) Corresponding gene nucleotide FASTA (.ffn), required for abundance.
*   `-1`, `-2`: (Optional) Reads for abundance quantification.
*   `-db`: Database folder containing `.hmm` and `.refpkg` files for mcrA, pmoA, mmoA.
*   `-o`: Output directory.

## Output

### Gene Module Output
*   `classification.tsv`: Detected genes and their functional classification.
*   `abundance.tsv`: TPM abundance of classified genes (if reads provided).
*   `MethanoHunt_report.html`: Interactive visualization.
*   `hmm/`, `pplacer/`: Intermediate results.

## Database

The `taxonomy` module uses `methanohunt_db.tsv` (included).
The `gene` module requires a database folder with:
*   `mcrA.hmm`, `pmoA.hmm`, `mmoA.hmm`
*   `mcrA.refpkg`, `pmoA.refpkg`, `mmoA.refpkg`

...🧙‍♂️🧬
