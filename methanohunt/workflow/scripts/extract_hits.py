import sys
import os

# Snakemake script preamble
# input: faa, tbl
# output: faa

def parse_tbl(tbl_file):
    hits = set()
    with open(tbl_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if not parts:
                continue
            # target name is usually col 0
            hits.add(parts[0])
    return hits

def extract_sequences(fasta_file, hits, output_file):
    found = 0
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        header = None
        keep = False
        for line in f_in:
            if line.startswith(">"):
                # Clean header to get ID (first word)
                header = line.strip()
                seq_id = header[1:].split()[0]
                if seq_id in hits:
                    keep = True
                    f_out.write(f">{seq_id}\n")
                    found += 1
                else:
                    keep = False
            else:
                if keep:
                    f_out.write(line)
    print(f"Extracted {found} sequences to {output_file}")

if __name__ == "__main__":
    hits = parse_tbl(snakemake.input.tbl)
    extract_sequences(snakemake.input.faa, hits, snakemake.output.faa)
