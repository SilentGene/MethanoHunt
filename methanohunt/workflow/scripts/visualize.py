import sys
import pandas as pd
import plotly.express as px
import plotly.io as pio

# Snakemake script
# input: data (abundance or classification)
# output: html

def visualize(input_file, output_file):
    try:
        df = pd.read_csv(input_file, sep="\t")
    except Exception:
        # Create empty HTML if file load fails
        with open(output_file, "w") as f:
            f.write("<html><body><h2>No data available for visualization.</h2></body></html>")
        return

    if df.empty:
        with open(output_file, "w") as f:
            f.write("<html><body><h2>No data found.</h2></body></html>")
        return

    # Determine plot type based on columns
    # Abundance: [Classification, TPM, Sample]
    # Classification only: [gene_id, marker, classification]
    
    if "TPM" in df.columns:
        # Bar chart of TPM by Classification
        fig = px.bar(df, x="Classification", y="TPM", color="Classification", 
                     title="Relative Abundance (TPM) of Methane Cyclers")
    elif "classification" in df.columns:
        # Count chart
        counts = df["classification"].value_counts().reset_index()
        counts.columns = ["Classification", "Count"]
        fig = px.bar(counts, x="Classification", y="Count", color="Classification",
                     title="Count of Detected Methane Cycler Genes")
    else:
        # Generic fallback
        fig = px.bar(title="Unknown Data Format")

    fig.update_layout(height=600)
    
    # Save HTML
    pio.write_html(fig, file=output_file, auto_open=False)

if __name__ == "__main__":
    visualize(snakemake.input.data, snakemake.output.html)
