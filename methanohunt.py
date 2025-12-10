#!/usr/bin/env python3

"""
Summarize relative abundance of methane cyclers from taxonomy tables.

Usage: python methanohunt.py -i *.tax.tsv [-db methanohunt_db.tsv] -o methanohunt_output.tsv

-i/--input: Input taxonomy tables (in long-format taxonomic abundance table from singleM or other profilers). Supports glob patterns.
-db/--database: Path to MethanoHunt database file (default: 'methanohunt_db.tsv').
-o/--output: Output TSV file to save relative abundance results. Also generates an interactive HTML chart with the same base name.

The script will look for 'methanohunt_db.tsv' in the script directory by default if -db is not provided.
"""

import argparse
import glob
import pandas as pd
from natsort import natsorted
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import os

def load_database(db_file):
    """
    Load MethanoHunt database.
    Extract keyword from the last '__' in GTDB_taxonomy for matching.
    """
    db = pd.read_csv(db_file, sep="\t")

    # Ensure Exception_taxonomy exists and normalize to list for later subtraction
    if "Exception_taxonomy" not in db.columns:
        db["Exception_taxonomy"] = ""
    db["Exception_taxonomy_list"] = (
        db["Exception_taxonomy"]
        .fillna("")
        .apply(lambda x: [t.strip() for t in str(x).split(",") if t.strip()])
    )

    keywords = []
    for tax in db["GTDB_taxonomy"]:
        if "__" in tax:
            keyword = tax.split("__")[-1]
        else:
            keyword = tax
        keywords.append(keyword)
    db["keyword"] = keywords
    return db


def load_singlem_files(input_patterns):
    """
    Load all sample.tax.tsv files using glob patterns.
    """
    files = []
    for p in input_patterns:
        files.extend(glob.glob(p))
    return sorted(files)


def get_sample_name(tax_file):
    """
    Extract sample name from the first column of singleM tax.tsv file.
    Assumes sample name is the value in the first row, first column.
    """
    df = pd.read_csv(tax_file, sep="\t")
    sample = df.iloc[0, 0]
    # if sample contains '_1' at the end, remove it
    if sample.endswith("_1"):
        sample = sample[:-2]
    return str(sample)


def compute_abundances(db, tax_files):
    """
    Compute relative abundance for each methane cycler group.
    """
    results = db.copy()
    sample_names = []

    # Map Subgroup values to their keyword so exceptions can look up coverage quickly
    subgroup_to_keyword = {}
    if "Subgroup" in results.columns:
        subgroup_to_keyword = dict(zip(results["Subgroup"], results["keyword"]))

    for tax_file in tax_files:
        df = pd.read_csv(tax_file, sep="\t")
        df.columns = ["sample", "coverage", "taxonomy"]

        sample_name = get_sample_name(tax_file)
        sample_names.append(sample_name)

        total_cov = df["coverage"].sum()
        if total_cov == 0:
            raise RuntimeError(f"Total coverage is zero in file: {tax_file}")

        # Compute coverage per keyword once so exceptions can reuse it
        coverage_by_keyword = {}
        for keyword in results["keyword"]:
            matched_cov = df[df["taxonomy"].str.contains(keyword, na=False)]["coverage"].sum()
            coverage_by_keyword[keyword] = matched_cov

        rel_abundances = []
        for idx, keyword in enumerate(results["keyword"]):
            matched_cov = coverage_by_keyword.get(keyword, 0)

            exception_cov = 0
            if "Exception_taxonomy_list" in results.columns:
                for exc in results.at[idx, "Exception_taxonomy_list"]:
                    exc_keyword = subgroup_to_keyword.get(exc)
                    if exc_keyword:
                        exception_cov += coverage_by_keyword.get(exc_keyword, 0)

            adjusted_cov = matched_cov - exception_cov
            rel_abundance = (max(adjusted_cov, 0) / total_cov) * 100
            rel_abundances.append(rel_abundance)

        results[sample_name] = rel_abundances

    # Sort sample columns using natsorted
    sorted_samples = natsorted(sample_names)
    # Keep original columns (except keyword) at the start
    base_cols = [c for c in db.columns if c not in ("keyword", "Exception_taxonomy_list")]

    return results[base_cols + sorted_samples]


def generate_stacked_bar_chart(result, output_file):
    """
    Generate an interactive grouped stacked bar chart using Plotly.
    X-axis: sample names (grouped by Classification)
    Y-axis: relative abundance (%)
    Stacked bars: colored by keyword (GTDB_taxonomy)
    
    Saves both HTML and JPG versions.
    """
    if "Classification" not in result.columns:
        print("Warning: 'Classification' column not found. Skipping chart generation.")
        return
    
    if "GTDB_taxonomy" not in result.columns:
        print("Warning: 'GTDB_taxonomy' column not found. Skipping chart generation.")
        return
    
    # Extract keywords from GTDB_taxonomy if not already present
    if "keyword" not in result.columns:
        keywords = []
        for tax in result["GTDB_taxonomy"]:
            if "__" in str(tax):
                keyword = str(tax).split(";")[-1].strip()
            else:
                keyword = str(tax)
            keywords.append(keyword)
        result = result.copy()
        result["keyword"] = keywords
    
    # Identify sample columns (those that are not metadata)
    metadata_cols = {"GTDB_taxonomy", "Subgroup", "Classification", "keyword", "Exception_taxonomy_list", "Exception_taxonomy"}
    sample_cols = [col for col in result.columns if col not in metadata_cols]
    
    if not sample_cols:
        print("Warning: No sample columns found. Skipping chart generation.")
        return
    
    # Prepare data for plotting
    plot_data = []
    colors = px.colors.qualitative.Plotly
    color_map = {}
    
    for idx, (_, row) in enumerate(result.iterrows()):
        keyword = row["keyword"]
        if keyword not in color_map:
            color_map[keyword] = colors[len(color_map) % len(colors)]
        
        classification = row["Classification"]
        
        for sample in sample_cols:
            rel_abundance = row[sample]
            # Keep zeros so samples remain on the axis even if abundance is zero
            plot_data.append({
                "Sample": sample,
                "Classification": classification,
                "Keyword": keyword,
                "Relative Abundance (%)": rel_abundance,
                "Color": color_map[keyword]
            })
    
    if not plot_data:
        print("Warning: No data to plot (all relative abundances are zero).")
        return
    
    plot_df = pd.DataFrame(plot_data)
    
    # Get classification order from original result to preserve database order
    # pd.unique preserves order of appearance
    db_classification_order = [c for c in result["Classification"].unique() if c in plot_df["Classification"].unique()]

    # Set Classification as categorical with prescribed order
    plot_df["Classification"] = pd.Categorical(plot_df["Classification"], categories=db_classification_order, ordered=True)
    
    # Sort by classification (now respects DB order) first, then by sample
    plot_df = plot_df.sort_values(["Classification", "Sample"])

    # Build individual figures per classification so each has its own legend and y-axis title
    unique_classifications = db_classification_order
    figs = []  # list of tuples (classification, fig)

    total_samples = len(plot_df["Sample"].unique())
    # Increase base width to accommodate fixed right margin for legend
    base_width = min(1600, max(1100, 1100 + total_samples * 40))

    for classification in unique_classifications:
        class_df = plot_df[plot_df["Classification"] == classification].copy()
        samples_in_class = natsorted(class_df["Sample"].unique())
        class_df.loc[:, "Sample"] = pd.Categorical(class_df["Sample"], categories=samples_in_class, ordered=True)
        class_df = class_df.sort_values("Sample")

        fig_c = go.Figure()
        for keyword in sorted(class_df["Keyword"].unique()):
            trace_data = class_df[class_df["Keyword"] == keyword]
            if trace_data["Relative Abundance (%)"].sum() == 0:
                continue  # Skip keywords with no signal in this classification
            fig_c.add_trace(
                go.Bar(
                    x=trace_data["Sample"],
                    y=trace_data["Relative Abundance (%)"],
                    name=keyword,
                    marker_color=color_map[keyword],
                    legendgroup=keyword,
                    showlegend=True,  # Always show legend entry even if only one keyword
                    hovertemplate="<b>Sample:</b> %{x}<br>" +
                                 "<b>Classification:</b> " + classification + "<br>" +
                                 f"<b>Keyword:</b> {keyword}<br>" +
                                 "<b>Rel. Abundance:</b> %{y:.2f}%<extra></extra>",
                )
            )

        fig_c.update_layout(
            title=f"{classification}",
            barmode="stack",
            height=420,
            width=base_width,
            hovermode="closest",
            legend=dict(
                title="Taxonomy",
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02
            ),
            # Set a large fixed right margin to accommodate the widest legend
            # This ensures the plot area width remains consistent across all figures
            margin=dict(t=60, b=90, l=80, r=350)
        )
        fig_c.update_xaxes(tickangle=-45, title_text="Sample")
        fig_c.update_yaxes(title_text="Relative Abundance (%)")
        figs.append((classification, fig_c))

    # Save HTML with all figures stacked
    html_file = output_file.replace(".tsv", "") + ".html"
    html_parts = []
    for idx, (_, fig_c) in enumerate(figs):
        include_js = "cdn" if idx == 0 else False
        html_parts.append(pio.to_html(fig_c, include_plotlyjs=include_js, full_html=False))
    html_content = "\n".join(html_parts)
    with open(html_file, "w", encoding="utf-8") as f:
        overall_title = "Relative Abundance of Methane Cyclers"
        f.write("<html><head><meta charset='utf-8'></head><body>\n")
        f.write(f"<h2 style='margin:10px 0 20px 0; font-family:Arial, sans-serif;'>{overall_title}</h2>\n")
        f.write(html_content)
        
        footnote = (
            "<div style='margin: 30px 10px 20px 10px; padding-top: 10px; border-top: 1px solid #ccc; "
            "font-family: Arial, sans-serif; font-size: 0.9em; color: #666; font-style: italic;'>"
            "<strong>Note:</strong> Functional classifications presented here are taxonomy-based and may carry a risk of false positives, "
            "as specific subgroups within a lineage may lack the predicted metabolic potential. "
            "Verification via functional gene analysis is recommended."
            "<br><br>"
            "&copy; 2025 Heyu Lin"
            "</div>"
        )
        f.write(footnote)
        f.write("\n</body></html>")
    print(f"Saved interactive chart to {html_file}")


def main():
    parser = argparse.ArgumentParser(description="Summarize relative abundance of methane cyclers from singleM tax.tsv files.")
    parser.add_argument("-i", "--input", nargs="+", required=True,
                        help="Input tax.tsv files (supports glob, e.g. '*.tax.tsv')")
    parser.add_argument("-db", "--database", required=False,
                        help="Database file methane_cycler.tsv (default: methane_cycler_db.tsv in script directory)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file")
    args = parser.parse_args()

    if args.database:
        db_path = args.database
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        db_path = os.path.join(script_dir, "methane_cycler_db.tsv")

    if not os.path.exists(db_path):
        import sys
        sys.exit(f"Error: Database file not found at '{db_path}'. Please provide it with -db or place 'methane_cycler_db.tsv' in the script directory.")

    db = load_database(db_path)
    tax_files = load_singlem_files(args.input)
    print(f"Found {len(tax_files)} singleM tax.tsv files to process.")

    if len(tax_files) == 0:
        raise RuntimeError("No tax.tsv files found. Please check -i argument.")

    result = compute_abundances(db, tax_files)
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Saved relative abundance results to {args.output}")
    
    # Generate visualization
    generate_stacked_bar_chart(result, args.output)


if __name__ == "__main__":
    main()
