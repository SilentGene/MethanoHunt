import glob
import pandas as pd
from natsort import natsorted
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import os
import sys

def load_database(db_file):
    """
    Load MethanoHunt database.
    Extract keyword from the last '__' in GTDB_taxonomy for matching.
    """
    if not os.path.exists(db_file):
        raise FileNotFoundError(f"Database file not found at '{db_file}'.")
        
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
    # If input_patterns is a string (single pattern), make it a list
    if isinstance(input_patterns, str):
        input_patterns = [input_patterns]
        
    for p in input_patterns:
        # Expand glob patterns
        matches = glob.glob(p)
        files.extend(matches)
    
    return sorted(list(set(files))) # Remove duplicates


def get_sample_name(tax_file):
    """
    Extract sample name from the first column of singleM tax.tsv file.
    Assumes sample name is the value in the first row, first column.
    """
    try:
        df = pd.read_csv(tax_file, sep="\t")
        if df.empty:
            return os.path.basename(tax_file).replace(".tax.tsv", "")
        sample = df.iloc[0, 0]
        # if sample contains '_1' at the end, remove it
        if str(sample).endswith("_1"):
            sample = str(sample)[:-2]
        return str(sample)
    except Exception as e:
        print(f"Warning: Could not read sample name from {tax_file}: {e}")
        return os.path.basename(tax_file)


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
        try:
            df = pd.read_csv(tax_file, sep="\t")
            # singleM/standard tax files usually have headers or specific columns.
            # The original script enforced columns=["sample", "coverage", "taxonomy"] logic implicitly if read without header?
            # Original script: df = pd.read_csv(tax_file, sep="\t"); df.columns = ["sample", "coverage", "taxonomy"]
            # Let's verify input format reliability. Assuming it works as per original script.
            if len(df.columns) < 3:
                 # Try reading without header if columns are few, or assume it has header.
                 # Original script logic suggests it just overwrites columns.
                 pass

            # Standardize columns for processing
            # note: blindly assigning columns can be risky if input varies, but we follow legacy logic
            df.columns = ["sample", "coverage", "taxonomy"] 

            sample_name = get_sample_name(tax_file)
            sample_names.append(sample_name)

            total_cov = df["coverage"].sum()
            if total_cov == 0:
                print(f"Warning: Total coverage is zero in file: {tax_file}. Skipping.")
                results[sample_name] = 0.0
                continue

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
        
        except Exception as e:
            print(f"Error processing file {tax_file}: {e}")

    # Sort sample columns using natsorted
    sorted_samples = natsorted(sample_names)
    # Keep original columns (except keyword) at the start
    base_cols = [c for c in db.columns if c not in ("keyword", "Exception_taxonomy_list")]
    
    # Ensure all sorted_samples exist in results (in case of skip)
    valid_samples = [s for s in sorted_samples if s in results.columns]

    return results[base_cols + valid_samples]


def generate_stacked_bar_chart(result, output_prefix):
    """
    Generate an interactive grouped stacked bar chart using Plotly.
    X-axis: sample names (grouped by Classification)
    Y-axis: relative abundance (%)
    Stacked bars: colored by keyword (GTDB_taxonomy)
    
    Saves HTML version using the provided prefix.
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
    
    # Set Sample as categorical with natural sorted order
    all_samples_sorted = natsorted(plot_df["Sample"].unique())
    plot_df["Sample"] = pd.Categorical(plot_df["Sample"], categories=all_samples_sorted, ordered=True)
    
    # Sort by classification (now respects DB order) first, then by sample (now uses natural sort)
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
    html_file = f"{output_prefix}.html"
    html_parts = []
    for idx, (_, fig_c) in enumerate(figs):
        include_js = "cdn" if idx == 0 else False
        html_parts.append(pio.to_html(fig_c, include_plotlyjs=include_js, full_html=False))
    html_content = "\n".join(html_parts)
    with open(html_file, "w", encoding="utf-8") as f:
        overall_title = "MethanoHunt Taxonomic Profiling Report"
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
            "<a href='https://github.com/SilentGene/MethanoHunt' style='color: #666; text-decoration: none;'>&copy; 2025 Heyu Lin MethanoHunt</a>"
            "</div>"
        )
        f.write(footnote)
        f.write("\n</body></html>")
    print(f"Saved interactive chart to {html_file}")


def run_taxonomy(input_patterns, database_path, output_dir):
    """
    Main entry point for taxonomy analysis.
    """
    if not database_path:
        # Fallback to default in package
        current_dir = os.path.dirname(os.path.abspath(__file__))
        # database is now in methanohunt/database/methanohunt_taxonomy_db.tsv
        # This file is in methanohunt/taxonomy.py
        database_path = os.path.join(current_dir, "database", "methanohunt_taxonomy_db.tsv")
    
    if not os.path.exists(database_path):
        raise FileNotFoundError(f"Database file not found at {database_path}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
        
    db = load_database(database_path)
    tax_files = load_singlem_files(input_patterns)
    print(f"Found {len(tax_files)} singleM tax.tsv files to process.")
    
    if not tax_files:
        print("No tax.tsv files found matching the pattern.")
        return

    result = compute_abundances(db, tax_files)

    # Define output paths with fixed names inside the output directory
    output_prefix = os.path.join(output_dir, "methanohunt_taxonomy")
    tsv_file = f"{output_prefix}.tsv"
    
    result.to_csv(tsv_file, sep="\t", index=False)
    print(f"Saved relative abundance results to {tsv_file}")
    
    generate_stacked_bar_chart(result, output_prefix)
