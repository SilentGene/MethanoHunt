import glob
import pandas as pd
from natsort import natsorted
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
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

    keywords = []
    for tax in db["GTDB_taxonomy"]:
        if "__" in tax:
            keyword = tax.split("__")[-1]
        else:
            keyword = tax
        keywords.append(keyword)
    db["keyword"] = keywords
    return db

def compute_abundances(db, df):
    """
    Compute relative abundance for each methane cycler group.
    """
    results = db.copy()
    sample_names = [col for col in df.columns if col != "taxonomy"]
    
    # Clean database taxonomy strings for exact matching
    db_tax = db["GTDB_taxonomy"].astype(str).str.strip()

    for sample in sample_names:
        try:
            rel_abundances = []
            for idx, tax_A in enumerate(db_tax):
                matched_cov = df[df["taxonomy"] == tax_A][sample].sum()
                
                exception_cov = 0
                for other_idx, tax_B in enumerate(db_tax):
                    if idx != other_idx and tax_A in tax_B:
                        exception_cov += df[df["taxonomy"] == tax_B][sample].sum()
                        
                adjusted_cov = matched_cov - exception_cov
                rel_abundances.append(max(adjusted_cov, 0))
                
            results[sample] = rel_abundances
        except Exception as e:
            print(f"Error processing sample {sample}: {e}")

    # Sort sample columns using natsorted
    sorted_samples = natsorted(sample_names)
    # Keep original columns (except keyword) at the start
    base_cols = [c for c in db.columns if c not in ("keyword",)]
    
    # Ensure all sorted_samples exist in results
    valid_samples = [s for s in sorted_samples if s in results.columns]

    return results[base_cols + valid_samples]

def generate_stacked_bar_chart(result, output_prefix, group_file=None):
    """
    Generate an interactive grouped stacked bar chart using Plotly.
    X-axis: sample names (grouped by Classification)
    Y-axis: relative abundance (%)
    Stacked bars: colored by keyword (GTDB_taxonomy)
    
    If group_file is provided, creates faceted subplots by group.
    
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
    metadata_cols = {"GTDB_taxonomy", "Subgroup", "Classification", "keyword", "Exception_subgroup"}
    sample_cols = [col for col in result.columns if col not in metadata_cols]
    
    if not sample_cols:
        print("Warning: No sample columns found. Skipping chart generation.")
        return
    
    # Prepare data for plotting
    plot_data = []
    colors = px.colors.qualitative.Plotly
    color_map = {}
    
    # Load groupings if provided
    sample_to_group = {}
    groups_ordered = []
    if group_file and os.path.exists(group_file):
        try:
            with open(group_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        s_name = parts[0].strip()
                        g_name = parts[1].strip()
                        sample_to_group[s_name] = g_name
                        if g_name not in groups_ordered:
                            groups_ordered.append(g_name)
            print(f"Loaded {len(sample_to_group)} sample groupings.")
            
            # --- Validation: Ensure all samples are in the group file ---
            missing_samples = [s for s in sample_cols if s not in sample_to_group]
            if missing_samples:
                print("\nError: The following samples were not found in the group file (-g):")
                for s in missing_samples:
                    print(f" - {s}")
                print("\nPlease ensure all samples are listed in your group file (tab-separated: SampleName\tGroupName).")
                sys.exit(1)
                
        except Exception as e:
            if isinstance(e, SystemExit):
                raise
            print(f"Error reading group file: {e}")
            group_file = None # Disable grouping on error
    else:
        group_file = None

    for idx, (_, row) in enumerate(result.iterrows()):
        keyword = row["keyword"]
        if keyword not in color_map:
            color_map[keyword] = colors[len(color_map) % len(colors)]
        
        classification = row["Classification"]
        
        for sample in sample_cols:
            rel_abundance = row[sample]
            
            # Determine group
            group = "Data" # Default if no group file
            if group_file:
                # If sample not in group file, you might want to exclude it or put in "Ungrouped"
                # Here we only include if in group file based on request? 
                # Request says: "subplots contains samples in each group". Implicitly, if not in group, maybe exclude?
                # Let's keep "Ungrouped" to be safe or skip?
                # Usually omitting samples not in group file is safer for clarity.
                if sample in sample_to_group:
                    group = sample_to_group[sample]
                else:
                    continue # Skip samples not in the group file if grouping is active
            
            # Keep zeros so samples remain on the axis even if abundance is zero
            plot_data.append({
                "Sample": sample,
                "Classification": classification,
                "Keyword": keyword,
                "Relative Abundance (%)": rel_abundance,
                "Color": color_map[keyword],
                "Group": group
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

    if group_file:
         # Use natsorted to order groups as requested (e.g., G31, G32, G33...)
         unique_groups = natsorted(plot_df["Group"].unique())
         # Set Group as categorical with natural sorted order for consistency
         plot_df["Group"] = pd.Categorical(plot_df["Group"], categories=unique_groups, ordered=True)
         # Sort by classification first, then group, then sample
         plot_df = plot_df.sort_values(["Classification", "Group", "Sample"])
    else:
        unique_groups = ["Data"]
        # Sort by classification first, then by sample
        plot_df = plot_df.sort_values(["Classification", "Sample"])

    # Build individual figures per classification so each has its own legend and y-axis title
    unique_classifications = db_classification_order
    figs = []  # list of tuples (classification, fig)

    total_samples = len(plot_df["Sample"].unique())
    # Increase base width to accommodate fixed right margin for legend
    base_width = min(1600, max(1100, 1100 + total_samples * 40))

    if group_file:
         unique_groups = groups_ordered
         # Filter groups that actually exist in the data
         existing_groups = plot_df["Group"].unique()
         unique_groups = [g for g in unique_groups if g in existing_groups]
         if not unique_groups: # Fallback
             unique_groups = natsorted(existing_groups)
    else:
        unique_groups = ["Data"]

    for classification in unique_classifications:
        class_df = plot_df[plot_df["Classification"] == classification].copy()
        
         # Create subplots only if we have groups
        if group_file:
             num_groups = len(unique_groups)
             fig_c = make_subplots(
                 rows=1, cols=num_groups,
                 subplot_titles=unique_groups,

                 shared_yaxes=True
             )
        else:
             fig_c = go.Figure()

        # Track added legend items to avoid duplicates across subplots
        added_legend_items = set()

        # Iterate over groups (or just single "Data" group)
        for i, group in enumerate(unique_groups):
            group_df = class_df[class_df["Group"] == group].copy()
            # Ensure we use all samples belonging to this group for the x-axis
            samples_in_group = natsorted(plot_df[plot_df["Group"] == group]["Sample"].unique())
            if not samples_in_group:
                continue

            group_df["Sample"] = pd.Categorical(group_df["Sample"].astype(str), categories=samples_in_group, ordered=True)
            group_df = group_df.sort_values("Sample")
            
            traces_added = 0
            for keyword in sorted(group_df["Keyword"].unique()):
                trace_data = group_df[group_df["Keyword"] == keyword]
                if trace_data["Relative Abundance (%)"].sum() == 0:
                    continue  # Skip keywords with no signal

                show_legend = True
                if keyword in added_legend_items:
                    show_legend = False
                else:
                    added_legend_items.add(keyword)
                
                trace = go.Bar(
                    x=trace_data["Sample"],
                    y=trace_data["Relative Abundance (%)"],
                    name=keyword,
                    marker_color=color_map[keyword],
                    legendgroup=keyword,
                    showlegend=show_legend, 
                    hovertemplate="<b>Sample:</b> %{x}<br>" +
                                 "<b>Classification:</b> " + classification + "<br>" +
                                 f"<b>Keyword:</b> {keyword}<br>" +
                                 "<b>Rel. Abundance:</b> %{y:.2f}%<extra></extra>",
                )
                
                if group_file:
                    fig_c.add_trace(trace, row=1, col=i+1)
                else:
                    fig_c.add_trace(trace)
                traces_added += 1

            if traces_added == 0:
                # Add a dummy trace to force subplot rendering (axes/titles)
                dummy_trace = go.Bar(
                    x=samples_in_group,
                    y=[0] * len(samples_in_group),
                    showlegend=False,
                    hoverinfo='skip',
                    marker_color='rgba(0,0,0,0)'
                )
                if group_file:
                    fig_c.add_trace(dummy_trace, row=1, col=i+1)
                else:
                    fig_c.add_trace(dummy_trace)

        layout_args = dict(
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
            margin=dict(t=60, b=90, l=80, r=350)
        )
        fig_c.update_layout(**layout_args)
        
        if not group_file:
             fig_c.update_xaxes(tickangle=-45)
             fig_c.update_yaxes(title_text="Relative Abundance (%)")
        else:
             # Update all xaxes in subplots
             fig_c.update_xaxes(tickangle=-45)
             fig_c.update_yaxes(title_text="Relative Abundance (%)", col=1)

        figs.append((classification, fig_c))
    
    # --- Grouped Box Plot ---
    box_fig = None
    if group_file:
        # Calculate totals per sample per classification
        # Group: Sample, Classification, Group. Sum: Relative Abundance (%)
        # Note: plot_df has all keywords. We sum them up.
        summary_df = plot_df.groupby(["Sample", "Group", "Classification"], observed=True)["Relative Abundance (%)"].sum().reset_index()
        
        box_fig = px.box(
            summary_df, 
            x="Group", 
            y="Relative Abundance (%)", 
            color="Classification",
            points="all", # show points
            hover_data=["Sample"],
            title="Total Relative Abundance by Group and Classification",
            category_orders={"Group": unique_groups, "Classification": db_classification_order}
        )
        box_fig.update_layout(
             boxmode='group',
             height=600,
             width=base_width,
             margin=dict(t=60, b=90, l=80, r=350)
        )

    # Save HTML with all figures stacked
    html_file = f"{output_prefix}.html"
    html_parts = []
    
    # Add Box Plot if groups are present
    if box_fig:
        figs.append(("Grouped Abundance Summary", box_fig))

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

def run_profile(input_wide, input_long, database_path, output_dir, group_file=None):
    """
    Main entry point for profile analysis.
    """
    if not database_path:
        # Fallback to default in package
        current_dir = os.path.dirname(os.path.abspath(__file__))
        database_path = os.path.join(current_dir, "database", "methanohunt_taxonomy_db.tsv")
    
    if not os.path.exists(database_path):
        raise FileNotFoundError(f"Database file not found at {database_path}")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
        
    db = load_database(database_path)

    if input_long:
        print(f"Loading long format input: {input_long}")
        df = pd.read_csv(input_long, sep="\t", comment="#")
        req_cols = {"sample", "taxonomy", "relative_abundance"}
        if not req_cols.issubset(df.columns):
            raise ValueError(f"Input long table must contain columns: {req_cols}")
        df = df[["sample", "taxonomy", "relative_abundance"]]
        # pivot to wide format
        df = df.pivot_table(index="taxonomy", columns="sample", values="relative_abundance", fill_value=0).reset_index()
    elif input_wide:
        print(f"Loading wide format input: {input_wide}")
        df = pd.read_csv(input_wide, sep="\t", comment="#")
        if df.empty or len(df.columns) < 2:
            raise ValueError("Input wide table must contain at least a taxonomy column and one sample column.")
        # Override first column name to 'taxonomy'
        cols = list(df.columns)
        cols[0] = "taxonomy"
        df.columns = cols
    else:
        raise ValueError("Must provide either input_wide or input_long")
        
    # Format taxonomy column perfectly according to requirements
    if "taxonomy" in df.columns:
        # Remove "Root; " prefix
        df["taxonomy"] = df["taxonomy"].astype(str).str.replace(r"^Root; ", "", regex=True)
        # Replace "|" with "; "
        df["taxonomy"] = df["taxonomy"].str.replace(r"\|", "; ", regex=True)
        # Clean whitespaces
        df["taxonomy"] = df["taxonomy"].str.strip()

    print(f"Processing abundances for {len(df.columns) - 1} samples...")
    result = compute_abundances(db, df)

    # Define output paths with fixed names inside the output directory
    output_prefix = os.path.join(output_dir, "methanohunt_profile")
    tsv_level_1 = f"{output_prefix}_level_1.tsv"
    
    result.to_csv(tsv_level_1, sep="\t", index=False)
    print(f"Saved relative abundance results to {tsv_level_1}")
    
    # Generate level 2 table
    cols_to_drop = ["Subgroup", "GTDB_taxonomy", "Exception_subgroup", "Prefered_name"]
    drop_cols = [c for c in cols_to_drop if c in result.columns]
    level_2_df = result.drop(columns=drop_cols)
    
    if "Classification" in level_2_df.columns:
        level_2_grouped = level_2_df.groupby("Classification").sum(numeric_only=True)
        level_2_transposed = level_2_grouped.T
        level_2_transposed.index.name = "Sample"
        
        tsv_level_2 = f"{output_prefix}_level_2.tsv"
        level_2_transposed.to_csv(tsv_level_2, sep="\t")
        print(f"Saved aggregated relative abundance results to {tsv_level_2}")

    generate_stacked_bar_chart(result, output_prefix, group_file)
