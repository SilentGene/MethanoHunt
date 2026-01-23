import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio

def plot_rpkg_distribution(class_file, subtype_file, output_html):
    figs = []

    # Helper function to create a figure
    def create_fig(input_file, title):
        try:
            df = pd.read_csv(input_file, sep='\t', index_col='Sample')
        except:
             return None
             
        if df.empty:
            return None
        
        # Structure: Rows=Sample, Cols=Categories.
        # Plotly goes by traces (categories).
        
        fig = go.Figure()
        
        # Get categories (columns)
        categories = df.columns.tolist()
        
        for cat in categories:
            fig.add_trace(go.Bar(
                name=cat,
                x=df.index,
                y=df[cat],
                hovertemplate=f"<b>Sample:</b> %{{x}}<br><b>{title}:</b> {cat}<br><b>RPKG:</b> %{{y:.4f}}<extra></extra>"
            ))
            
        fig.update_layout(
            title=f"RPKG Distribution by {title}",
            barmode='stack',
            xaxis_title="Sample",
            yaxis_title="RPKG",
            legend_title=title,
            hovermode="closest",
            height=700
        )
        return fig

    # Create charts
    fig_class = create_fig(class_file, "Classification")
    if fig_class: figs.append(fig_class)
    
    fig_subtype = create_fig(subtype_file, "Subtype")
    if fig_subtype: figs.append(fig_subtype)
    
    # Save to HTML
    with open(output_html, 'w', encoding="utf-8") as f:
        f.write("<html><head><meta charset='utf-8'><title>MethanoHunt Report</title></head><body>\n")
        f.write("<h1 style='font-family: Arial, sans-serif; text-align: left;'>MethanoHunt Gene Profiling Report</h1>\n")
        f.write("<p style='font-family: Arial, sans-serif; font-size: 14px;'>This report displays the RPKG (Reads Per Kilobase per Genome equivalent) abundance of different functional gene categories.</p>\n")
        
        if not figs:
            f.write("<p>No data available for visualization.</p>")
        else:
            for i, fig in enumerate(figs):
                include_js = "cdn" if i == 0 else False
                f.write(pio.to_html(fig, include_plotlyjs=include_js, full_html=False))
                f.write("<br>")
        
        footnote = (
            "<div style='margin: 30px 10px 20px 10px; padding-top: 10px; border-top: 1px solid #ccc; "
            "font-family: Arial, sans-serif; font-size: 0.9em; color: #666; font-style: italic;'>"
            "<br>"
            "<a href='https://github.com/SilentGene/MethanoHunt' style='color: #666; text-decoration: none;'>&copy; 2025 Heyu Lin MethanoHunt</a>"
            "</div>"
        )
        f.write(footnote)
        
        f.write("</body></html>")

if __name__ == "__main__":
    class_tsv = snakemake.input.class_tsv
    subtype_tsv = snakemake.input.subtype_tsv
    html_out = snakemake.output.html
    
    plot_rpkg_distribution(class_tsv, subtype_tsv, html_out)
