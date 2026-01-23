import pandas as pd
import os

def aggregate_rpkg(rpkg_files, classification_file, out_combined, out_class, out_subtype):
    # 1. Load classification mapping
    cls_df = pd.read_csv(classification_file, sep='\t')
    
    if 'name' not in cls_df.columns:
        raise ValueError("Classification file missing 'name' column.")
    
    if 'classification' not in cls_df.columns: cls_df['classification'] = 'Unknown'
    if 'subtype' not in cls_df.columns: cls_df['subtype'] = 'Unknown'

    gene_to_class = dict(zip(cls_df['name'], cls_df['classification']))
    gene_to_subtype = dict(zip(cls_df['name'], cls_df['subtype']))

    # 2. Combine Samples
    combined_df = pd.DataFrame() 
    
    for f in rpkg_files:
        sample_name = os.path.basename(f).replace('_rpkg.tsv', '')
        # Read: GeneID, RPKG
        df = pd.read_csv(f, sep='\t')
        
        if df.empty:
             continue
             
        df = df.set_index('GeneID')
        df.columns = [sample_name]
        
        if combined_df.empty:
            combined_df = df
        else:
            combined_df = combined_df.join(df, how='outer')
            
    combined_df = combined_df.fillna(0)
    combined_df.index.name = 'GeneID'
    
    # Save combined
    combined_df_T = combined_df.T
    combined_df_T.index.name = 'Sample'
    combined_df_T.to_csv(out_combined, sep='\t')

    # 3. Aggregate
    # Map index to Classification -> Sum
    class_groups = combined_df.groupby(lambda x: gene_to_class.get(x, 'Unclassified')).sum()
    class_groups_T = class_groups.T
    class_groups_T.index.name = 'Sample'
    class_groups_T.to_csv(out_class, sep='\t')
    
    subtype_groups = combined_df.groupby(lambda x: gene_to_subtype.get(x, 'Unclassified')).sum()
    subtype_groups_T = subtype_groups.T
    subtype_groups_T.index.name = 'Sample'
    subtype_groups_T.to_csv(out_subtype, sep='\t')

if __name__ == "__main__":
    rpkg_files = snakemake.input.rpkgs
    classification_file = snakemake.input.classification
    
    aggregate_rpkg(rpkg_files, classification_file, 
                  snakemake.output.combined, 
                  snakemake.output.class_tsv, 
                  snakemake.output.subtype_tsv)
