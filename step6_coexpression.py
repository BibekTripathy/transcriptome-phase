import pandas as pd
import numpy as np
import os

# Phase 5 - Step 6 & 7: Build Coexpression Network and Edge List
# This script computes gene-gene relationships (WGCNA-like) from expression.

def run_coexpression(s1_norm, s2_norm, out_file):
    """
    Implements Step 6: Coexpression Network.
    Logic:
    1. Merge normalized expression from S1 and S2.
    2. Compute Pearson Correlation Matrix across all samples.
    3. Convert to Edge List (Adjacency).
    """
    print("--- Running Step 6: Coexpression Network ---")
    
    # Load normalized data
    df1 = pd.read_csv(s1_norm, sep="\t", index_col=0)
    df2 = pd.read_csv(s2_norm, sep="\t", index_col=0)
    
    # Pre-processing: Average duplicate genes in each study
    if df1.index.duplicated().any():
        print(f"  S1: Averaging {df1.index.duplicated().sum()} duplicates.")
        df1 = df1.groupby(df1.index).mean()
    if df2.index.duplicated().any():
        print(f"  S2: Averaging {df2.index.duplicated().sum()} duplicates.")
        df2 = df2.groupby(df2.index).mean()

    # Align on common genes
    common_genes = sorted(list(set(df1.index) & set(df2.index)))
    df_merged = pd.concat([df1.loc[common_genes], df2.loc[common_genes]], axis=1)
    print(f"  Merged {df_merged.shape[1]} samples across {len(common_genes)} common genes.")
    
    # Step 6 sub-step: Compute Pearson Correlation
    # (Using a subset of high-variance genes to speed up for this example, or all if feasible)
    # Let's take the top 2000 most variable genes to keep it efficient and relevant.
    print("  Selecting top 2000 high-variance genes for network...")
    gene_vars = df_merged.var(axis=1)
    top_genes = gene_vars.sort_values(ascending=False).head(2000).index
    df_top = df_merged.loc[top_genes].T # Transpose: samples as rows, genes as columns
    
    corr_matrix = df_top.corr(method="pearson")
    
    # Step 6 sub-step: Convert to Edge List
    # We keep absolute correlation as weight (as suggested in the toy example in doc)
    print("  Converting correlation matrix to edge list...")
    edges = corr_matrix.stack().reset_index()
    edges.columns = ["Gene1", "Gene2", "Weight"]
    
    # Filter: 
    # 1. Remove self-loops
    # 2. Keep only Gene1 < Gene2 to avoid duplicates
    # 3. Threshold at 0.5 for a clean network
    edges = edges[edges["Gene1"] < edges["Gene2"]]
    edges["Weight"] = edges["Weight"].abs()
    edges = edges[edges["Weight"] > 0.5]
    
    # Step 6 & 7 output: Save coexpression network
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    edges.to_csv(out_file, sep="\t", index=False)
    print(f"  Step 6/7 Complete: Saved {len(edges)} edges to {out_file}\n")

# Execute Coexpression
run_coexpression("results/transcriptome/S1/expression_normalized.tsv", "results/transcriptome/S2/expression_normalized.tsv", "results/networks/coexp_wgcna_AD.tsv")
