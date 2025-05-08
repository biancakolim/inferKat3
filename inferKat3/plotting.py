
def plot_hmm_cnv_heatmap(
    adata,
    cell_type: str,
    chromosome: str,
    obsm_key: str = "hmm_cnv_states",
    groupby: str = "simulated_cnvs",
    figsize=(12, 6),
    cmap="coolwarm"
):
    '''
    Plots a heatmap of HMM CNV calls for a specific cell type and chromosome.

    Parameters:
    - adata: AnnData object with hmm_cnv_states in .obsm
    - cell_type: cell type to subset (e.g. 'CD14 monocyte')
    - chromosome: chromosome to visualize (e.g. '22' or 'chr22')
    - obsm_key: key in adata.obsm with HMM CNV calls
    - groupby: column in adata.obs to use for row annotation (e.g. 'simulated_cnvs')
    '''

    # Subset to cells of interest
    ad_sub = adata[adata.obs["cell_type"] == cell_type].copy()

    # Pull out bin-level CNV calls
    cnv_df = ad_sub.obsm[obsm_key]

    # Filter columns (bins) for the specified chromosome
    chrom_bins = [col for col in cnv_df.columns if col.startswith(f"{chromosome}:") or col.startswith(f"chr{chromosome}:")]
    if not chrom_bins:
        raise ValueError(f"No bins found for chromosome {chromosome} in {obsm_key}.")
    cnv_df = cnv_df[chrom_bins]

    # Optionally convert categorical states to numeric for heatmap
    state_map = {"loss": -1, "neutral": 0, "gain": 1}
    cnv_numeric = cnv_df.replace(state_map)

    # Optional: sort by group
    if groupby in ad_sub.obs.columns:
        sorted_idx = ad_sub.obs[groupby].sort_values().index
        cnv_numeric = cnv_numeric.loc[sorted_idx]
        row_colors = pd.Categorical(ad_sub.obs.loc[sorted_idx, groupby]).codes
    else:
        row_colors = None

    # Plot heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(
        cnv_numeric,
        cmap=cmap,
        center=0,
        xticklabels=True,
        yticklabels=False,
        cbar_kws={'label': 'CNV State'},
    )
    plt.title(f"HMM CNV Heatmap: {cell_type} cells, Chr{chromosome}")
    plt.xlabel("Genomic Bins")
    plt.ylabel("Cells")
    plt.tight_layout()
    plt.show()


def plot_heatmap_from_adata(adata, use_layer=None, log=True, standard_scale=None):
    group_key = 'copykat_prediction'
    if group_key not in adata.obs:
        raise ValueError(f"'{group_key}' not found in adata.obs.")

    sc.pl.heatmap(
        adata,
        var_names=adata.var_names,
        groupby=group_key,
        layer=use_layer,
        log=log,
        standard_scale=standard_scale,
        cmap='RdBu_r',
        dendrogram=False
    )

# Create cell-level summaries
def summarize_state(row, target):
    # return ", ".join([region for region, call in zip(bin_coords, row) if call == target])

    adata.obs["hmm_gain_regions"] = cnv_df.apply(lambda x: summarize_state(x, "gain"), axis=1)
    adata.obs["hmm_loss_regions"] = cnv_df.apply(lambda x: summarize_state(x, "loss"), axis=1)

    return adata
