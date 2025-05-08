
def smooth_matrix(expr, gaussian_sigma=2, use_kalman=False):
    '''
    Fast Gaussian smoothing along the gene axis (rows).
    expr: np.ndarray (genes × cells)
    Returns smoothed matrix (genes × cells)
    '''
    # Transpose to cells × genes for efficient row access
    expr_T = expr.T

    # Fast Gaussian smoothing along genes
    smoothed = gaussian_filter(expr_T, sigma=(0, gaussian_sigma))  # only smooth across genes

    # Center each row (cell)
    smoothed -= smoothed.mean(axis=1, keepdims=True)

    if use_kalman:
        # Optional Kalman polish, slower
        from pykalman import KalmanFilter

        def kalman_polish(x, dV=0.5, dW=0.01):
            kf = KalmanFilter(
                transition_matrices=[1],
                observation_matrices=[1],
                initial_state_mean=x[0],
                initial_state_covariance=1,
                observation_covariance=dV,
                transition_covariance=dW
            )
            smoothed_state_means, _ = kf.smooth(x)
            return smoothed_state_means[:, 0] - np.mean(smoothed_state_means[:, 0])

        from joblib import Parallel, delayed
        smoothed = Parallel(n_jobs=-1)(
            delayed(kalman_polish)(smoothed[i]) for i in range(smoothed.shape[0])
        )
        smoothed = np.vstack(smoothed)

    return smoothed.T  # Return to original shape: genes × cells


def synthetic_baseline(norm_mat: numpy.ndarray, gene_names: list[str], cell_names: list[str], min_cells: int = 10, max_k: int = 6, seed: int = 123):
    '''
    Estimate synthetic baseline using intra-normal GMM-like approach.

    Parameters:
        norm_mat (np.ndarray): smoothed matrix (genes × cells)
        gene_names (list): gene identifiers (rows of norm_mat)
        cell_names (list): cell identifiers (columns of norm_mat)
        min_cells (int): minimum number of cells per cluster
        max_k (int): initial number of clusters
        seed (int): seed for reproducibility

    Returns:
        expr_relat_df (pd.DataFrame): relative expression matrix (genes × cells)
        syn_df (pd.DataFrame): synthetic baseline (genes × clusters)
        cl_labels (np.ndarray): cluster labels per cell (same order as input)
    '''
    numpy.random.seed(seed)
    norm_df = pd.DataFrame(norm_mat, index=gene_names, columns=cell_names)

    # Compute pairwise distances between cells
    dist = squareform(pdist(norm_df.T, metric='euclidean'))

    # Hierarchical clustering
    linkage_matrix = linkage(dist, method="ward")
    k = max_k
    cl_labels = cut_tree(linkage_matrix, n_clusters=k).flatten()

    # Reduce k until all clusters have min_cells
    while numpy.any(numpy.bincount(cl_labels) < min_cells):
        k -= 1
        if k < 2:
            break
        cl_labels = cut_tree(linkage_matrix, n_clusters=k).flatten()

    # Prepare outputs
    expr_relat = []
    syn_cols = []
    valid_clusters = numpy.unique(cl_labels)

    for i in valid_clusters:
        cell_mask = cl_labels == i
        cluster_data = norm_df.iloc[:, cell_mask]

        if cluster_data.shape[1] < min_cells:
            continue

        # Per-gene standard deviation
        sd = cluster_data.std(axis=1)

        # Sample synthetic baseline from N(0, sd)
        syn_norm = numpy.random.normal(loc=0, scale=sd)
        syn_cols.append(pd.Series(syn_norm, index=norm_df.index, name=f"cluster_{i}"))

        # Subtract from cluster cells
        expr_cluster = cluster_data.subtract(syn_norm, axis=0)
        expr_relat.append(expr_cluster)

    # Concatenate outputs
    expr_relat_df = pd.concat(expr_relat, axis=1)
    syn_df = pd.concat(syn_cols, axis=1)

    return expr_relat_df, syn_df, cl_labels


def convert_to_genomic_bins(segmented_cna, adata, bin_size=220_000, chrom_key="chromosome", start_key="start", verbose=True):
  '''
  Convert a segmented CNA matrix to fixed-size genomic bins (e.g. 220kb).

  Parameters:
  - segmented_cna: np.ndarray, shape (n_cells, n_genes), CNA matrix.
  - adata: AnnData object with gene metadata in adata.var.
  - bin_size: int, size of each genomic bin in base pairs.
  - chrom_key: str, column in adata.var for chromosome.
  - start_key: str, column in adata.var for start positions.
  - verbose: bool, whether to print progress.

  Returns:
  - binned_cna: np.ndarray, shape (n_cells, n_bins), binned CNA matrix.
  - bin_annotations: pd.DataFrame, metadata for each bin.
  '''

  chroms = adata.var[chrom_key]
  starts = adata.var[start_key]

  # Filter out genes with missing annotations
  valid_mask = chroms.notna() & starts.notna()
  valid_chroms = chroms[valid_mask]
  valid_starts = starts[valid_mask]
  valid_expr = segmented_cna[:, valid_mask.values]

  # Initialize structures
  bins = []
  bin_matrix = []

  grouped = pd.DataFrame({
      "chrom": valid_chroms.values,
      "start": valid_starts.values
  })

  grouped["gene_idx"] = np.arange(valid_expr.shape[1])

  bin_expr = []
  bin_meta = []

  for chrom in sorted(grouped["chrom"].unique()):
      sub_df = grouped[grouped["chrom"] == chrom].sort_values("start")
      max_pos = sub_df["start"].max()
      for start in range(0, int(max_pos) + bin_size, bin_size):
          end = start + bin_size
          in_bin = sub_df[(sub_df["start"] >= start) & (sub_df["start"] < end)]
          if in_bin.empty:
              continue
          idx = in_bin["gene_idx"].values
          avg_expr = np.median(valid_expr[:, idx], axis=1)
          bin_expr.append(avg_expr)
          bin_meta.append((chrom, start, end))

  binned_cna = np.vstack(bin_expr).T  # shape: (n_cells, n_bins)
  bin_annotations = pd.DataFrame(bin_meta, columns=["chromosome", "start", "end"])

  return binned_cna, bin_annotations


def adjust_baseline(uber_mat_adj, preN, adata, distance='correlation'):
    '''
    Adjusts the CNA matrix using a baseline derived from diploid cells.

    Parameters:
    - uber_mat_adj: DataFrame of shape (bins × cells), adjusted CNA values
    - preN: list or set of putative diploid cell IDs (column names in uber_mat_adj)
    - distance: distance metric for clustering ("correlation" or "euclidean")

    Returns:
    - mat_adj: numpy array of adjusted matrix (same shape as input)
    - com_pred: Series mapping each cell to 'diploid' or 'aneuploid'
    - hcc: linkage matrix from hierarchical clustering
    '''
    adata = adata[:, uber_mat_adj.index].copy()
    # Step 1: Hierarchical clustering
    dists = pdist(uber_mat_adj.T, metric=distance)
    hcc = linkage(dists, method='ward')
    cluster_labels = fcluster(hcc, 2, criterion='maxclust')
    cluster_series = pd.Series(cluster_labels, index=uber_mat_adj.columns)

    # Step 2: Identify diploid cluster (most overlap with preN)
    proportions = [
        len(set(cluster_series[cluster_series == i].index) & set(preN)) / sum(cluster_series == i)
        for i in [1, 2]
    ]
    diploid_cluster = proportions.index(max(proportions)) + 1
    copykat_prediction = cluster_series.map(lambda x: 'diploid' if x == diploid_cluster else 'aneuploid')

    # Step 3: Subtract diploid mean, center matrix
    diploid_cells = copykat_prediction[copykat_prediction == 'diploid'].index
    diploid_mean = uber_mat_adj[diploid_cells].mean(axis=1)
    results_com_rat = uber_mat_adj.subtract(diploid_mean, axis=0)
    results_com_rat = results_com_rat.subtract(results_com_rat.mean(axis=0), axis=1)

    # Step 4: Identify baseline signal
    diploid_only = results_com_rat[diploid_cells]
    base = diploid_only.mean(axis=1)
    cf_h = diploid_only.std(axis=1)

    # Step 5: Smoothing CNA calls
    mask = (results_com_rat.sub(base, axis=0)).abs().le(0.25 * cf_h, axis=0)
    adj_results = results_com_rat.mask(mask, results_com_rat.mean(axis=0), axis=1)
    mat_adj = adj_results.subtract(adj_results.mean(axis=0), axis=1)

    print(copykat_prediction.value_counts())

    # Step 6: Annotate AnnData
    adata.obs['copykat_prediction'] = adata.obs_names.map(copykat_prediction).astype('category')

    threshold = 1.05  # try a more sensitive threshold than 0.1
    base_min = base.min()
    base_max = base.max()

    if np.any(base > threshold*base_max) or np.any(base < threshold*base_min):
        cnv_class = np.where(base > threshold*base_max, 'gain', np.where(base < threshold*base_min, 'loss', 'neutral'))
    else:
        print(f"No CNVs found with threshold={threshold}. Base range: {base_min} to {base_max}")
        cnv_class = np.array(['neutral'] * len(base))  # fallback

    adata.var['cnv_type'] = cnv_class

    print("Length of cnv_class:", len(cnv_class))
    print("Length of adata.var:", adata.shape[1])
    non_neutral = sum(c != 'neutral' for c in cnv_class)
    print("Non-neutral CNV bins:", non_neutral)

    cnv_regions = [
        {'chromosome': adata.var.at[adata.var_names[i], 'chromosome'],
          'start': adata.var.at[adata.var_names[i], 'start'],
          'end': adata.var.at[adata.var_names[i], 'end'],
          'type': cnv_type}
        for i, cnv_type in enumerate(cnv_class) if cnv_type != 'neutral'
    ]
    adata.uns['cnv_regions'] = cnv_regions

    return mat_adj, copykat_prediction, hcc, adata
