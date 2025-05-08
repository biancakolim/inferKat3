
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


def synthetic_baseline(norm_mat: np.ndarray, gene_names: list[str], cell_names: list[str], min_cells: int = 10, max_k: int = 6, seed: int = 123):
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
    np.random.seed(seed)
    norm_df = pd.DataFrame(norm_mat, index=gene_names, columns=cell_names)

    # Compute pairwise distances between cells
    dist = squareform(pdist(norm_df.T, metric='euclidean'))

    # Hierarchical clustering
    linkage_matrix = linkage(dist, method="ward")
    k = max_k
    cl_labels = cut_tree(linkage_matrix, n_clusters=k).flatten()

    # Reduce k until all clusters have min_cells
    while np.any(np.bincount(cl_labels) < min_cells):
        k -= 1
        if k < 2:
            break
        cl_labels = cut_tree(linkage_matrix, n_clusters=k).flatten()

    # Prepare outputs
    expr_relat = []
    syn_cols = []
    valid_clusters = np.unique(cl_labels)

    for i in valid_clusters:
        cell_mask = cl_labels == i
        cluster_data = norm_df.iloc[:, cell_mask]

        if cluster_data.shape[1] < min_cells:
            continue

        # Per-gene standard deviation
        sd = cluster_data.std(axis=1)

        # Sample synthetic baseline from N(0, sd)
        syn_norm = np.random.normal(loc=0, scale=sd)
        syn_cols.append(pd.Series(syn_norm, index=norm_df.index, name=f"cluster_{i}"))

        # Subtract from cluster cells
        expr_cluster = cluster_data.subtract(syn_norm, axis=0)
        expr_relat.append(expr_cluster)

    # Concatenate outputs
    expr_relat_df = pd.concat(expr_relat, axis=1)
    syn_df = pd.concat(syn_cols, axis=1)

    return expr_relat_df, syn_df, cl_labels