# Segmentation and predicting CNV calls

def filter_by_length(state_seq, min_len=3):
    """Suppress CNV calls shorter than `min_len` consecutive bins."""
    filtered = state_seq.copy()
    current_state = state_seq[0]
    start = 0
    for i in range(1, len(state_seq)):
        if state_seq[i] != current_state:
            if current_state != "neutral" and (i - start) < min_len:
                filtered[start:i] = ["neutral"] * (i - start)
            current_state = state_seq[i]
            start = i
    # Check the last run
    if current_state != "neutral" and (len(state_seq) - start) < min_len:
        filtered[start:] = ["neutral"] * (len(state_seq) - start)
    return filtered


def run_hmm_cnv_segmentation(
    adata,
    bin_size: int = 100,
    n_states: int = 3,
    use_layer: str = "smoothed",  # Fallback to X_cnv if unavailable
    obsm_key: str = "hmm_cnv_states",
    min_cnv_length: int = 3  # <--- New parameter for CNV filtering
):
    """
    Performs HMM-based CNV segmentation on smoothed gene expression data and
    stores bin-level CNV states and summary region calls in the AnnData object.

    Parameters:
    - adata: AnnData object with gene-level expression and genome annotations.
    - bin_size: Number of genes per genomic bin.
    - n_states: Number of HMM states (typically 3: loss, neutral, gain).
    - use_layer: Use 'smoothed' layer if it exists; fallback to 'X_cnv'.
    - obsm_key: Key under adata.obsm where bin-level CNV state calls will be saved.
    - min_cnv_length: Minimum length (in bins) to call a gain/loss.

    Returns:
    - Modified AnnData object with CNV states in `.obsm[obsm_key]` and summaries in `.obs`.
    """

    # Determine input matrix
    if use_layer in adata.layers:
        X_expr = adata.layers[use_layer]
    elif "X_cnv" in adata.obsm:
        X_expr = adata.obsm["X_cnv"]
    else:
        raise ValueError("No suitable expression data found in layers['smoothed'] or obsm['X_cnv'].")

    # Sort genes by chromosome and position
    var = adata.var.copy()
    var = var.sort_values(["chromosome", "start"])
    adata_sorted = adata[:, var.index].copy()
    valid_indices = [adata_sorted.var_names.get_loc(gene) for gene in var.index]
    X_expr = X_expr[:, valid_indices]

    # Binning
    bins = [
        (i, i + bin_size)
        for i in range(0, adata_sorted.n_vars, bin_size)
    ]

    bin_expr, bin_coords = [], []
    for start, end in bins:
        end = min(end, adata_sorted.n_vars)  # prevent overflow
        bin_vals = X_expr[:, start:end].mean(axis=1)
        coords = var.iloc[start:end]
        region = f"{coords['chromosome'].values[0]}:{coords['start'].min()}-{coords['end'].max()}"
        bin_expr.append(bin_vals)
        bin_coords.append(region)

    X_bins = np.vstack(bin_expr).T  # shape: cells × bins

    # HMM segmentation
    model = hmm.GaussianHMM(n_components=n_states, covariance_type="diag", n_iter=100)
    model.fit(X_bins.T.flatten().reshape(-1, 1))  # train on all cells jointly

    # Assign CNV states per cell
    cnv_calls = []
    for i in range(X_bins.shape[0]):
        obs = X_bins[i].reshape(-1, 1)
        state_seq = model.predict(obs)

        # Map numeric states to labels based on sorted means
        state_means = model.means_.flatten()
        sorted_states = np.argsort(state_means)
        state_map = {sorted_states[0]: "loss", sorted_states[1]: "neutral", sorted_states[2]: "gain"}
        labeled_seq = np.vectorize(state_map.get)(state_seq)

        # Apply length filter
        filtered_seq = filter_by_length(labeled_seq, min_len=min_cnv_length)
        cnv_calls.append(filtered_seq)

    # Save bin-level CNV calls
    cnv_df = pd.DataFrame(cnv_calls, columns=bin_coords, index=adata.obs_names)
    adata.obsm[obsm_key] = cnv_df

    # Create cell-level summaries
    def summarize_state(row, target):
        return ", ".join([region for region, call in zip(bin_coords, row) if call == target])

    adata.obs["hmm_gain_regions"] = cnv_df.apply(lambda x: summarize_state(x, "gain"), axis=1)
    adata.obs["hmm_loss_regions"] = cnv_df.apply(lambda x: summarize_state(x, "loss"), axis=1)

    return adata
