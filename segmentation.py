

def segment_cna_clusterwise(
  adata,
  layer="relative",
  window_size=50,
  ks_cutoff=0.2,
  dr_threshold=0.05,
  up_dr=0.98,
  cluster_key="cluster_labels",
  verbose=True
):
  '''
  Cluster-wise segmentation of CNAs using KS tests and sliding windows.

  Parameters:
  - adata: AnnData object with preprocessed gene expression.
  - layer: Layer from which to take the normalized matrix.
  - window_size: Number of genes per window.
  - ks_cutoff: KS statistic threshold to define a breakpoint.
  - dr_threshold: Minimum gene detection rate to include a gene.
  - up_dr: Maximum gene detection rate to filter out highly-expressed ubiquitous genes.
  - cluster_key: Key in adata.obs containing cluster assignments.
  - verbose: Print progress info.

  Returns:
  - segmented_cna: np.ndarray, shape = (n_cells, n_genes), smoothed CNA matrix.
  - breakpoints: list of detected gene breakpoints.
  '''
  X = adata.layers[layer]
  if isinstance(X, np.ndarray):
      expr = X
  else:
      expr = X.toarray()

  # Filter genes by detection rate
  detection_rate = np.sum(expr > 0, axis=0) / expr.shape[0]
  keep_genes = (detection_rate > dr_threshold) & (detection_rate < up_dr)
  expr = expr[:, keep_genes]
  if verbose:
      print(f"Filtered to {expr.shape[1]} genes based on detection rate thresholds.")

  # Sort genes by chromosome and start coordinate
  chroms = adata.var.loc[keep_genes, "chromosome"]
  starts = adata.var.loc[keep_genes, "start"]

  valid_mask = chroms.notna() & starts.notna()
  chroms_valid = chroms[valid_mask]
  starts_valid = starts[valid_mask]
  sort_idx = np.lexsort((starts_valid.values, chroms_valid.values))

  expr = expr[:, sort_idx]
  gene_idx_sorted = np.where(keep_genes)[0][sort_idx]

  clusters = adata.obs[cluster_key].astype(str)
  cluster_labels = sorted(clusters.unique())
  cluster_medians = []

  # Compute median expression per cluster
  for label in cluster_labels:
      mask = clusters == label
      cluster_expr = expr[mask.values, :]
      cluster_medians.append(np.median(cluster_expr, axis=0))

  cluster_medians = np.array(cluster_medians)

  # Use KS to find breakpoints
  n_genes = cluster_medians.shape[1]
  breakpoints = []

  for start in tqdm(range(0, n_genes - 2 * window_size, window_size), disable=not verbose):
      i1, i2 = start, start + window_size
      i3 = i2
      i4 = i2 + window_size
      group1 = cluster_medians[:, i1:i2].flatten()
      group2 = cluster_medians[:, i3:i4].flatten()
      ks_stat, _ = ks_2samp(group1, group2)
      if ks_stat > ks_cutoff:
          breakpoints.append(i2)

  breakpoints = sorted(set([0] + breakpoints + [n_genes]))

  # Apply segmentation: average expression within segments
  segmented = np.zeros_like(expr)
  for i in range(len(breakpoints) - 1):
      start, end = breakpoints[i], breakpoints[i + 1]
      mean_vals = np.mean(expr[:, start:end], axis=1)
      segmented[:, start:end] = mean_vals[:, None]

  # Map back to full gene space
  segmented_full = np.zeros((adata.n_obs, adata.n_vars))
  segmented_full[:, gene_idx_sorted] = segmented

  return segmented_full, breakpoints


# use this or hmm segmentation
def annotate_cnv_regions(
    adata,
    cnv_key="X_cnv",
    threshold=0.3,
    normal_cn=2,
    output_key="inferred_cnvs"
):
    '''
    Annotate cells with CNV events in the format:
    'chr:start-end (CN copy_number)', stored in adata.obs[output_key].

    Parameters:
    - adata: AnnData object with infercnvpy outputs
    - cnv_key: key in adata.obsm for CNV smoothed values (default: 'X_cnv')
    - threshold: minimum deviation from diploid (2) to count as CNV (default: 0.3)
    - normal_cn: normal (diploid) copy number to anchor around (default: 2)
    - output_key: name of adata.obs column to store output
    '''

    cnv = adata.obsm[cnv_key]
    var = adata.var.loc[:, ["chromosome", "start", "end"]].reset_index(drop=True)

    # Convert 'start' and 'end' to numeric, handling NaNs
    var["start"] = pd.to_numeric(var["start"], errors="coerce").fillna(0).astype(int)
    var["end"] = pd.to_numeric(var["end"], errors="coerce").fillna(0).astype(int)
    var["chromosome"] = var["chromosome"].astype(str)

    var["chr_numeric"] = var["chromosome"].apply(
        lambda x: int(x) if x.isdigit() else (23 if x == "X" else 24 if x == "Y" else 25)
    )
    var_sorted = var.sort_values(by=["chr_numeric", "start"]).reset_index(drop=True)
    shared_indices = var_sorted.index.intersection(pd.RangeIndex(0, cnv.shape[1])) # find common indices
    cnv = cnv[:, shared_indices]  # Subset cnv based on shared indices
    var_sorted = var_sorted.iloc[shared_indices] # Subset var_sorted based on shared indices

    all_events = []
    for cell_idx in range(cnv.shape[0]):
        row = cnv[cell_idx].toarray()[0]  # Convert to dense and get the first row (since it's 1D)
        events = []

        start_idx = 0
        current_cn = round(normal_cn + row[0])
        current_chr = var_sorted.loc[0, "chromosome"]

        for i in range(1, len(row)):
            val = row[i]
            cn = round(normal_cn + val)
            chrom = var_sorted.loc[i, "chromosome"]

            # Breakpoint if CN changes or chromosome changes
            if abs(row[i] - row[i - 1]) > threshold or chrom != current_chr:
                prev_cn = round(normal_cn + row[i - 1])
                if abs(row[i - 1]) >= threshold and prev_cn != normal_cn:
                    start = var_sorted.loc[start_idx, "start"]
                    end = var_sorted.loc[i - 1, "end"]
                    event_str = f"{current_chr}:{start}-{end} (CN {prev_cn})"
                    events.append(event_str)
                start_idx = i
                current_chr = chrom
                current_cn = cn

        # Add final segment if it's a CNV
        final_cn = round(normal_cn + row[-1])
        if abs(row[-1]) >= threshold and final_cn != normal_cn:
            start = var_sorted.loc[start_idx, "start"]
            end = var_sorted.loc[len(row) - 1, "end"]
            chrom = var_sorted.loc[len(row) - 1, "chromosome"]
            event_str = f"{chrom}:{start}-{end} (CN {final_cn})"
            events.append(event_str)

        all_events.append(", ".join(events) if events else "")

    adata.obs[output_key] = all_events
    print(f"Added CNV annotations to adata.obs['{output_key}'].")

    return adata

# another option for segemention which is supposed to be better
def run_hmm_cnv_segmentation(
    adata,
    bin_size: int = 100,
    n_states: int = 3,
    use_layer: str = "smoothed",  # Fallback to X_cnv if unavailable
    obsm_key: str = "hmm_cnv_states",
    min_cnv_length: int = 3  # <--- New parameter for CNV filtering
):
    '''
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
    '''

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

