
# driver function for inferKat
def inferKat(adata, set_min_genes=500, low_dr = 0.02, high_dr = 0.98, min_genes_per_chr = 5, ks_cutoff = 0.2):

  # preprocess
  print("Step 1: read and filter data ...")
  ad_clean = cleanAd(adata, set_min_genes)
  ad_norm = normalizeAd(ad_clean)
  ad_filt, high_dr = filter_dr(ad_norm, low_dr, high_dr)

  # annotate genes and secondary filtering
  print("Step 2: annotate gene coordinates and secondary filtering genes...")
  cnv.io.genomic_position_from_biomart(ad_filt, biomart_gene_id="hgnc_symbol", species="hsapiens", inplace = True)
  ad_filt2 = remove_cell_cycle(ad_filt)
  ad_filt3 = remove_hla(ad_filt2)
  ad_filt4 = second_filter(ad_filt3, min_genes_per_chr)
  # sort by position
  ad_filt4.var = ad_filt4.var.sort_values(by=["chromosome", "start"])
  ad_filt4 = ad_filt4[:, ad_filt4.var.index]

  # smoothening
  print("Step 3: smoothing data with Gaussian filter...")
  expr = ad_filt4.layers['counts'] if 'counts' in ad_filt4.layers else ad_filt4.X.toarray()
  norm_mat_smooth = smooth_matrix(expr)
  ad_filt4.layers["smoothed"] = norm_mat_smooth # put back in adata

  print("Step 4: select and calculate baseline expression...")

  # Load smoothed expression (genes × cells)
  norm_mat_smooth = ad_filt4.layers["smoothed"].T  # shape: genes × cells
  gene_names = ad_filt4.var_names.tolist()
  cell_names = ad_filt4.obs_names.tolist()

  # Run synthetic baseline estimation
  expr_relat_df, syn_df, cl_labels = synthetic_baseline(norm_mat_smooth, gene_names, cell_names, min_cells=10)

  # Store synthetic clusters in AnnData
  ad_filt4.obs["synthetic_cluster"] = pd.Categorical([f"cluster_{i}" for i in cl_labels], ordered=True)

  # Optionally store relative expression and baseline
  ad_filt4.layers["relative"] = expr_relat_df.loc[:, ad_filt4.obs_names].values.T  # back to cells × genes
  ad_filt4.uns["synthetic_baseline"] = syn_df  # genes × clusters

  return ad_filt4
