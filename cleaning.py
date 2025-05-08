# Preprocessing and filtering step

# Limit noise, doublets, and low quality reads for more efficient downstream analysis
def cleanAd(ad, filter_mt=True, set_pct_counts_mt=20, filter_counts=True, set_min_genes=500, set_max_counts=30000, set_min_cells=3):
  ad.var['mt'] = ad.var_names.str.startswith('MT-')
  ribo_prefix = ("RPS","RPL")
  ad.var['ribo'] = ad.var_names.str.startswith(ribo_prefix)
  sc.pp.calculate_qc_metrics(ad, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=True)

  adClean = ad.copy()

  fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,4), gridspec_kw={'wspace':0.25}, constrained_layout=True)
  ax1_dict = sc.pl.scatter(adClean, x='total_counts', y='pct_counts_mt', ax=ax1, show=False)
  ax2_dict = sc.pl.scatter(adClean, x='total_counts', y='n_genes_by_counts',ax=ax2, show=False)
  ax3_dict = sc.pl.scatter(adClean, x='pct_counts_ribo', y='n_genes_by_counts',ax=ax3, show=False)
  fig.suptitle('Before filtering')
  plt.show()

  if filter_mt:
    adClean = adClean[adClean.obs['pct_counts_mt']<set_pct_counts_mt,:].copy()
  if filter_counts:
    sc.pp.filter_cells(adClean, min_genes=set_min_genes)
    sc.pp.filter_cells(adClean, max_counts=set_max_counts)
    sc.pp.filter_genes(adClean, min_cells=set_min_cells) # Remove outlier genes

  fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,4), gridspec_kw={'wspace':0.25}, constrained_layout=True)
  ax1_dict = sc.pl.scatter(adClean, x='total_counts', y='pct_counts_mt', ax=ax1, show=False)
  ax2_dict = sc.pl.scatter(adClean, x='total_counts', y='n_genes_by_counts',ax=ax2, show=False)
  ax3_dict = sc.pl.scatter(adClean, x='pct_counts_ribo', y='n_genes_by_counts',ax=ax3, show=False)
  fig.suptitle('After filtering')
  plt.show()

  print(adClean.shape)
  return adClean

  # Normalize to correct for cell-to-cell variation
def normalizeAd(ad):
  adNorm = ad.copy()
  adNorm.layers['counts'] = adNorm.X.copy()
  sc.pp.normalize_total(adNorm , target_sum=1e4)
  sc.pp.log1p(adNorm)
  return adNorm

# Detection rate
def filter_dr(adata, low_dr=0.02, high_dr=0.98):
  expression_matrix = adata.layers['counts'].toarray() if 'counts' in adata.layers else adata.X.toarray()
  detection_rate = np.sum(expression_matrix > 0, axis=0) / expression_matrix.shape[0]
  filtered_genes = detection_rate > low_dr
  ad_filtered = adata[:, filtered_genes].copy()
  print(f"{filtered_genes.sum()} genes passed low detection rate filtering")

  # quality check
  num_genes = ad_filtered.shape[1]
  if num_genes < 7000:
      print("WARNING: low data quality; assigned low_dr to high_dr...")
      high_dr = low_dr
      warning = "low data quality"
  else:
      warning = "data quality is ok"
      high_dr = high_dr
  print(warning)

  return ad_filtered, high_dr