
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
