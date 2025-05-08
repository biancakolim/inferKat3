
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

# removes cell cycle genes
def remove_cell_cycle(adata):
  adata = adata.copy()

  # Define S-phase genes
  s_genes = [
    'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1',
    'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1',
    'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76',
    'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51',
    'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM',
    'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'
  ]

# Define G2/M-phase genes
  g2m_genes = [
    'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A',
    'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
    'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB',
    'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
    'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5',
    'CDCA3', 'HN1', 'CENPA'
    ]

# Combine all cell cycle genes
  cc_genes = s_genes + g2m_genes

  adata = adata[:, ~adata.var_names.isin(cc_genes)].copy()
  print(f"Shape after removing cell cycle genes: {adata.shape}")
  return adata

def remove_hla(adata):
  hla_genes = [gene for gene in adata.var_names if gene.startswith('HLA-')]
  print(f"Found {len(hla_genes)} HLA genes to remove.")
  adata = adata[:, ~adata.var_names.isin(hla_genes)].copy()
  return adata


def filter_by_length(state_seq, min_len=3):
    '''Suppress CNV calls shorter than `min_len` consecutive bins.'''
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


def second_filter(adata, min_genes_per_chr):
  print("Filtering low-quality cells based on chromosomal coverage...")
  genes = adata.var
  expr = adata.layers['counts'].toarray() if 'counts' in adata.layers else adata.X.toarray()

  good_cells = []
  for i in range(expr.shape[0]):
      cell_expr = expr[i, :]
      nonzero_idx = cell_expr > 0
      cell_chr = genes['chromosome'][nonzero_idx]
      chr_rle = cell_chr.value_counts()
      if (chr_rle.min() >= min_genes_per_chr):
          good_cells.append(i)

  adata = adata[good_cells, :].copy()
  print(f"Remaining {adata.n_obs} cells after filtering")

  # Variance stabilization and mean-centering
  print("Applying sqrt(x)+sqrt(x+1) transform and mean-centering...")
  expr = adata.layers['counts'].toarray() if 'counts' in adata.layers else adata.X.toarray()
  expr = np.log1p(np.sqrt(expr) + np.sqrt(expr+1))
  expr = expr - expr.mean(axis=1, keepdims=True)
  adata.layers['counts'] = expr

  return adata
