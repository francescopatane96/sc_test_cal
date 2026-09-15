# Configurazione confermata dall'utente. Analisi principale solo donne.
cfg <- list(
  male_ids = c("3611", "3553", "1743", "3608", "3627"),
  subject = "sample", condition = "Condition", tissue = "Tissue",
  age = "Age", edss = "EDSS",
  azimuth_l1 = "predicted.celltype.l1", azimuth_l2 = "predicted.celltype.l2",
  myeloid_cluster = "myeloid_cluster",
  condition_map = c(pira="pira", rrms="rrms", rmss="rrms", ctrl="ctrl"),
  tissue_map = c(csf="CSF", blood="blood"),
  counts_scale = "counts", # counts RNA autentici/corretti; MAI dati gia normalizzati
  min_subject_cells = 30, min_cluster_cells = 10, min_subjects_group = 4,
  min_residual_df = 4, fdr = 0.05, lfc = 0.25,
  min_ORA_genes = 10, gene_set_min = 15, gene_set_max = 500,
  UCell_maxRank = 1500,
  run_scProportionTest = TRUE, scprop_permutations = 10000,
  run_single_cell_DE = TRUE, single_cell_cap_subject = 200,
  run_DESeq2 = FALSE, # Sensibilita opzionale; soltanto counts interi autentici
  run_SCPA = TRUE, scpa_cap_condition = 300, scpa_repeats = 3,
  run_LOPO = TRUE, loo_top = 20,
  seed = 123, output = "PIRA_results_female_v3",
  custom_signatures_rds = NULL # lista nominata di gene-set; documentare fonte/versione
)
