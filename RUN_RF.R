# Eseguire DOPO RUN_ALL.R, nella stessa sessione R.
if(!exists("pira_data")||!exists("run_RF"))stop("Eseguire prima RUN_ALL.R.")
rf_result <- run_RF(pira_data,cfg,n_repeats=10,n_permutations=0,compute_shap=TRUE)
# Opzionale: molto costoso, ripete tutta la validazione con etichette permutate.
# rf_perm <- run_RF(pira_data,cfg,n_repeats=10,n_permutations=199,compute_shap=FALSE)
