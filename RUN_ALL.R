# Aprire in RStudio ed eseguire Source: richiede scc e myeloid in memoria.
script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)),error=function(e)getwd())
if(!file.exists(file.path(script_dir,"00_config.R")))stop("Aprire RUN_ALL.R dalla cartella estratta della pipeline.")
for(f in c("00_config.R","01_helpers.R","02_prepare.R","03_analysis.R","04_RF.R"))source(file.path(script_dir,f))
if(!exists("scc")||!exists("myeloid"))stop("Caricare gli oggetti scc e myeloid prima di Source.")
set.seed(cfg$seed)
pira_data<-prepare_data(scc,myeloid,cfg)
saveRDS(cfg,file.path(cfg$output,"config.rds"))
contrasts<-list(
  list(name="PIRA_vs_RRMS",ref="rrms",target="pira"),
  list(name="PIRA_vs_ctrl_CSF",ref="ctrl",target="pira"),
  list(name="RRMS_vs_ctrl_CSF",ref="ctrl",target="rrms"),
  list(name="ctrl_CSF_vs_blood",ref="blood",target="CSF")
)
statuses<-list()
for(view in c("L1","L2","myeloid"))for(con in contrasts) {
  st<-tryCatch({run_contrast(pira_data,view,con,cfg);"completed"},error=function(e)paste("FAILED:",conditionMessage(e)))
  statuses[[paste(view,con$name)]]<-data.frame(view=view,contrast=con$name,status=st)
  message(view," / ",con$name," — ",st)
  write.csv(frame_bind(statuses),file.path(cfg$output,"MODULE_STATUS.csv"),row.names=FALSE)
}
capture.output(sessionInfo(),file=file.path(cfg$output,"sessionInfo.txt"))
method_files<-list.files(cfg$output,pattern="^methods[.]md$",recursive=TRUE,full.names=TRUE)
output_status<-frame_bind(lapply(method_files,function(f) {
  z<-readLines(f,warn=FALSE);status<-z[grepl("^Status:",z)]
  data.frame(module=dirname(f),status=paste(status,collapse="; "))
}))
write.csv(output_status,file.path(cfg$output,"MODULE_OUTPUT_STATUS.csv"),row.names=FALSE)
message("Core completato. Controllare MODULE_STATUS.csv: FAILED non significa risultato negativo.")
# Eseguire separatamente i moduli costosi; pira_data contiene SOLO donne.
# rf_result <- run_RF(pira_data,cfg,n_repeats=10,n_permutations=0,compute_shap=TRUE)
# rf_perm <- run_RF(pira_data,cfg,n_repeats=10,n_permutations=199,compute_shap=FALSE)
