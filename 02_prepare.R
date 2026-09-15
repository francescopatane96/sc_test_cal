prepare_data <- function(scc,myeloid,cfg) {
  need(c("Seurat","SeuratObject","Matrix","dplyr","ggplot2","msigdbr","UCell","edgeR","limma","statmod"))
  if(cfg$run_scProportionTest) need("scProportionTest")
  if(cfg$run_SCPA) need("SCPA")
  if(cfg$run_DESeq2) need("DESeq2")
  out <- cfg$output; dir.create(out,recursive=TRUE,showWarnings=FALSE)
  m <- scc[[]]
  required <- c(cfg$subject,cfg$condition,cfg$tissue,cfg$age,cfg$edss,cfg$azimuth_l1,cfg$azimuth_l2)
  if(!all(required %in% names(m))) stop("Metadati mancanti: ",paste(setdiff(required,names(m)),collapse=", "))
  m$subject <- trimws(as.character(m[[cfg$subject]]))
  m$group_raw <- tolower(trimws(as.character(m[[cfg$condition]])))
  m$tissue_raw <- tolower(trimws(as.character(m[[cfg$tissue]])))
  m$group <- unname(cfg$condition_map[m$group_raw]); m$tissue <- unname(cfg$tissue_map[m$tissue_raw])
  if(anyNA(m$subject) || any(!nzchar(m$subject))) stop("Subject mancante.")
  audit <- as.data.frame(table(subject=m$subject,Condition=m$group_raw,Tissue=m$tissue_raw));audit<-audit[audit$Freq>0,]
  audit$excluded_male <- audit$subject %in% cfg$male_ids
  write.csv(audit,file.path(out,"subject_exclusions.csv"),row.names=FALSE)
  writeLines(paste("Male IDs non presenti in scc:",paste(setdiff(cfg$male_ids,m$subject),collapse=", ")),file.path(out,"male_ID_audit.txt"))
  # Esclusione unica su scc prima di ogni analisi, anche per controlli/blood
  m <- m[keep_female_subjects(m$subject,cfg),,drop=FALSE]
  if(anyNA(m$group)||anyNA(m$tissue)) stop("Valori Condition/Tissue non riconosciuti. Aggiornare le mappe in config.")
  for(nm in c("Age","EDSS")) {
    src <- if(nm=="Age") cfg$age else cfg$edss
    z <- trimws(as.character(m[[src]])); z[z %in% c("","NA","N/A")]<-NA
    m[[nm]] <- suppressWarnings(as.numeric(gsub(",",".",z,fixed=TRUE)))
    if(any(!is.na(z)&!is.finite(m[[nm]]))) stop("Valori non numerici in ",nm)
  }
  # Controlli: modello indipendente, verifica assenza di ID condivisi CSF/blood
  ctr <- unique(m[m$group=="ctrl",c("subject","tissue")])
  if(anyDuplicated(ctr$subject)) stop("Controlli CSF/blood con ID condiviso: confermare pairing prima di usare modello indipendente.")
  # Un soggetto non deve avere piu condizioni; visite replicate non deducibili dai barcode
  for(id in unique(m$subject)) if(length(unique(m$group[m$subject==id]))!=1)
    stop("Soggetto con piu condizioni: ",id)
  m$unit <- paste(m$subject,m$tissue,sep="__")
  patient <- frame_bind(lapply(split(seq_len(nrow(m)),m$unit),function(ii) {
    for(nm in c("Age","EDSS")) if(length(unique(na.omit(m[[nm]][ii])))>1)
      stop("Valori clinici discordanti/visite ripetute: ",m$unit[ii[1]]," ",nm)
    get <- function(nm) {a<-unique(na.omit(m[[nm]][ii]));if(length(a))a[1] else NA_real_}
    data.frame(unit=m$unit[ii[1]],subject=m$subject[ii[1]],group=m$group[ii[1]],tissue=m$tissue[ii[1]],
      Age=get("Age"),EDSS=get("EDSS"),n_cells=length(ii))
  }))
  rownames(patient)<-patient$unit
  write.csv(patient,file.path(out,"patient_metadata_female.csv"),row.names=FALSE)
  m$L1<-as.character(m[[cfg$azimuth_l1]]);m$L2<-as.character(m[[cfg$azimuth_l2]])
  m$myeloid <- NA_character_
  if(is.null(myeloid)) stop("Caricare anche l'oggetto myeloid per le annotazioni manuali.")
  mm <- myeloid[[]]
  if(!cfg$myeloid_cluster %in% names(mm)) stop("myeloid_cluster assente in myeloid.")
  if(!cfg$subject %in% names(mm))stop("Colonna soggetto assente in myeloid.")
  retained_mm <- rownames(mm)[keep_female_subjects(mm[[cfg$subject]],cfg)]
  if(!all(retained_mm %in% rownames(m))) stop("Barcode myeloid non trovati in scc: correggere il mapping senza assegnazioni arbitrarie.")
  pos <- match(retained_mm,rownames(m));m$myeloid[pos]<-as.character(mm[retained_mm,cfg$myeloid_cluster])
  m$myeloid[pos[is.na(m$myeloid[pos])]] <- "unassigned"
  for(nm in c("L1","L2")) m[[nm]][is.na(m[[nm]])|!nzchar(m[[nm]])] <- "unassigned"
  # Score preesistenti: copia barcode-correlata, nessuna nuova inferenza su firme non documentate
  for(nm in grep("^score_",names(mm),value=TRUE)) {
    m[[nm]]<-NA_real_;m[[nm]][pos]<-as.numeric(mm[retained_mm,nm])
  }
  obj <- subset(scc,cells=rownames(m))
  if(inherits(obj[["RNA"]],"Assay5")) {
    obj<-SeuratObject::JoinLayers(obj,assay="RNA")
    counts<-SeuratObject::LayerData(obj,assay="RNA",layer="counts")
  } else counts<-SeuratObject::GetAssayData(obj,assay="RNA",slot="counts")
  counts<-counts[,rownames(m),drop=FALSE]
  if(!nrow(counts)||anyDuplicated(rownames(counts)))stop("Counts vuoti/ID genici duplicati.")
  v<-if(inherits(counts,"sparseMatrix"))summary(counts)$x else as.vector(counts)
  if(any(!is.finite(v))||any(v<0))stop("Counts negativi o non finiti.")
  if(cfg$counts_scale!="counts")stop("Per DE/pseudobulk servono counts autentici, non log-data.")
  fractional<-any(abs(v-round(v))>1e-8);rm(v)
  writeLines(if(fractional)"Counts frazionari: verificare provenienza; edgeR ammesso, DESeq2 non eseguito." else "Counts interi.",file.path(out,"counts_audit.txt"))
  lib<-Matrix::colSums(counts);if(any(lib<=0))stop("Cellule con counts totali zero.")
  logdata<-log1p(counts %*% Matrix::Diagonal(x=10000/lib));dimnames(logdata)<-dimnames(counts)
  # Normalizzazione calcolata per cellula da RNA counts, senza integrazione
  obj<-SeuratObject::SetAssayData(obj,assay="RNA",layer="data",new.data=logdata)
  obj$analysis_unit<-m$unit;obj$analysis_group<-m$group
  sets_file<-file.path(out,"gene_sets_reference.rds")
  if(file.exists(sets_file)) sets<-readRDS(sets_file) else {
    get_sets<-function(collection,sub=NULL) {
      args<-list(species="Homo sapiens")
      modern<-"collection" %in% names(formals(msigdbr::msigdbr))
      args[[if(modern)"collection" else "category"]]<-collection
      if(!is.null(sub))args[[if(modern)"subcollection" else "subcategory"]]<-sub
      do.call(msigdbr::msigdbr,args)
    }
    h<-get_sets("H");r<-get_sets("C2","CP:REACTOME")
    sets<-list(Hallmark=split(h$gene_symbol,h$gs_name),Reactome=split(r$gene_symbol,r$gs_name),reference=frame_bind(list(h,r)))
    saveRDS(sets,sets_file)
  }
  if(sum(rownames(counts)%in%unique(unlist(sets$Hallmark)))<100)stop("ID geni incompatibili con simboli HGNC Hallmark.")
  sig<-sets$Hallmark
  if(!is.null(cfg$custom_signatures_rds)) sig<-c(sig,readRDS(cfg$custom_signatures_rds))
  if(anyDuplicated(names(sig)))stop("Nomi firme duplicati.")
  covg<-data.frame(signature=names(sig),n_genes=lengths(sig),
    mapped=vapply(sig,function(z)length(intersect(z,rownames(counts))),integer(1)))
  write.csv(covg,file.path(out,"signature_coverage.csv"),row.names=FALSE)
  sig<-sig[covg$mapped>=cfg$gene_set_min]
  scores<-UCell::ScoreSignatures_UCell(counts,features=sig,maxRank=cfg$UCell_maxRank,name="",ncores=1)
  scores<-as.data.frame(scores)
  if(!all(rownames(m)%in%rownames(scores)))stop("UCell non allineato ai barcode.")
  scores<-scores[rownames(m),,drop=FALSE]
  saveRDS(list(scores=scores,meta=m),file.path(out,"UCell_cell_scores.rds"))
  # Audit e atlante: UMAP se disponibile, altrimenti figura di numerosita
  qc<-as.data.frame(table(m$group,m$tissue));names(qc)<-c("group","tissue","cells")
  p<-ggplot2::ggplot(qc,ggplot2::aes(group,cells,fill=tissue))+ggplot2::geom_col(position="dodge")+ggplot2::theme_bw()
  export_module(file.path(out,"00_cohort"),patient,p,
    paste("All analyses excluded subjects",paste(cfg$male_ids,collapse=", "),
      ". Samples were grouped by subject and tissue. PIRA had occurred before sampling. Controls were treated as independent.",
      "Subjects are assumed female unless in the supplied male list; this is not independently verified from sex metadata."))
  reds<-names(obj@reductions)
  if("umap"%in%reds) {
    um<-SeuratObject::Embeddings(obj,"umap")[rownames(m),1:2,drop=FALSE]
    for(view in c("L1","L2","myeloid")) {
      dd<-data.frame(UMAP1=um[,1],UMAP2=um[,2],label=m[[view]],tissue=m$tissue)
      dd<-dd[!is.na(dd$label),]
      p<-ggplot2::ggplot(dd,ggplot2::aes(UMAP1,UMAP2,color=label))+ggplot2::geom_point(size=.15,alpha=.5)+
        ggplot2::facet_wrap(~tissue)+ggplot2::theme_bw()
      export_module(file.path(out,paste0("atlas_",view)),dd,p,
        "Existing embeddings and annotations were retained after subject exclusion. Embeddings are descriptive; no outcome-driven reclustering was performed.")
    }
  }
  for(nm in grep("^score_",names(m),value=TRUE)) {
    z<-m[is.finite(m[[nm]]),,drop=FALSE]
    if(!nrow(z))next
    dd<-aggregate(z[[nm]],list(subject=z$subject,group=z$group,tissue=z$tissue),median)
    names(dd)[4]<-"score"
    p<-ggplot2::ggplot(dd,ggplot2::aes(group,score,color=group))+ggplot2::geom_jitter(width=.1,height=0)+ggplot2::facet_wrap(~tissue)+ggplot2::theme_bw()
    export_module(file.path(out,"existing_scores",nm),dd,p,"Previously calculated score; descriptive only. Gene list and scoring provenance require confirmation.")
  }
  list(meta=m,patient=patient,counts=counts,logdata=logdata,obj=obj,scores=scores,sets=sets,fractional=fractional)
}
