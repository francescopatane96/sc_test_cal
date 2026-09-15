need <- function(p) {
  x <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)]
  if(length(x)) stop("Pacchetti mancanti: ",paste(x,collapse=", "),". Vedere INSTALL.R.")
}
keep_female_subjects <- function(ids,cfg) {
  ids<-trimws(as.character(ids))
  if(anyNA(ids)||any(!nzchar(ids)))stop("Subject mancante.")
  !ids %in% cfg$male_ids
}
empty_plot <- function(text) ggplot2::ggplot() + ggplot2::theme_void() +
  ggplot2::annotate("text",x=0,y=0,label=text,size=4)
export_module <- function(path, tab, plot, methods, note="completed") {
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  write.csv(tab,file.path(path,"results.csv"),row.names=FALSE)
  if(is.null(plot)) plot <- empty_plot(note)
  ggplot2::ggsave(file.path(path,"figure.png"),plot,width=11,height=7,dpi=250,limitsize=FALSE)
  ggplot2::ggsave(file.path(path,"figure.pdf"),plot,width=11,height=7,limitsize=FALSE)
  writeLines(c("# Methods",methods,"",paste("Status:",note)),file.path(path,"methods.md"))
  invisible(tab)
}
frame_bind <- function(x) {
  if(!length(x)) return(data.frame())
  dplyr::bind_rows(x)
}
make_D <- function(d, adjust, cfg) {
  # Ricostruzione esplicita per evitare EDSS nei confronti con controlli
  terms <- c(if(adjust) "Age",if(adjust && all(d$group %in% c("pira","rrms"))) "EDSS","group")
  if(anyNA(d[,terms,drop=FALSE])) stop("Covariate mancanti.")
  const <- terms[vapply(d[,terms,drop=FALSE],function(x) length(unique(x))<2,logical(1))]
  if("group" %in% const) stop("Un solo gruppo.")
  D <- model.matrix(reformulate(setdiff(terms,const)),droplevels(d))
  if(qr(D)$rank<ncol(D)) stop("Disegno collineare/non identificabile.")
  if(nrow(D)-ncol(D)<cfg$min_residual_df) stop("Gradi di liberta insufficienti.")
  attr(D,"omitted_constants") <- const
  D
}
fit_edge <- function(pb,D) {
  if(any(colSums(pb)<=0)) stop("Libreria vuota.")
  y <- edgeR::DGEList(pb)
  keep <- edgeR::filterByExpr(y,design=D)
  if(sum(keep)<20) stop("Meno di 20 geni testabili.")
  y <- y[keep,,keep.lib.sizes=FALSE]
  if(any(y$samples$lib.size<=0)) stop("Libreria vuota dopo filtro.")
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y,D,robust=TRUE)
  f <- edgeR::glmQLFit(y,D,robust=TRUE)
  coef <- grep("^group",colnames(D))
  stopifnot(length(coef)==1)
  t <- edgeR::glmQLFTest(f,coef=coef)
  tab <- edgeR::topTags(t,n=Inf,sort.by="none")$table
  tab$gene <- rownames(tab)
  list(y=y,table=tab,coef=coef)
}
ora <- function(genes,universe,sets,cfg) {
  genes <- intersect(unique(genes),universe)
  if(length(genes)<cfg$min_ORA_genes) return(data.frame())
  rows <- lapply(names(sets),function(nm) {
    g <- intersect(unique(sets[[nm]]),universe)
    if(length(g)<cfg$gene_set_min || length(g)>cfg$gene_set_max) return(NULL)
    overlap <- intersect(g,genes)
    # Tutti i set eleggibili testati, incluso overlap nullo
    data.frame(pathway=nm,overlap=length(overlap),list_size=length(genes),
      set_size=length(g),universe_size=length(universe),
      fold_enrichment=(length(overlap)/length(genes))/(length(g)/length(universe)),
      PValue=phyper(length(overlap)-1,length(g),length(universe)-length(g),length(genes),lower.tail=FALSE),
      genes=paste(overlap,collapse=";"))
  })
  frame_bind(rows)
}
volcano <- function(t,title,cfg) {
  if(!nrow(t)) return(empty_plot("Nessun test disponibile"))
  t$highlight <- t$FDR_global<cfg$fdr & abs(t$logFC)>=cfg$lfc
  t$y <- -log10(pmax(t$FDR_global,.Machine$double.xmin))
  ggplot2::ggplot(t,ggplot2::aes(logFC,y,color=highlight)) +
    ggplot2::geom_point(size=.8,alpha=.5) + ggplot2::facet_wrap(~cluster,scales="free") +
    ggplot2::geom_hline(yintercept=-log10(cfg$fdr),linetype="dashed") +
    ggplot2::theme_bw() + ggplot2::labs(title=title,x="log2FC: gruppo target - riferimento",y="-log10 BH gene x cluster")
}
rank_plot <- function(t,label="pathway") {
  if(!nrow(t)) return(empty_plot("Nessun set eleggibile; non equivale ad assenza biologica"))
  v <- t[order(t$FDR_global),,drop=FALSE]; v <- head(v,30)
  v$item <- paste(v$cluster,v[[label]],sep=": ")
  v$item <- factor(v$item,levels=rev(unique(v$item)))
  v$score <- -log10(pmax(v$FDR_global,.Machine$double.xmin))
  ggplot2::ggplot(v,ggplot2::aes(score,item,color=Direction)) +
    ggplot2::geom_point(size=2) + ggplot2::theme_bw() +
    ggplot2::labs(x="-log10 BH; primi 30 risultati anche se non significativi",y=NULL)
}
select_balanced <- function(md,cap_subject,seed) {
  set.seed(seed)
  unlist(lapply(split(rownames(md),md$subject),function(z)
    z[sample.int(length(z),min(length(z),cap_subject))]),use.names=FALSE)
}
