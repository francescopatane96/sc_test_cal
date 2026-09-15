run_contrast <- function(dat,view,contrast,cfg) {
  root<-file.path(cfg$output,view,contrast$name)
  dir.create(root,recursive=TRUE,showWarnings=FALSE)
  m<-dat$meta; pat<-dat$patient
  if(contrast$name=="ctrl_CSF_vs_blood") {
    pat<-pat[pat$group=="ctrl",,drop=FALSE];pat$group<-pat$tissue
  } else pat<-pat[pat$tissue=="CSF",,drop=FALSE]
  pat<-pat[pat$group%in%c(contrast$ref,contrast$target),,drop=FALSE]
  pat$group<-factor(pat$group,levels=c(contrast$ref,contrast$target))
  # Complete-case coorte condivisa per modello primario e non aggiustato
  req<-c("Age",if(contrast$name=="PIRA_vs_RRMS")"EDSS")
  pat$included<-complete.cases(pat[,req,drop=FALSE])
  write.csv(pat,file.path(root,"cohort_inclusion.csv"),row.names=FALSE)
  pat<-pat[pat$included,,drop=FALSE]
  m<-m[m$unit%in%pat$unit & !is.na(m[[view]]),,drop=FALSE]
  total<-table(m$unit);pat<-pat[pat$unit%in%names(total)[total>=cfg$min_subject_cells],,drop=FALSE]
  m<-m[m$unit%in%pat$unit,,drop=FALSE]
  if(any(table(pat$group)<cfg$min_subjects_group))stop("Pazienti insufficienti per ",view," ",contrast$name)
  pat<-pat[order(pat$unit),,drop=FALSE];rownames(pat)<-pat$unit
  D<-make_D(pat,TRUE,cfg)
  write.csv(D,file.path(root,"design.csv"))
  m$cluster<-m[[view]]; m$group<-pat$group[match(m$unit,pat$unit)]
  ct<-table(factor(m$cluster),factor(m$unit,levels=pat$unit))
  pr<-sweep(as.matrix(ct),2,colSums(ct),"/")
  method_base<-paste("Comparison:",contrast$target,"versus",contrast$ref,"; annotation:",view,
    "; subjects:",sum(pat$group==contrast$target),"target and",sum(pat$group==contrast$ref),"reference.",
    "Primary analyses used complete cases for age",if(contrast$name=="PIRA_vs_RRMS")"and baseline EDSS." else ".",
    "All five prespecified male subjects were excluded. Positive effects indicate the target group.",
    "BH families were defined within this contrast and annotation level; resolutions are not independent replication.")
  # Proporzioni per individuo e scProportionTest richiesto
  long<-as.data.frame(as.table(pr));names(long)<-c("cluster","unit","proportion")
  long$group<-pat$group[match(long$unit,pat$unit)]
  p<-ggplot2::ggplot(long,ggplot2::aes(group,proportion,color=group))+ggplot2::geom_boxplot(outlier.shape=NA)+
    ggplot2::geom_jitter(width=.1,height=0)+ggplot2::facet_wrap(~cluster,scales="free_y")+ggplot2::theme_bw()
  vary<-apply(pr,1,var)>0
  if(sum(vary)>=2) {
    f<-limma::eBayes(limma::lmFit(asin(sqrt(pr[vary,,drop=FALSE])),D),robust=TRUE)
    tt<-limma::topTable(f,coef=grep("^group",colnames(D)),number=Inf,sort.by="none");tt$cluster<-rownames(tt)
    tt$delta_percentage_points<-100*(rowMeans(pr[tt$cluster,pat$group==contrast$target,drop=FALSE])-rowMeans(pr[tt$cluster,pat$group==contrast$ref,drop=FALSE]))
    export_module(file.path(root,"proportions_subject"),tt,p,paste(method_base,
      "Subject-level asin-sqrt proportions were analyzed using moderated limma models. Denominator: all cells in this annotation view and specimen.",
      "The limma logFC column is an asin-sqrt effect, not a log2 fold change. Percentage-point differences are unadjusted descriptive estimates."))
  } else export_module(file.path(root,"proportions_subject"),long,p,method_base,"Fewer than two variable populations; no test.")
  write.csv(long,file.path(root,"proportions_subject","source_data.csv"),row.names=FALSE)
  if(cfg$run_scProportionTest) {
    tryCatch({
    o<-subset(dat$obj,cells=rownames(m));o$analysis_cluster<-m$cluster;o$analysis_group<-as.character(m$group)
    sc<-scProportionTest::sc_utils(o)
    sc<-scProportionTest::permutation_test(sc,cluster_identity="analysis_cluster",sample_identity="analysis_group",
      sample_1=contrast$ref,sample_2=contrast$target,n_permutations=cfg$scprop_permutations)
    tt<-as.data.frame(sc@results$permutation)
    psc<-scProportionTest::permutation_plot(sc)
    export_module(file.path(root,"scProportionTest"),tt,psc,paste(method_base,
      "Standard scProportionTest pooled cells; permutations/bootstrap:",cfg$scprop_permutations,
      ". Its cell-level resampling does not model independent subjects. Inference across subjects is reported separately; these p-values are exploratory."))
    },error=function(e)export_module(file.path(root,"scProportionTest"),data.frame(),NULL,
      method_base,paste("FAILED:",conditionMessage(e))))
  }
  de<-scde<-path<-scoretabs<-lopo<-list(); cache<-list();scorelong<-list();cluster_audit<-list()
  audit_cluster<-function(cl,d,status) {
    cluster_audit[[cl]]<<-data.frame(cluster=cl,n_reference=sum(d$group==contrast$ref),
      n_target=sum(d$group==contrast$target),status=status)
    write.csv(frame_bind(cluster_audit),file.path(root,"cluster_eligibility.csv"),row.names=FALSE)
  }
  for(cl in setdiff(rownames(ct),"unassigned")) {
    message(view," / ",contrast$name," / ",cl)
    ids<-colnames(ct)[ct[cl,]>=cfg$min_cluster_cells]
    d<-pat[match(ids,pat$unit),,drop=FALSE]
    if(any(table(d$group)<cfg$min_subjects_group)){audit_cluster(cl,d,"Insufficient subjects with enough cells");next}
    DD<-tryCatch(make_D(d,TRUE,cfg),error=function(e){audit_cluster(cl,d,conditionMessage(e));NULL});if(is.null(DD))next
    mi<-m[m$cluster==cl&m$unit%in%ids,,drop=FALSE]
    ix<-match(rownames(mi),colnames(dat$counts))
    M<-Matrix::sparseMatrix(i=seq_len(nrow(mi)),j=match(mi$unit,ids),x=1,dims=c(nrow(mi),length(ids)))
    pb<-as.matrix(dat$counts[,ix,drop=FALSE]%*%M);colnames(pb)<-ids
    ans<-tryCatch(fit_edge(pb,DD),error=function(e){audit_cluster(cl,d,conditionMessage(e));NULL});if(is.null(ans))next
    audit_cluster(cl,d,"Pseudobulk eligible")
    t<-ans$table;t$cluster<-cl;t$model<-"adjusted";de[[paste(cl,"adjusted")]]<-t
    DU<-make_D(d,FALSE,cfg)
    u<-fit_edge(pb,DU)$table;u$cluster<-cl;u$model<-"unadjusted";de[[paste(cl,"unadjusted")]]<-u
    cache[[cl]]<-list(pb=pb,d=d,D=DD,y=ans$y,mi=mi)
    saveRDS(list(counts=pb,metadata=d,design=DD),file.path(root,paste0("pb_",make.names(cl),".rds")))
    lc<-edgeR::cpm(ans$y,log=TRUE,prior.count=2)
    pc<-prcomp(t(lc),scale.=FALSE)
    dd<-data.frame(subject=d$subject,group=d$group,PC1=pc$x[,1],PC2=pc$x[,2])
    pp<-ggplot2::ggplot(dd,ggplot2::aes(PC1,PC2,color=group,label=subject))+ggplot2::geom_point()+
      ggplot2::geom_text(vjust=-.5,size=2.5)+ggplot2::theme_bw()
    export_module(file.path(root,"PCA",make.names(cl)),dd,pp,paste(method_base,"PCA on TMM logCPM; descriptive, no group-wise test on cells."))
    # Score UCell mediano per paziente/populazione; no smoothed scores
    smat<-sapply(ids,function(id)apply(dat$scores[rownames(mi)[mi$unit==id],,drop=FALSE],2,median))
    if(is.null(dim(smat)))smat<-matrix(smat,ncol=length(ids))
    rownames(smat)<-colnames(dat$scores);colnames(smat)<-ids
    vf<-apply(smat,1,var)>0
    if(sum(vf)>=2) {
      fs<-limma::eBayes(limma::lmFit(smat[vf,,drop=FALSE],DD),robust=TRUE)
      st<-limma::topTable(fs,coef=grep("^group",colnames(DD)),number=Inf,sort.by="none")
      st$signature<-rownames(st);st$cluster<-cl;st$PValue<-st$P.Value;scoretabs[[cl]]<-st
    }
    sl<-data.frame(signature=rep(rownames(smat),each=ncol(smat)),unit=rep(ids,nrow(smat)),score=as.vector(t(smat)))
    sl$cluster<-cl;sl$subject<-d$subject[match(sl$unit,d$unit)];sl$group<-d$group[match(sl$unit,d$unit)]
    sl$Age<-d$Age[match(sl$unit,d$unit)];sl$EDSS<-d$EDSS[match(sl$unit,d$unit)]
    scorelong[[cl]]<-sl
    for(coll in c("Hallmark","Reactome")) {
      idx<-limma::ids2indices(dat$sets[[coll]],rownames(ans$y));idx<-idx[lengths(idx)>=cfg$gene_set_min&lengths(idx)<=cfg$gene_set_max]
      if(length(idx)) {
        vv<-limma::voom(ans$y,DD,plot=FALSE)
        ca<-limma::camera(vv,idx,DD,contrast=grep("^group",colnames(DD)),inter.gene.cor=.01,sort=FALSE)
        ca$pathway<-rownames(ca);ca$cluster<-cl;ca$collection<-coll;path[[paste(cl,coll)]]<-ca
      }
    }
    if(cfg$run_single_cell_DE) {
      sel<-select_balanced(mi,cfg$single_cell_cap_subject,cfg$seed)
      oo<-subset(dat$obj,cells=sel);oo$analysis_group<-as.character(mi[sel,"group"])
      ss<-Seurat::FindMarkers(oo,group.by="analysis_group",ident.1=contrast$target,ident.2=contrast$ref,
        assay="RNA",test.use="wilcox",min.pct=.1,logfc.threshold=0)
      ss$gene<-rownames(ss);ss$cluster<-cl;scde[[cl]]<-ss
    }
    if(cfg$run_DESeq2 && !dat$fractional && all(pb==round(pb))) {
      ds<-DESeq2::DESeqDataSetFromMatrix(countData=pb,colData=d,design=formula(reformulate(c(
        if("Age"%in%colnames(DD))"Age",if("EDSS"%in%colnames(DD))"EDSS","group"))))
      ds<-ds[rowSums(DESeq2::counts(ds)>=10)>=cfg$min_subjects_group,]
      ds<-DESeq2::DESeq(ds,quiet=TRUE,sfType="poscounts")
      rr<-as.data.frame(DESeq2::results(ds,contrast=c("group",contrast$target,contrast$ref)))
      rr$gene<-rownames(rr);rr$cluster<-cl
      export_module(file.path(root,"DESeq2",make.names(cl)),rr,
        ggplot2::ggplot(rr,ggplot2::aes(log2FoldChange,-log10(padj)))+ggplot2::geom_point(alpha=.3)+ggplot2::theme_bw(),
        paste(method_base,"DESeq2 sensitivity analysis using authentic integer counts; no rounding or method selection by significance."))
    }
    if(cfg$run_LOPO) {
      target<-head(t$gene[order(t$PValue)],cfg$loo_top)
      for(id in d$unit) {
        dd<-d[d$unit!=id,,drop=FALSE]
        rr<-tryCatch({if(any(table(dd$group)<3))stop("Few subjects")
          fit_edge(pb[,dd$unit,drop=FALSE],make_D(dd,TRUE,cfg))$table},error=function(e)NULL)
        lopo[[paste(cl,id)]]<-data.frame(cluster=cl,omitted=id,gene=target,
          original_logFC=t$logFC[match(target,t$gene)],
          logFC=if(is.null(rr))rep(NA_real_,length(target))else rr$logFC[match(target,rr$gene)])
      }
    }
  }
  # Famiglie BH esplicite; nessuna selezione di modello in base alla significativita
  all<-frame_bind(de)
  if(nrow(all)) {
    all$FDR_global<-NA_real_
    for(mo in unique(all$model)){ii<-all$model==mo;all$FDR_global[ii]<-p.adjust(all$PValue[ii],"BH")}
    prim<-all[all$model=="adjusted",]
    export_module(file.path(root,"pseudobulk_DE"),all,volcano(prim,contrast$name,cfg),paste(method_base,
      "Counts summed per subject/tissue/population; TMM, filterByExpr, robust edgeR QL.",
      "Minimum cells:",cfg$min_cluster_cells,"; subjects/group:",cfg$min_subjects_group,
      ". BH across genes and clusters separately for each model. Adjusted and unadjusted models use the same subjects."))
    ort<-list();scpat<-list()
    for(cl in names(cache)) {
      aa<-prim[prim$cluster==cl,];ca<-cache[[cl]]
      for(direction in c("Up","Down")) {
        direction_ok<-if(direction=="Up")aa$logFC>=cfg$lfc else aa$logFC<=-cfg$lfc
        genes<-aa$gene[aa$FDR_global<cfg$fdr & direction_ok]
        for(coll in c("Hallmark","Reactome")) {
          ot<-ora(genes,aa$gene,dat$sets[[coll]],cfg)
          if(nrow(ot)){ot$cluster<-cl;ot$Direction<-direction;ot$collection<-coll;ort[[paste(cl,direction,coll)]]<-ot}
        }
      }
      if(cfg$run_SCPA && sum(aa$FDR_global<cfg$fdr & abs(aa$logFC)>=cfg$lfc)<cfg$min_ORA_genes) {
        # Analisi esplorativa pianificata quando DEG list insufficiente
        for(r in seq_len(cfg$scpa_repeats)) {
          matrices<-lapply(c(contrast$target,contrast$ref),function(g) {
            mm<-ca$mi[ca$mi$group==g,,drop=FALSE]
            cap<-max(1,floor(cfg$scpa_cap_condition/length(unique(mm$subject))))
            selected<-select_balanced(mm,cap,cfg$seed+r)
            as.matrix(dat$logdata[,selected,drop=FALSE])
          })
          scpa_sets<-lapply(names(dat$sets$Hallmark),function(nm)
            data.frame(Pathway=nm,Genes=unique(dat$sets$Hallmark[[nm]])))
          sp<-tryCatch(SCPA::compare_pathways(samples=matrices,pathways=scpa_sets,
            downsample=cfg$scpa_cap_condition,min_genes=cfg$gene_set_min,max_genes=cfg$gene_set_max),error=function(e){message(conditionMessage(e));NULL})
          if(!is.null(sp)) {
            st<-as.data.frame(sp);st$cluster<-cl;st$repeat_id<-r;scpat[[paste(cl,r)]]<-st
          }
        }
      }
    }
    oo<-frame_bind(ort)
    if(nrow(oo)) {oo$FDR_global<-NA_real_;for(coll in unique(oo$collection)){i<-oo$collection==coll;oo$FDR_global[i]<-p.adjust(oo$PValue[i],"BH")}}
    export_module(file.path(root,"ORA"),oo,rank_plot(oo),paste(method_base,
      "ORA on up/down edgeR genes with global BH<",cfg$fdr,"and absolute logFC>=",cfg$lfc,
      ". Background: all tested genes in the same population. Minimum mapped DEG list:",cfg$min_ORA_genes,
      ". Both directions and clusters corrected jointly within each collection. No eligible list is not evidence for absent biology."))
    sp<-frame_bind(scpat)
    ps<-NULL
    if(nrow(sp)) {
      qn<-names(sp)[tolower(names(sp))%in%c("qval","qvalue")]
      pn<-names(sp)[tolower(names(sp))%in%c("pathway","pathways")]
      if(length(qn)==1&&length(pn)==1) {
        sp$item<-paste(sp$cluster,sp[[pn]],sep=": ")
        top<-head(unique(sp$item[order(sp[[qn]],decreasing=TRUE)]),25);vv<-sp[sp$item%in%top,]
        ps<-ggplot2::ggplot(vv,ggplot2::aes(x=.data[[qn]],y=item,color=factor(repeat_id)))+ggplot2::geom_point()+ggplot2::theme_bw()+
          ggplot2::labs(x="SCPA Qval (higher = stronger perturbation; NOT Storey q)",y=NULL,color="Subsample")
      }
    }
    export_module(file.path(root,"SCPA"),sp,ps,paste(method_base,
      "Exploratory SCPA used only when fewer than",cfg$min_ORA_genes,"DEGs passed the fixed thresholds.",
      "Hallmark; balanced subject contributions; cap",cfg$scpa_cap_condition,"cells/condition; repeats",cfg$scpa_repeats,
      ". Qval is not a Storey q-value; native cell-level statistics do not establish patient-level significance. No patient-adjusted causal interpretation."),
      if(nrow(sp))"completed: exploratory cellular distribution analysis" else "Not run/no eligible result; see DEG availability and package messages.")
  } else export_module(file.path(root,"pseudobulk_DE"),data.frame(),NULL,method_base,"No eligible cluster/model.")
  pt<-frame_bind(path)
  if(nrow(pt)){pt$FDR_global<-NA_real_;for(coll in unique(pt$collection)){i<-pt$collection==coll;pt$FDR_global[i]<-p.adjust(pt$PValue[i],"BH")}}
  export_module(file.path(root,"CAMERA"),pt,rank_plot(pt),paste(method_base,
    "Competitive CAMERA tests on voom pseudobulk; inter-gene correlation=0.01; BH set x cluster within collection. Direction is relative enrichment, not functional activation."))
  st<-frame_bind(scoretabs);sl<-frame_bind(scorelong)
  if(nrow(st)) {
    st$FDR_global<-p.adjust(st$PValue,"BH");st$Direction<-ifelse(st$logFC>0,"Up","Down")
  }
  export_module(file.path(root,"UCell"),st,rank_plot(st,"signature"),paste(method_base,
    "UCell computed on original RNA ranks; maxRank",cfg$UCell_maxRank,
    ". Unsmoothed medians per subject/population, limma adjusted model. logFC is a score difference, not gene log2FC. BH across signatures x clusters."))
  if(nrow(sl))write.csv(sl,file.path(root,"UCell","source_data.csv"),row.names=FALSE)
  # Correlazioni within-group: evita correlazioni dovute soltanto alla differenza PIRA/ctrl
  cr<-list()
  if(nrow(sl))for(cl in unique(sl$cluster))for(g in unique(sl$group))for(sig in unique(sl$signature))for(cv in c("Age","EDSS")) {
    if(cv=="EDSS" && !(g%in%c("pira","rrms")))next
    z<-sl[sl$cluster==cl&sl$group==g&sl$signature==sig,,drop=FALSE]
    z<-z[is.finite(z[[cv]])&is.finite(z$score),]
    if(nrow(z)<6 || length(unique(z[[cv]]))<3 || length(unique(z$score))<3)next
    tt<-cor.test(z$score,z[[cv]],method="spearman",exact=FALSE)
    cr[[paste(cl,g,sig,cv)]]<-data.frame(cluster=cl,group=g,signature=sig,covariate=cv,n=nrow(z),rho=unname(tt$estimate),PValue=tt$p.value)
  }
  cr<-frame_bind(cr);cp<-NULL
  if(nrow(cr)) {
    cr$FDR_global<-p.adjust(cr$PValue,"BH");cr$item<-paste(cr$cluster,cr$group,cr$signature,cr$covariate,sep=": ")
    z<-head(cr[order(cr$FDR_global),],30)
    cp<-ggplot2::ggplot(z,ggplot2::aes(rho,reorder(item,rho),color=-log10(FDR_global)))+ggplot2::geom_point()+ggplot2::theme_bw()+ggplot2::labs(y=NULL)
  }
  export_module(file.path(root,"correlations"),cr,cp,paste(method_base,
    "Within-group Spearman correlations of subject-level signature medians with age/EDSS; at least 6 complete subjects and 3 unique values. Descriptive, unadjusted; BH across tests in this module."))
  sc<-frame_bind(scde);scp<-NULL
  if(nrow(sc)) {
    sc$FDR_global<-p.adjust(sc$p_val,"BH")
    colfc<-intersect(c("avg_log2FC","avg_logFC"),names(sc))[1]
    if(!is.na(colfc)){sc$logFC<-sc[[colfc]];scp<-volcano(sc,"Single-cell screening: not independent subjects",cfg)}
  }
  export_module(file.path(root,"single_cell_DE"),sc,scp,paste(method_base,
    "Seurat Wilcoxon on log-normalized RNA; min.pct=0.1, no FC prefilter, cap",cfg$single_cell_cap_subject,
    "cells/subject. Cell-level p-values are exploratory and are not independent patient-level evidence."))
  lo<-frame_bind(lopo);lp<-NULL
  if(nrow(lo)) {
    lo$item<-paste(lo$cluster,lo$gene,sep=": ")
    lp<-ggplot2::ggplot(lo,ggplot2::aes(logFC,item))+ggplot2::geom_point(alpha=.4)+ggplot2::geom_vline(xintercept=0,linetype="dashed")+
      ggplot2::theme_bw()+ggplot2::labs(y=NULL,x="LOPO log2FC; missing refits retained as NA")
  }
  export_module(file.path(root,"LOPO"),lo,lp,paste(method_base,
    "Top genes ranked by adjusted p-value were refitted after omission of each subject. This is influence analysis, not independent validation; failures remain NA."))
  invisible(root)
}
