# Classificazione retrospettiva a livello SOGGETTO, CSF PIRA vs RRMS.
# Non stima rischio futuro o tempo alla PIRA. Non usa ctrl per addestrare C1.
run_RF <- function(dat,cfg,n_repeats=10,n_permutations=0,compute_shap=TRUE) {
  need(c("caret","ranger","pROC"));if(compute_shap)need("fastshap")
  out<-file.path(cfg$output,"RF_subject_CSF");dir.create(out,recursive=TRUE,showWarnings=FALSE)
  d<-dat$patient
  d<-d[d$tissue=="CSF"&d$group%in%c("pira","rrms")&complete.cases(d[,c("Age","EDSS")])&d$n_cells>=cfg$min_subject_cells,,drop=FALSE]
  if(any(table(factor(d$group,levels=c("pira","rrms")))<4))stop("RF: meno di 4 soggetti/gruppo dopo esclusione maschi e missing.")
  stopifnot(!any(d$subject%in%cfg$male_ids),!anyDuplicated(d$subject))
  m<-dat$meta[dat$meta$unit%in%d$unit,,drop=FALSE]
  M<-Matrix::sparseMatrix(i=seq_len(nrow(m)),j=match(m$unit,d$unit),x=1,dims=c(nrow(m),nrow(d)))
  pb<-as.matrix(dat$counts[,rownames(m),drop=FALSE]%*%M)
  genes<-t(log1p(sweep(pb,2,colSums(pb),"/")*1e6))
  colnames(genes)<-make.names(paste0("gene_",rownames(pb)),unique=TRUE)
  features_map<-data.frame(feature=colnames(genes),gene=rownames(pb));write.csv(features_map,file.path(out,"gene_names.csv"),row.names=FALSE)
  sig<-t(sapply(d$unit,function(id)apply(dat$scores[rownames(m)[m$unit==id],,drop=FALSE],2,median)))
  colnames(sig)<-make.names(paste0("signature_",colnames(dat$scores)),unique=TRUE)
  prop<-table(factor(m$unit,levels=d$unit),factor(m$L2))
  prop<-matrix(as.numeric(prop),nrow=nrow(prop),dimnames=dimnames(prop))
  prop<-prop/rowSums(prop)
  colnames(prop)<-make.names(paste0("proportion_",colnames(prop)),unique=TRUE)
  clinical<-d[,c("Age","EDSS"),drop=FALSE]
  model_names<-c("EDSS","Age","EDSS_Age","genes","signatures","composition","genes_covariates","combined")
  k<-min(5,min(table(d$group)))
  pfun<-function(object,newdata)as.numeric(predict(object,newdata=newdata,type="prob")$pira)
  run_once<-function(labels,details=FALSE,shap=FALSE) {
    preds<-imps<-selected<-list();aucs<-matrix(NA_real_,n_repeats,length(model_names),dimnames=list(NULL,model_names))
    for(r in seq_len(n_repeats)) {
      set.seed(cfg$seed+r);fold<-integer(nrow(d))
      for(g in c("pira","rrms")){ii<-which(labels==g);fold[ii]<-sample(rep(seq_len(k),length.out=length(ii)))}
      P<-matrix(NA_real_,nrow(d),length(model_names))
      for(f in seq_len(k)) {
        tr<-which(fold!=f);te<-which(fold==f)
        stopifnot(!any(d$subject[tr]%in%d$subject[te]))
        v<-apply(genes[tr,,drop=FALSE],2,var)
        good<-which(is.finite(v)&v>0&colMeans(genes[tr,,drop=FALSE]>log1p(1))>=.2)
        chosen<-head(good[order(v[good],decreasing=TRUE)],100)
        if(length(chosen)<2)stop("RF: pochi geni training.")
        allX<-list(clinical["EDSS"],clinical["Age"],clinical,as.data.frame(genes[,chosen,drop=FALSE]),
          as.data.frame(sig),as.data.frame(prop),cbind(clinical,genes[,chosen,drop=FALSE]),
          cbind(clinical,genes[,chosen,drop=FALSE],sig,prop))
        for(j in seq_along(model_names)) {
          A<-as.data.frame(allX[[j]][tr,,drop=FALSE]);B<-as.data.frame(allX[[j]][te,,drop=FALSE])
          varying<-vapply(A,function(x)length(unique(x))>1,logical(1));A<-A[,varying,drop=FALSE];B<-B[,varying,drop=FALSE]
          if(anyNA(A)||anyNA(B))stop("Feature RF mancanti: non codificare cluster assenti come attivita zero.")
          if(!ncol(A)){P[te,j]<-.5;next}
          set.seed(cfg$seed+1000*r+10*f+j)
          w<-1/as.numeric(table(labels[tr])[labels[tr]])
          fit<-caret::train(x=A,y=factor(labels[tr],levels=c("pira","rrms")),weights=w,
            method="ranger",trControl=caret::trainControl(method="none",classProbs=TRUE),
            tuneGrid=data.frame(mtry=max(1,floor(sqrt(ncol(A)))),splitrule="gini",min.node.size=3),
            num.trees=300,num.threads=2,importance="none")
          P[te,j]<-pfun(fit,B)
          if(details)selected[[paste(r,f,j)]]<-data.frame(repetition=r,fold=f,model=model_names[j],feature=names(A))
          if(shap && model_names[j]%in%c("genes","genes_covariates","combined")) {
            S<-fastshap::explain(fit,X=A,newdata=B,pred_wrapper=pfun,nsim=30,adjust=TRUE)
            imps[[paste(r,f,j)]]<-data.frame(repetition=r,fold=f,model=model_names[j],
              subject=rep(d$subject[te],each=ncol(S)),feature=rep(colnames(S),times=nrow(S)),
              SHAP=as.vector(t(S)),expression=as.vector(t(as.matrix(B))))
          }
        }
      }
      for(j in seq_along(model_names)) {
        aucs[r,j]<-as.numeric(pROC::auc(pROC::roc(labels,P[,j],levels=c("rrms","pira"),direction="<",quiet=TRUE)))
        if(details)preds[[paste(r,j)]]<-data.frame(repetition=r,subject=d$subject,group=labels,fold=fold,model=model_names[j],probability=P[,j])
      }
      message("RF repetition ",r,"/",n_repeats)
    }
    list(aucs=aucs,stat=colMeans(aucs),predictions=frame_bind(preds),SHAP=frame_bind(imps),selection=frame_bind(selected))
  }
  result<-run_once(d$group,TRUE,compute_shap)
  saveRDS(result,file.path(out,"RF_results.rds"))
  pp<-result$predictions
  metrics<-frame_bind(lapply(split(pp,list(pp$repetition,pp$model),drop=TRUE),function(z) {
    positive<-z$group=="pira";ord<-order(z$probability,decreasing=TRUE)
    # Average precision a blocchi di score uguali, baseline = prevalenza PIRA
    zz<-aggregate(cbind(tp=as.numeric(positive),n=rep(1,nrow(z))),list(score=z$probability),sum);zz<-zz[order(zz$score,decreasing=TRUE),]
    ap<-sum((zz$tp/sum(positive))*(cumsum(zz$tp)/cumsum(zz$n)))
    data.frame(repetition=z$repetition[1],model=z$model[1],
      AUC=as.numeric(pROC::auc(pROC::roc(z$group,z$probability,levels=c("rrms","pira"),direction="<",quiet=TRUE))),
      average_precision=ap,Brier=mean((z$probability-as.numeric(positive))^2),
      balanced_accuracy=mean(c(mean(z$probability[positive]>=.5),mean(z$probability[!positive]<.5))))
  }))
  p<-ggplot2::ggplot(metrics,ggplot2::aes(model,AUC,group=repetition))+ggplot2::geom_line(alpha=.3)+ggplot2::geom_point()+
    ggplot2::geom_hline(yintercept=.5,linetype="dashed")+ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1))
  export_module(out,metrics,p,paste("Retrospective PIRA/RRMS classification in female CSF subjects.",n_repeats,
    "repetitions of subject-stratified",k,"fold CV. Fixed ranger forest 300 trees, min.node.size=3.",
    "Gene selection (100 variable genes) within training only. Clinical features: age and EDSS; Hallmark medians and L2 proportions.",
    "Gene pseudobulk aggregates all CSF cells and therefore also reflects composition. Average precision and Brier are descriptive OOF metrics.",
    "Repeated splits are dependent, not independent replicates. Scores classify existing status, not future PIRA."))
  write.csv(pp,file.path(out,"OOF_predictions.csv"),row.names=FALSE)
  write.csv(result$selection,file.path(out,"feature_selection.csv"),row.names=FALSE)
  patient_scores<-frame_bind(lapply(split(pp,list(pp$model,pp$subject),drop=TRUE),function(z)
    data.frame(model=z$model[1],subject=z$subject[1],group=z$group[1],
      mean=mean(z$probability),min=min(z$probability),max=max(z$probability))))
  p<-ggplot2::ggplot(patient_scores,ggplot2::aes(mean,reorder(subject,mean),color=group))+
    ggplot2::geom_segment(ggplot2::aes(x=min,xend=max,yend=subject))+
    ggplot2::geom_point()+ggplot2::geom_vline(xintercept=.5,linetype="dashed")+
    ggplot2::facet_wrap(~model)+ggplot2::theme_bw()+ggplot2::labs(x="Score PIRA OOF: media e range tra ripetizioni",y="Soggetto")
  export_module(file.path(out,"patient_scores"),patient_scores,p,
    "Each dot is a subject's mean held-out score over repeated splits. Bars are min-max ranges, not confidence intervals. Threshold 0.5 is descriptive; scores have not been clinically calibrated.")
  deltas<-data.frame(repetition=seq_len(n_repeats),
    genes_added=result$aucs[,"genes_covariates"]-result$aucs[,"EDSS_Age"],
    all_added=result$aucs[,"combined"]-result$aucs[,"EDSS_Age"])
  dl<-data.frame(repetition=rep(deltas$repetition,2),model=rep(c("genes_added","all_added"),each=n_repeats),
    delta_AUC=c(deltas$genes_added,deltas$all_added))
  p<-ggplot2::ggplot(dl,ggplot2::aes(model,delta_AUC,group=repetition))+
    ggplot2::geom_line(alpha=.3)+ggplot2::geom_point()+ggplot2::geom_hline(yintercept=0,linetype="dashed")+ggplot2::theme_bw()
  export_module(file.path(out,"increment_over_clinical"),dl,p,
    "Paired AUC differences on identical subject splits relative to Age+EDSS. Repeated splits are dependent; no t-test across repetitions. This is descriptive incremental discrimination, not a formal conditional permutation test.")
  if(nrow(result$SHAP)) {
    ss<-result$SHAP
    # Se una feature non e selezionata, non le si assegna SHAP zero
    av<-aggregate(abs(ss$SHAP),ss[,c("model","subject","feature")],mean);names(av)[4]<-"mean_abs_SHAP"
    imp<-aggregate(mean_abs_SHAP~model+feature,av,mean)
    # Copertura necessaria per leggere la media condizionata alla selezione
    cov<-aggregate(subject~model+feature,av,length);names(cov)[3]<-"subjects_explained"
    imp<-merge(imp,cov,by=c("model","feature"))
    imp<-imp[order(-imp$mean_abs_SHAP),];top<-frame_bind(lapply(split(imp,imp$model),head,20))
    p<-ggplot2::ggplot(top,ggplot2::aes(mean_abs_SHAP,reorder(feature,mean_abs_SHAP)))+ggplot2::geom_col()+
      ggplot2::facet_wrap(~model,scales="free_y")+ggplot2::theme_bw()+ggplot2::labs(y=NULL)
    export_module(file.path(out,"SHAP"),imp,p,"SHAP on held-out subjects, training-only background, 30 Monte Carlo simulations. Means conditional on feature selection; subject coverage and feature-selection tables are provided. These are explanations, not causal effects.")
    saveRDS(ss,file.path(out,"SHAP","SHAP_long.rds"))
    keys<-paste(top$model,top$feature)
    sv<-ss[paste(ss$model,ss$feature)%in%keys,,drop=FALSE]
    sv$scaled_value<-ave(sv$expression,interaction(sv$model,sv$feature,drop=TRUE),FUN=function(x) {
      rr<-range(x);if(diff(rr)==0)rep(.5,length(x)) else (x-rr[1])/diff(rr)
    })
    # Plot subsample only; numerical importance uses every available explanation.
    set.seed(cfg$seed)
    if(nrow(sv)>20000)sv<-sv[sample.int(nrow(sv),20000),]
    p<-ggplot2::ggplot(sv,ggplot2::aes(SHAP,feature,color=scaled_value))+
      ggplot2::geom_jitter(width=0,height=.2,alpha=.5,size=.7)+
      ggplot2::geom_vline(xintercept=0,linetype="dashed")+
      ggplot2::scale_color_gradient(low="#3C6E9E",high="#C6613F")+
      ggplot2::facet_wrap(~model,scales="free_y")+ggplot2::theme_bw()+
      ggplot2::labs(x="SHAP: negativo verso RRMS | positivo verso PIRA",y=NULL,color="Valore relativo")
    export_module(file.path(out,"SHAP_direction"),sv,p,
      "Held-out subject explanations across repetitions; positive SHAP increases the model's PIRA probability. Colour is feature-wise min-max scaled expression or clinical value, used only for plotting. Points across repetitions are dependent.")
  }
  if(n_permutations>0) {
    null<-matrix(NA_real_,n_permutations,length(model_names),dimnames=list(NULL,model_names))
    for(b in seq_len(n_permutations)) {
      message("RF permutation ",b,"/",n_permutations)
      set.seed(cfg$seed+100000+b);lab<-sample(d$group)
      null[b,]<-run_once(lab)$stat
      saveRDS(null,file.path(out,"permutation_progress.rds"))
    }
    pval<-vapply(seq_along(model_names),function(j)(1+sum(null[,j]>=result$stat[j]))/(n_permutations+1),numeric(1))
    tab<-data.frame(model=model_names,observed=result$stat,p=pval,FDR=p.adjust(pval,"BH"))
    ll<-data.frame(model=rep(model_names,each=n_permutations),AUC=as.vector(null))
    p<-ggplot2::ggplot(ll,ggplot2::aes(AUC))+ggplot2::geom_histogram(bins=20)+ggplot2::geom_vline(data=tab,ggplot2::aes(xintercept=observed),color="red")+
      ggplot2::facet_wrap(~model)+ggplot2::theme_bw()
    export_module(file.path(out,"permutation"),tab,p,"Patient-label permutation, entire repeated CV rerun; (b+1)/(B+1). Global association test, not conditional evidence for genes independent of covariates.")
  }
  invisible(result)
}
