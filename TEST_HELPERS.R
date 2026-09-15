# Test sintetici delle invarianti, non esecuzione del dataset reale.
source("00_config.R");source("01_helpers.R")
stopifnot(identical(keep_female_subjects(factor(c("3611","3553","1743","3608","3627","9999")),cfg),
  c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)))
stopifnot(identical(keep_female_subjects(c(3611,9999),cfg),c(FALSE,TRUE)))
stopifnot(identical(keep_female_subjects(c(" 3611 ","9999"),cfg),c(FALSE,TRUE)))
d<-data.frame(Age=c(30,44,35,39,52,46,38,29,47,41,33,58),
  EDSS=c(1,2,1,3,2,1,2,3,1,2,4,3),group=factor(rep(c("rrms","pira"),each=6),levels=c("rrms","pira")))
D<-make_D(d,TRUE,cfg)
stopifnot(all(c("Age","EDSS","grouppira")%in%colnames(D)),qr(D)$rank==ncol(D))
du<-make_D(d,FALSE,cfg);stopifnot(identical(colnames(du),c("(Intercept)","grouppira")))
dc<-d;dc$group<-factor(rep(c("ctrl","pira"),each=6));dc$EDSS<-NA_real_
stopifnot(!"EDSS"%in%colnames(make_D(dc,TRUE,cfg)))
dx<-d;dx$Age<-as.numeric(dx$group)
stopifnot(inherits(try(make_D(dx,TRUE,cfg),silent=TRUE),"try-error"))
# Aggregazione counts gene x cellula -> gene x soggetto conserva ordine e somme.
A<-Matrix::Matrix(matrix(1:24,nrow=4),sparse=TRUE)
ids<-c("B","A","B","C","A","C");target<-c("C","A","B")
M<-Matrix::sparseMatrix(i=seq_along(ids),j=match(ids,target),x=1,dims=c(6,3))
pb<-as.matrix(A%*%M)
expect<-sapply(target,function(id)Matrix::rowSums(A[,ids==id,drop=FALSE]))
stopifnot(isTRUE(all.equal(unname(pb),unname(expect))))
# ORA include tutti i set eleggibili, anche senza overlap.
cf<-cfg;cf$gene_set_min<-2;cf$min_ORA_genes<-2
o<-ora(c("a","b"),letters,list(hit=c("a","b"),nohit=c("x","y")),cf)
stopifnot(nrow(o)==2,o$PValue[o$pathway=="nohit"]==1)
cat("PASS: esclusioni, design, collinearita, aggregazione, universo ORA.\n")
