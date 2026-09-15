# Eseguire una volta, in una sessione pulita se R richiede aggiornamenti.
install.packages(c("Seurat","dplyr","ggplot2","msigdbr","caret","ranger","pROC","fastshap","remotes","BiocManager"))
BiocManager::install(c("edgeR","limma","statmod","UCell"))
remotes::install_github("rpolicastro/scProportionTest")
remotes::install_github("jackbibby1/SCPA")
# DESeq2 opzionale: BiocManager::install("DESeq2")
