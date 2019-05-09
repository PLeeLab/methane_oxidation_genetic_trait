#install.packages('seqinr')
#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("forcats")
#install.packages("robustbase")
#install.packages("vegan")
library('seqinr')
library("vegan")

genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
genes_metadata$gene_category <- gsub(pattern = "pmo/amo", replacement = "pmo", x = genes_metadata$gene_category)
ALL_query_CDSs <- read.fasta(file = "query_CH4_clean.fasta", seqtype = 'DNA')
dim(genes_metadata)[1] == length(ALL_query_CDSs)

rscu_all <- sapply(ALL_query_CDSs, uco, index = "rscu")
rscu_all <- t(rscu_all)
rscu_all[is.na(rscu_all)] <- 0
rscu_all <- rscu_all[,-c(15, 49,51,57,59)] # removing Met, Trp and the 3 Stop codons

#---- Type Ia ----#
rscu_Ia <- subset(x = rscu_all, subset = genes_metadata$group  =="Ia")

# PCA
perm_rscu_df_Ia <- data.frame(rscu_Ia)
perm_rscu_df_Ia$CDS <- rownames(rscu_Ia)
rownames(perm_rscu_df_Ia) <- NULL
perm_rscu_integrated_df_Ia <- merge(x = genes_metadata, y = perm_rscu_df_Ia, by="CDS")
#rscu_c <- scale(perm_rscu_integrated_df_Ia[ ,23:81])
rscu_c <- (perm_rscu_integrated_df_Ia[ ,23:81]) # select the columns of codons
pca <- rda(rscu_c)

# plot
plot(pca, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green', 'purple', 'grey', 'orange')
ordispider(pca, groups = perm_rscu_integrated_df_Ia$gene_category, label = TRUE, col = "grey")
points(pca, display='sites', bg = cols[as.factor(perm_rscu_integrated_df_Ia$gene_category)], pch = 21)
ordihull(pca, groups=perm_rscu_integrated_df_Ia$gene_category)

# PerMANOVA - partitioning the euclidean distance matrix by genes (DF = 5)
adonis(rscu_c ~ gene_category, data = perm_rscu_integrated_df_Ia, method='euclidean')
# The test is one of no effect of `gene_category`. We have found
# evidence against this hypothesis and thus could reject this hypothesis,
# instead accepting the alternative hypothesis that there is an effect of
# `gene_category`. As to the strength of this effect? ~37.162% of the sums of
# squares can be explained by `gene_category`.

# PerMANOVA - testing for pmo and others (DF = 1)
#perm_rscu_integrated_df_Ia$gene_perm_test <- ifelse(test = perm_rscu_integrated_df_Ia$gene_category == 'pmo', yes = 'pmo', no = 'other')
#adonis(formula = rscu_c ~ gene_perm_test, data = perm_rscu_integrated_df_Ia, method='euclidean')
# PerMANOVA - testing for pmo and xox separate  (DF = 2)
#perm_rscu_integrated_df_Ia$gene_perm_test <- ifelse(test = perm_rscu_integrated_df_Ia$gene_category == 'pmo', yes = 'pmo', no = ifelse(test = perm_rscu_integrated_df_Ia$gene_category == 'xox', yes = 'xox', no = 'other'))
#adonis(rscu_c ~ gene_perm_test, data = perm_rscu_integrated_df_Ia, method='euclidean')


#---- Type Ib ----#
rscu_Ib <- subset(x = rscu_all, subset = genes_metadata$group  =="Ib")
rownames(rscu_Ib)[c(15,17,19)] <- paste0(rownames(rscu_Ib)[c(15,17,19)], "_2")

# PCA
perm_rscu_df_Ib <- data.frame(rscu_Ib)
perm_rscu_df_Ib$CDS <- rownames(rscu_Ib)
rownames(perm_rscu_df_Ib) <- NULL
perm_rscu_integrated_df_Ib <- merge(x = genes_metadata, y = perm_rscu_df_Ib, by="CDS")
rscu_c_Ib <- (perm_rscu_integrated_df_Ib[ ,23:81])
pca_Ib <- rda(rscu_c_Ib)

# plot
plot(pca_Ib, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green', 'purple', 'grey', 'orange')
ordispider(pca_Ib, groups = perm_rscu_integrated_df_Ib$gene_category, label = TRUE, col = "grey")
points(pca_Ib, display='sites', bg = cols[as.factor(perm_rscu_integrated_df_Ib$gene_category)], pch = 21)
ordihull(pca_Ib, groups=perm_rscu_integrated_df_Ib$gene_category)

# PerMANOVA - partitioning the euclidean distance matrix by genes
adonis(rscu_c_Ib ~ gene_category, data = perm_rscu_integrated_df_Ib, method='euclidean')



#---- Type IIa ----#
rscu_IIa <- subset(x = rscu_all, subset = genes_metadata$group  =="IIa")

# PCA
perm_rscu_df_IIa <- data.frame(rscu_IIa)
perm_rscu_df_IIa$CDS <- rownames(rscu_IIa)
rownames(perm_rscu_df_IIa) <- NULL
perm_rscu_integrated_df_IIa <- merge(x = genes_metadata, y = perm_rscu_df_IIa, by="CDS")
rscu_c_IIa <- (perm_rscu_integrated_df_IIa[ ,23:81])
pca_IIa <- rda(rscu_c_IIa)

# plot
plot(pca_IIa, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green', 'purple', 'grey', 'orange')
ordispider(pca_IIa, groups = perm_rscu_integrated_df_IIa$gene_category, label = TRUE, col = "grey")
points(pca_IIa, display='sites', bg = cols[as.factor(perm_rscu_integrated_df_IIa$gene_category)], pch = 21)
ordihull(pca_IIa, groups=perm_rscu_integrated_df_IIa$gene_category)

adonis(rscu_c_IIa ~ gene_category, data = perm_rscu_integrated_df_IIa, method='euclidean')


#---- Type IIb ----#
rscu_IIb <- subset(x = rscu_all, subset = genes_metadata$group  =="IIb")

# PCA
perm_rscu_df_IIb <- data.frame(rscu_IIb)
perm_rscu_df_IIb$CDS <- rownames(rscu_IIb)
rownames(perm_rscu_df_IIb) <- NULL
perm_rscu_integrated_df_IIb <- merge(x = genes_metadata, y = perm_rscu_df_IIb, by="CDS")
rscu_c_IIb <- (perm_rscu_integrated_df_IIb[ ,23:81])
pca_IIb <- rda(rscu_c_IIb)

# plot
plot(pca_IIb, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green', 'purple', 'grey', 'orange')
ordispider(pca_IIb, groups = perm_rscu_integrated_df_IIb$gene_category, label = TRUE, col = "grey")
points(pca_IIb, display='sites', bg = cols[as.factor(perm_rscu_integrated_df_IIb$gene_category)], pch = 21)
ordihull(pca_IIb, groups=perm_rscu_integrated_df_IIb$gene_category)

adonis(rscu_c_IIb ~ gene_category, data = perm_rscu_integrated_df_IIb, method='euclidean')




#---- Type III ----#
rscu_III <- subset(x = rscu_all, subset = genes_metadata$group  =="III")

# PCA
perm_rscu_df_III <- data.frame(rscu_III)
perm_rscu_df_III$CDS <- rownames(rscu_III)
rownames(perm_rscu_df_III) <- NULL
perm_rscu_integrated_df_III <- merge(x = genes_metadata, y = perm_rscu_df_III, by="CDS")
rscu_c_III <- (perm_rscu_integrated_df_III[ ,23:81])
pca_III <- rda(rscu_c_III)

# plot
plot(pca_III, type = 'n', display = 'sites')
cols <- c('red', 'blue', 'green', 'purple', 'grey', 'orange')
ordispider(pca_III, groups = perm_rscu_integrated_df_III$gene_category, label = TRUE, col = "grey")
points(pca_III, display='sites', bg = cols[as.factor(perm_rscu_integrated_df_III$gene_category)], pch = 21)
ordihull(pca_III, groups=perm_rscu_integrated_df_III$gene_category)

adonis(rscu_c_III ~ gene_category, data = perm_rscu_integrated_df_III, method='euclidean')
