#install.packages('seqinr')
#install.packages("robustbase")
#install.packages("gplots")
#install.packages("reshape")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("RColorBrewer")
#install.packages("ggrepel")
library('seqinr')
library("robustbase")
library("gplots")
library("reshape")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("ggrepel")


genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
ALL_query_CDSs <- read.fasta(file = "query_ALL_clean.fasta", seqtype = 'DNA')

rscu_all <- sapply(ALL_query_CDSs, uco, index = "rscu")
rscu_all <- t(rscu_all)
rscu_all[is.na(rscu_all)] <- 0
rscu_all <- rscu_all[,-c(15, 49,51,57,59)] # removing Met, Trp and Stop codons

#----------------- Amino acid usage
AA_matrix <- matrix(nrow = length(ALL_query_CDSs), ncol = length(aaa()))
for (j in 1:length(ALL_query_CDSs)) {
  prot_temp <- AAstat(seqinr::translate(ALL_query_CDSs[[j]]), plot = FALSE)$Compo/ (length(ALL_query_CDSs[[j]])/3)
  AA_matrix[j,] <- prot_temp
}
colnames(AA_matrix) <- aaa()
rownames(AA_matrix) <- names(ALL_query_CDSs)
AA_df <- as.data.frame(AA_matrix)
AA_df$CDS <- rownames(AA_matrix)
AA_integrated_df <- merge(x = genes_metadata, y = AA_df, by="CDS", sort = FALSE)
AA_melted <- melt(data = AA_integrated_df, measure.vars = 25:45)
AA_melted <- AA_melted[, c(-24)]
colnames(AA_melted)[c(24, 25)] <- c("AA", "Fraction")
AA_melted <- AA_melted[AA_melted$AA != "Stp",]
AA_melted$nFraction <- AA_melted$Fraction*(AA_melted$Length/3)
#-----------------

df_rscu_all <- as.data.frame(rscu_all)
df_rscu_all$CDS <- rownames(df_rscu_all)
df_rscu_integrated_all <- merge(x = genes_metadata, y = df_rscu_all, by="CDS")

# New function
s2c_n_translate_fun <- function(x) {translate(s2c(x))}

colnames(df_rscu_integrated_all)[25:83] <- paste(aaa(sapply(X = colnames(df_rscu_integrated_all)[25:83], FUN = s2c_n_translate_fun)),
                                                 "_",
                                                 toupper(colnames(df_rscu_integrated_all)[25:83]),
                                                 sep = "")
# New function
sort_df = function (data, vars = names(data), decreasing = F){
  if (length(vars) == 0 || is.null(vars))
    return(data)
  data[do.call("order", list(what = data[, vars, drop = FALSE], decreasing = decreasing)), , drop = FALSE]
}

# New function
colSdApply <- function(x, ...)apply(X=x, MARGIN=2, FUN=sd, ...)


# linear model RSCU vs AA usage
linear_model_rscu_aafraction <- function(pick_group = "Ia", pick_gene = "pmo/amo", pick_colour = "blue"){
  sub_medians <- colMedians(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group==pick_group & (df_rscu_integrated_all$gene_category == pick_gene), c(25:83)]))
  sub_sds     <- colSdApply(x = as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group==pick_group & (df_rscu_integrated_all$gene_category == pick_gene), c(25:83)]))
  
  #sub_medians_ordered <- sub_medians[order(names(sub_medians))]
  pmo_rscu_df <- data.frame(aaa_codon = names(sub_medians), rscu_median = sub_medians, rscu_sd = sub_sds)
  pmo_rscu_df$aa <- gsub(pattern = "_...", replacement = "", x = pmo_rscu_df$aaa_codon)
  pmo_rscu_df$codon <- gsub(pattern = "..._", replacement = "", x = pmo_rscu_df$aaa_codon)
  rownames(pmo_rscu_df) <- NULL
  pmo_rscu_df <- sort_df(data = pmo_rscu_df, vars = "rscu_median", decreasing = TRUE)
  pmo_rscu_df <- pmo_rscu_df[!duplicated(pmo_rscu_df[,c('aa')]),]
  
  # aa fraction
  pmo_aa_fraction <- AA_integrated_df[AA_integrated_df$group == pick_group, ]
  pmo_aa_fraction <- pmo_aa_fraction[, c(9, 26:45)]
  pmo_aa_fraction <- pmo_aa_fraction[pmo_aa_fraction$gene_category == pick_gene, ]
  pmo_aa_fraction_vector <- as.numeric(x = colMedians(x = as.matrix(pmo_aa_fraction[,c(2:21)])))
  names(pmo_aa_fraction_vector) <- colnames(pmo_aa_fraction[c(2:21)])
  pmo_aa_fraction_df <- data.frame(aa = names(pmo_aa_fraction_vector), aa_fraction_median = pmo_aa_fraction_vector)
  rownames(pmo_aa_fraction_df) <- NULL
  pmo_aa_fraction_df$aa_fraction_sd <- colSdApply(x = pmo_aa_fraction[,c(2:21)])
  
  # rscu and aa fraction together with plot
  pmo_final_df <- merge(x = pmo_rscu_df, y = pmo_aa_fraction_df, by = "aa")
  prebiotic_amino_acids <- c("Ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Pro", "Ser", "Thr", "Val")# Longo, et al. (2013) PNAS, 110(6), 2135-2139.
  pmo_final_df$prebiotic <- ifelse(test = pmo_final_df$aa %in% prebiotic_amino_acids, yes = "prebiotic",no = "modern")
  
  plot_model <- ggplot(data = pmo_final_df, mapping = aes(x = rscu_median, y = aa_fraction_median, label=aa)) +
    geom_errorbar(aes(ymin=aa_fraction_median-aa_fraction_sd, ymax=aa_fraction_median+aa_fraction_sd), width=.1, colour="grey") +
    geom_errorbarh(aes(xmin=rscu_median-rscu_sd, xmax=rscu_median+rscu_sd), height = .002, colour="grey") +
    geom_point(aes(fill=prebiotic, shape=prebiotic), size=2) +
    scale_shape_manual(values=c(21, 24))+
    scale_fill_manual(values = c("white", "black")) +
    geom_smooth(method = "lm", colour = pick_colour, lwd=2) +
    geom_text_repel()+
    #xlim(c(0,max(pmo_final_df$rscu_median)))+
    ggtitle(label = "", subtitle = paste(pick_group, pick_gene, "| P =", signif(summary(lm(formula = aa_fraction_median ~ rscu_median, data = pmo_final_df))$coef[2,4], 2))) +
    xlab(label = "RSCU") +
    ylab(label = "Amino acid usage fraction") +
    theme_bw()+
    theme(strip.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    NULL
  
  return(plot_model)
}

ggsave(filename = "Fig3F.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "Ia", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "Ia", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "Ia", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        linear_model_rscu_aafraction(pick_group = "Ia", pick_gene = "xox", pick_colour = "#FF7F00"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)
ggsave(filename = "Fig3FS1.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "Ib", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "Ib", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "Ib", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        linear_model_rscu_aafraction(pick_group = "Ib", pick_gene = "xox", pick_colour = "#FF7F00"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)

ggsave(filename = "Fig3FS2.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "Ic", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "Ic", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "Ic", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        linear_model_rscu_aafraction(pick_group = "Ic", pick_gene = "xox", pick_colour = "#FF7F00"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)

ggsave(filename = "Fig3FS3.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "IIa", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "IIa", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "IIa", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)

ggsave(filename = "Fig3FS4.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "IIb", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "IIb", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "IIb", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)

ggsave(filename = "Fig3FS5.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "III", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "III", pick_gene = "xox", pick_colour = "#FF7F00"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)

ggsave(filename = "Fig3FS6.pdf",
       plot = ggarrange(linear_model_rscu_aafraction(pick_group = "NC10", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        linear_model_rscu_aafraction(pick_group = "NC10", pick_gene = "mmo", pick_colour = "#377EB8"),
                        linear_model_rscu_aafraction(pick_group = "NC10", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        ncol = 2,
                        nrow = 2),
       device = "pdf",
       width = 11,
       height = 9,
       useDingbats=FALSE)







# Prebiotuc AA enrichment
prebiotic_AA_enrichment <- function(pick_group = "Ia", pick_gene = "pmo/amo", pick_colour = "blue"){
  
  # aa fraction
  pmo_aa_fraction <- AA_integrated_df[AA_integrated_df$group == pick_group, ]
  pmo_aa_fraction <- pmo_aa_fraction[, c(9, 26:45)]
  pmo_aa_fraction <- pmo_aa_fraction[pmo_aa_fraction$gene_category == pick_gene, ]
  pmo_aa_fraction_vector <- as.numeric(x = colMedians(x = as.matrix(pmo_aa_fraction[,c(2:21)])))
  names(pmo_aa_fraction_vector) <- colnames(pmo_aa_fraction[c(2:21)])
  pmo_aa_fraction_df <- data.frame(aa = names(pmo_aa_fraction_vector), aa_fraction_median = pmo_aa_fraction_vector)
  rownames(pmo_aa_fraction_df) <- NULL
  pmo_aa_fraction_df$aa_fraction_sd <- colSdApply(x = pmo_aa_fraction[,c(2:21)])
  
  prebiotic_amino_acids <- c("Ala", "Asp", "Glu", "Gly", "Ile", "Leu", "Pro", "Ser", "Thr", "Val") # Longo, et al. (2013) PNAS, 110(6), 2135-2139.
  pmo_aa_fraction_df$prebiotic <- ifelse(test = pmo_aa_fraction_df$aa %in% prebiotic_amino_acids, yes = "prebiotic",no = "modern")
  
  plot_prebiotic_enrich <- ggplot(data = pmo_aa_fraction_df, mapping = aes(x = prebiotic, y = aa_fraction_median, label=aa)) +
    geom_boxplot(outlier.alpha = 0)+
    geom_point(aes(shape=prebiotic, fill=prebiotic), size=2) +
    stat_compare_means(comparisons = list(c("modern", "prebiotic")), method = "t.test")+
    scale_shape_manual(values=c(21, 24))+
    scale_fill_manual(values = c("white", "black")) +
    geom_text_repel()+
    ggtitle(label = "", subtitle = paste(pick_group, pick_gene)) +
    xlab(label = "Amino acid") +
    ylab(label = "Median amino acid usage fraction") +
    theme_classic2()+
    theme(strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "none")+
    NULL
  
  return(plot_prebiotic_enrich)
}

ggsave(filename = "Fig3G.pdf",
       plot = ggarrange(prebiotic_AA_enrichment(pick_group = "Ia", pick_gene = "pmo/amo", pick_colour = "#984EA3"),
                        prebiotic_AA_enrichment(pick_group = "Ia", pick_gene = "mmo", pick_colour = "#377EB8"),
                        prebiotic_AA_enrichment(pick_group = "Ia", pick_gene = "mxa", pick_colour = "#4DAF4A"),
                        prebiotic_AA_enrichment(pick_group = "Ia", pick_gene = "xox", pick_colour = "#FF7F00"),
                        ncol = 4,
                        nrow = 1),
       device = "pdf",
       width = 10,
       height = 5,
       useDingbats=FALSE)
