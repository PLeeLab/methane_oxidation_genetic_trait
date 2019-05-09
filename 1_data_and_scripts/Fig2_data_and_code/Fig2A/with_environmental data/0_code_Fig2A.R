#install.packages('seqinr')
#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("forcats")
#install.packages("robustbase")
#install.packages("RColorBrewer")
library('seqinr')
library('ggplot2')
library("ggpubr")
library("forcats")
library("RColorBrewer")
library("robustbase")


genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
genes_metadata$gene_category <- gsub(pattern = "pmo/amo", replacement = "pmo", x = genes_metadata$gene_category)
ALL_query_CDSs <- read.fasta(file = "query_CH4_clean.fasta", seqtype = 'DNA')
dim(genes_metadata)[1] == length(ALL_query_CDSs)

rscu_all <- sapply(ALL_query_CDSs, uco, index = "rscu")
rscu_all <- t(rscu_all)
rscu_all[is.na(rscu_all)] <- 0
rscu_all <- rscu_all[,-c(15, 49,51,57,59)] # removing Met, Trp and Stop codons

rscu_Ia <- subset(x = rscu_all, subset = genes_metadata$group  =="Ia")
rscu_Ib <- subset(x = rscu_all, subset = genes_metadata$group  =="Ib")
rscu_IIa <- subset(x = rscu_all, subset = genes_metadata$group =="IIa")
rscu_IIb <- subset(x = rscu_all, subset = genes_metadata$group =="IIb")
rscu_III <- subset(x = rscu_all, subset = genes_metadata$group =="III")

#----------- Figure 2A -----------#
df_pca_rscu_Ia <- prcomp(rscu_Ia)
df_out_rscu_Ia <- as.data.frame(df_pca_rscu_Ia$x)
df_out_rscu_Ia$CDS <- rownames(df_out_rscu_Ia)
integrated_df_Ia <- merge(x = genes_metadata, y = df_out_rscu_Ia, by="CDS")
dim(rscu_Ia)[1]==dim(integrated_df_Ia)[1]
percentage_Ia <- round(df_pca_rscu_Ia$sdev / sum(df_pca_rscu_Ia$sdev) * 100, 2)
percentage_Ia <- paste( colnames(df_out_rscu_Ia), "(", paste( as.character(percentage_Ia), "%", ")", sep="") )



env_data <- read.csv(file = "7_environments_data.csv")
#env_data <- env_data[, c(3, 31)] # subset colums: "strain" and "Habitat_standard"

df_X <- merge(x = integrated_df_Ia, y = env_data, by="strain", sort = FALSE, all.x = TRUE)

ggplot(df_X, aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=gene_category), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL


ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=gene_category_subunit), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL


ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=Habitat_standard), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL


ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=gene_category_subunit), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL


ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=sample.x), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL



ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=Sequencing_Center), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL



ggplot(df_X[df_X$gene_category=="pmo", ], aes(x=PC1, y=PC2)) + 
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=Status), alpha=0.8, size=3, color="white", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.8, size=2, color="white")+ 
  #scale_shape_manual(values=c(21, 24))+
  xlab(percentage_Ia[1]) +
  ylab(percentage_Ia[2]) +
  #theme_classic2()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_blank(),
        legend.text = element_text(face = "italic"))+
  #ggtitle("RSCU Type Ia") +
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL

print(unique(df_X$strain[df_X$gene_category=="pmo" & df_X$PC1 > 0]), max.levels = 0)
print(unique(df_X$gene_category_subunit[df_X$gene_category=="pmo" & df_X$PC1 > 0]), max.levels = 0)
print(unique(df_X$Habitat_standard[df_X$gene_category=="pmo" & df_X$PC1 > 0]), max.levels = 0)




