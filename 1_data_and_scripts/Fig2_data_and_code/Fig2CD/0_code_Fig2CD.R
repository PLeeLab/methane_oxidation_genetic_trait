# If all codons are used, the Nc value will be 61. 
# If only one codon is used for each amino acid the Nc value will be 20. 
# Low values therefore indicate a strong codon bias,
# and high values indicate a low bias (and possibly a non-coding region).
# http://www.bioinformatics.nl/cgi-bin/emboss/help/chips

#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("EnvStats")
#install.packages("ggpubr")
library('ggplot2')
library("ggpubr")
library("EnvStats")
library("ggpubr")


genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
Nc_all_df <- read.delim(file = "2_Nc_chips_CH4.tsv", header = FALSE)
colnames(Nc_all_df) <- c("CDS", "Nc")
dim(genes_metadata)[1] == dim(Nc_all_df)[1]
integrated_Nc_df <- merge(x = genes_metadata, y = Nc_all_df, by="CDS", all.x = TRUE)
length(Nc_all_df)==dim(integrated_Nc_df)[1]
integrated_Nc_df$CDS_unambiguous <- paste(integrated_Nc_df$CDS, integrated_Nc_df$curated_gene, sep = "_")
integrated_Nc_df <- integrated_Nc_df[!duplicated(integrated_Nc_df$CDS_unambiguous), ]
dim(Nc_all_df)[1]==dim(integrated_Nc_df)[1]
integrated_Nc_df$gene_category <- gsub(pattern = "pmo/amo", replacement = "pmo", x = integrated_Nc_df$gene_category)


my_comparisons <- list( c("pmo", "mmo"), c("pmo", "mxa"), c("pmo", "xox"), c("pmo", "hao"), c("pmo", "pxm") )
p_fig2c <- ggplot(integrated_Nc_df[integrated_Nc_df$group=="Ia",], aes(x=gene_category, y=Nc, fill=gene_category)) +
  geom_boxplot(outlier.alpha = 0)+
  geom_point(aes(fill=gene_category), shape=21, alpha = 0.5, position = position_jitter(width=0.2, height = 0)) +
  scale_x_discrete(limits=c("pmo", "mmo","mxa", "xox", "hao", "pxm")) +
  scale_fill_manual(values = c("red", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow")) +
  scale_colour_manual(values = c("red", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "yellow")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, method = "wilcox")+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        legend.position="none")+
  xlab(label = "Gene")+
  ylab(label = "ENC")+
  NULL
ggsave(filename = "Fig2C.pdf", plot = p_fig2c, device = "pdf", width = 5, height = 5, useDingbats=FALSE)


p_fig2d <- ggplot(integrated_Nc_df[integrated_Nc_df$group=="Ia" | integrated_Nc_df$group=="Ib", ], aes(x=GC3, y=Nc)) +
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=gene_category), alpha=0.8, size=2, color="black", shape=21)+ #geom_point(aes(fill=gene_category, shape=sample), alpha=0.5, size=2, color="black")+ 
  #scale_shape_manual(values=c(21, 24))+
  stat_smooth(method = "lm", se = FALSE, aes(colour=gene_category), show.legend = FALSE)+
  geom_rug(aes(colour=gene_category), show.legend = FALSE) +
  facet_grid(group ~ gene_category)+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        legend.position="none")+
  xlab(label = expression(GC[3]))+
  ylab(label = "ENC")+
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL
ggsave(filename = "Fig2D.pdf", plot = p_fig2d, device = "pdf", width = 11, height = 4, useDingbats=FALSE)


p_Fig2DS1 <- ggplot(integrated_Nc_df, aes(x=GC3, y=Nc)) +
  scale_fill_brewer(palette="Set1") +
  scale_colour_brewer(palette="Set1") +
  geom_point(aes(fill=gene_category), alpha=0.8, size=2, color="black", shape=21)+
  stat_smooth(method = "lm", se = FALSE, aes(colour=gene_category), show.legend = FALSE)+
  geom_rug(aes(colour=gene_category), show.legend = FALSE) +
  facet_grid(group ~ gene_category)+
  theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(face = "italic"),
        legend.text = element_text(face = "italic"),
        legend.position="none")+
  xlab(label = expression(GC[3]))+
  ylab(label = "ENC")+
  guides(fill = guide_legend(title="Gene", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black")))+
  NULL
ggsave(filename = "Fig2DS1.pdf", plot = p_Fig2DS1, device = "pdf", width = 9, height = 10, useDingbats=FALSE)
