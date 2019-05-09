#install.packages('seqinr')
#install.packages("reshape")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("RColorBrewer")
library('seqinr')
library("reshape")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")

genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
ALL_query_CDSs <- read.fasta(file = "query_ALL_clean.fasta", seqtype = 'DNA')

# Amino acid usage
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
AA_melted <- melt(data = AA_integrated_df, measure.vars = 23:43)
colnames(AA_melted)[c(23, 24)] <- c("AA", "Fraction")
AA_melted <- AA_melted[AA_melted$AA != "Stp",]
AA_melted$nFraction <- AA_melted$Fraction*(AA_melted$Length/3)


brewer.pal(6, "Set1")   #to pick manual colors

p_AA_pmo_Ia <- ggplot(data = AA_melted[AA_melted$group=="Ia" & (AA_melted$gene_category == "pmo/amo"), ], aes(x = reorder(AA, Fraction, FUN = median), y = Fraction)) +
  geom_point(alpha = 0.5, position = position_jitter(width=0.1, height = 0), color="white", pch=21, fill="#984EA3") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  coord_flip()+
  ggtitle("PMMO")+
  xlab("Amino acid")+ylab("Fraction")+ylim(c(0,0.2))+theme_classic2()

p_AA_mmo_Ia <- ggplot(data = AA_melted[AA_melted$group=="Ia" & (AA_melted$gene_category == "mmo"), ], aes(x = reorder(AA, Fraction, FUN = median), y = Fraction)) +
  geom_point(alpha = 0.5, position = position_jitter(width=0.1, height = 0), color="white", pch=21, fill="#377EB8") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  coord_flip()+
  ggtitle("SMMO")+
  xlab("")+ylab("")+ylim(c(0,0.2))+theme_classic2()

p_AA_xox_Ia <- ggplot(data = AA_melted[AA_melted$group=="Ia" & (AA_melted$gene_category == "xox"), ], aes(x = reorder(AA, Fraction, FUN = median), y = Fraction)) +
  geom_point(alpha = 0.5, position = position_jitter(width=0.1, height = 0), color="white", pch=21, fill="#FF7F00") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  coord_flip()+
  ggtitle("XOX")+
  xlab("")+ylab("")+ylim(c(0,0.2))+theme_classic2()

p_AA_mxa_Ia <- ggplot(data = AA_melted[AA_melted$group=="Ia" & (AA_melted$gene_category == "mxa"), ], aes(x = reorder(AA, Fraction, FUN = median), y = Fraction)) +
  geom_point(alpha = 0.5, position = position_jitter(width=0.1, height = 0), color="white", pch=21, fill="#4DAF4A") +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  coord_flip()+
  ggtitle("MXA")+
  xlab("")+ylab("")+ylim(c(0,0.2))+theme_classic2()

ggsave(plot = ggarrange(p_AA_pmo_Ia, p_AA_mmo_Ia, p_AA_mxa_Ia, p_AA_xox_Ia, ncol = 4, nrow = 1), filename = "Fig4B.pdf", device = "pdf", width = 20, height = 5)
