#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("reshape")
#install.packages("ggpubr")
library("ggplot2")
library("dplyr")
library("reshape")
library("ggpubr")

rscu_trna_network <- read.delim(file = "RSCU_complete_network.txt")[,c(7,8,5,1,3)]
#rscu_trna_network <- rscu_trna_network[rscu_trna_network$gene != "genome" & rscu_trna_network$gene != "Ribo", ]
#ggplot(data = rscu_trna_network, mapping = aes(x = gene, y = tRNA_avail)) + stat_summary(fun.y = sum, geom = "bar", width = 0.7, color="black")

#read_by_trna_antiCodon
rscu_trna_network_unique_tRNA_anticodon <- group_by(.data = rscu_trna_network, gene)
rscu_trna_network_unique_tRNA_anticodon <- unique(rscu_trna_network_unique_tRNA_anticodon[c("gene", "read_by_trna_antiCodon", "tRNA_avail")])
relative_tRNA <- rscu_trna_network_unique_tRNA_anticodon %>%
  group_by(gene) %>%
  summarise(tRNA_avail = sum(tRNA_avail))
ggplot(data = relative_tRNA, mapping = aes(x = gene, y = tRNA_avail)) + stat_summary(fun.y = sum, geom = "bar", width = 0.7, color="black")


relative_tRNA <- data.frame(gene = relative_tRNA$gene)
absolute_tRNA <- data.frame(gene = relative_tRNA$gene)

j <- 2
for (i in seq(from=0, to = 2.0, by = 0.01)) {
  rscu_trna_network <- rscu_trna_network[rscu_trna_network$RSCU >= i, ]
  rscu_trna_network_unique_tRNA_anticodon <- group_by(.data = rscu_trna_network, gene)
  rscu_trna_network_unique_tRNA_anticodon <- unique(rscu_trna_network_unique_tRNA_anticodon[c("gene", "read_by_trna_antiCodon", "tRNA_avail")])
  temp_relative_tRNA <- rscu_trna_network_unique_tRNA_anticodon %>%
    group_by(gene) %>%
    summarise(total_tRNA = sum(tRNA_avail))
  relative_tRNA <- cbind(relative_tRNA, temp_relative_tRNA$total_tRNA/sum(temp_relative_tRNA$total_tRNA))
  colnames(relative_tRNA)[j] <- paste0("RSCU_", i)
  absolute_tRNA <- cbind(absolute_tRNA, temp_relative_tRNA$total_tRNA)
  colnames(absolute_tRNA)[j] <- paste0("RSCU_", i)
  j <- j + 1
}


relative_tRNA_df <- melt(data = relative_tRNA, id.vars = "gene")
colnames(relative_tRNA_df)[2:3] <- c("RSCU", "tRNA_avail")
relative_tRNA_df$RSCU_numeric <- as.numeric(gsub(pattern = "RSCU_", replacement = "", x = relative_tRNA_df$RSCU))

# Relative
p_Fig4E1 <- ggplot(data = relative_tRNA_df, mapping = aes(x = RSCU_numeric, colour = gene, y = tRNA_avail, group=gene))+
  geom_vline(xintercept = c(1.0))+
  #geom_line()+
  #geom_point()+
  geom_smooth(method = "loess")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.position = "none",
        axis.text.x = element_text(colour = "black", angle = 0),
        axis.text.y = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks=element_blank())+
  xlab(label = "RSCU cutoff")+
  ylab(label = "Fraction of tRNA pool")+
  #scale_y_continuous(breaks = seq(0.1, 0.25, 0.025), limits = c(0.1, 0.23))+
  #scale_colour_discrete(c("blue","#377EB8", "#4DAF4A", "#984EA3", "black", "#FF7F00"))+     
  #geom_text(aes(x = 1.5, y = 0.208, label = "pmo", color = "pmo"), fontface="italic") + 
  #geom_text(aes(x = 1.5, y = 0.220, label = "xox", color = "xox"), fontface="italic")+
  #geom_text(aes(x = 1.5, y = 0.179, label = "mxa", color = "mxa"), fontface="italic") +
  #geom_text(aes(x = 1.5, y = 0.160, label = "mmo", color = "mmo"), fontface="italic")+
  #geom_text(aes(x = 1.5, y = 0.125, label = "ribosomal", color = "Ribo"))+
  #geom_text(aes(x = 1.5, y = 0.112, label = "genome", color = "genome"))+
  NULL

# Absolute
absolute_tRNA_df <- melt(data = absolute_tRNA, id.vars = "gene")
colnames(absolute_tRNA_df)[2:3] <- c("RSCU", "tRNA_avail")
absolute_tRNA_df$RSCU_numeric <- as.numeric(gsub(pattern = "RSCU_", replacement = "", x = absolute_tRNA_df$RSCU))

p_Fig4E2 <- ggplot(data = absolute_tRNA_df, mapping = aes(x = RSCU_numeric, colour = gene, y = tRNA_avail, group=gene))+
  geom_vline(xintercept = c(1.0))+
  #geom_line()+
  #geom_point()+
  geom_smooth(method = "loess")+
  theme_bw()+theme(strip.background = element_blank(),
                   strip.text = element_text(face = "italic"),
                   legend.position = "top",
                   axis.text.x = element_text(colour = "black", angle = 0),
                   axis.text.y = element_text(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.ticks=element_blank())+
  xlab(label = "")+
  ylab(label = "tRNA isoacceptors")+
  guides(colour=guide_legend(nrow = 1)) +
  NULL


ggsave(filename = "Fig4E.pdf",
       plot   = ggarrange(p_Fig4E2, p_Fig4E1, ncol = 1, nrow = 2, align = "hv", common.legend = T),
       device = "pdf",
       width  = 5,
       height = 5,
       useDingbats=FALSE)

ggsave(filename = "Fig4E_square.pdf",
       plot   = ggarrange(p_Fig4E2, p_Fig4E1, ncol = 1, nrow = 2, align = "hv", common.legend = T),
       device = "pdf",
       width  = 3,
       height = 6,
       useDingbats=FALSE)

# Test on module
chisq.test(data_frame(expected = absolute_tRNA$RSCU_0[c(2:4,6)], observed = absolute_tRNA$RSCU_2[c(2:4,6)])) # module only
