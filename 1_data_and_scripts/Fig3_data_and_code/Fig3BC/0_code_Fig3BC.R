#install.packages('ggplot2')
#install.packages('ggrepel')
#install.packages('gridExtra')
#install.packages('scales')
#install.packages("EnvStats")
#install.packages("ggpubr")
#install.packages('ggridges')
library('ggplot2')
library('ggrepel')
library('gridExtra')
library('scales')
library("EnvStats")
library("ggpubr")
library('ggridges')

df <- read.delim(file = '1_QC_V_manuallycurated.txt', header = TRUE)
all_tAI_df <- read.delim(file = '2_ALL_CH4_tAI_dataframe.txt', header = TRUE)
groups_assignment <- read.csv(file="3_list_isolates.csv")
all_tAI_df_with_group <- merge(x = all_tAI_df, y = groups_assignment, by = "strain")

#------ New 3BC ------#
my_comparisonsB <- list( c("pmoC", "mxaF"), c("pmoC", "mxaI") )
p_Fig3B <- ggplot(df[(df$gene_category=='pmo/amo' | df$gene_category=='mxa' | df$gene_category_subunit=='xoxF') & df$group == "Ia", ],
                  aes(x=gene_category_subunit, y=ptAI, fill=gene_category, colour=gene_category)) +
  geom_point(size=2, alpha = 1/3, position = position_jitter(width=0.1, height = 0)) +
  geom_boxplot(fill="transparent", outlier.colour = "transparent")+
  #scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("#4DAF4A", "#984EA3", "orange")) +
  scale_colour_manual(values = c("#4DAF4A", "#984EA3", "orange")) +
  scale_x_discrete(limits=c("pmoC", "pmoB", "pmoA", "xoxF", "mxaF", "mxaI")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisonsB, method = "t.test", method.args = list(conf.level = 0.99, alternative = "two.sided") )+
  theme_classic2()+ 
  stat_n_text(y.pos=0, geom="text")+
  theme(axis.text=element_text(size=12),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(face="italic", colour = "black"),
        axis.title.x = element_blank(),
        axis.title=element_text(size=14, face="bold"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab(label = "Percentile rank of tAI")+
  NULL


my_comparisonsC <- list( c("pmoC", "pmoA"), c("pmoC", "hao") )
p_Fig3C <- ggplot(df[(df$gene_category=='pmo/amo' | df$gene_category=='pxm' | df$gene_category_subunit=='hao') & df$group=='Ia', ],
                  aes(x=gene_category_subunit, y=ptAI, fill=gene_category, colour=gene_category)) +
  geom_point(size=2, alpha = 1/3, position = position_jitter(width=0.1, height = 0)) +
  geom_boxplot(fill="transparent", outlier.colour = "transparent")+
  #scale_y_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = c("grey70", "#984EA3", "grey70")) +
  scale_colour_manual(values = c("grey70", "#984EA3", "grey70")) +
  scale_x_discrete(limits=c("pmoC", "pmoB", "pmoA", "hao", "pxmA", "pxmC", "pxmB")) +
  stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisonsC, method = "t.test", method.args = list(conf.level = 0.99, alternative = "two.sided") )+
  theme_classic2()+ 
  stat_n_text(y.pos=0, geom="text")+
  theme(axis.text=element_text(size=12),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(face="italic", colour = "black"),
        axis.title.x = element_blank(),
        axis.title=element_text(size=14, face="bold"),
        legend.position="none",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab(label = "")+
  NULL

ggsave(filename = "Fig3BC.pdf",
       plot = ggarrange(p_Fig3B, p_Fig3C, ncol = 2, nrow = 1),
       device = "pdf",
       width = 10,
       height = 5,
       useDingbats=FALSE)


p_Fig3BCS1 <- ggplot(df, aes(x=gene_category_subunit, y=ptAI)) +
  geom_point(alpha = 0.5, position = position_jitter(width=0.1, height = 0), color="grey") +
  facet_wrap(~group, ncol = 3, nrow = 3, scales = "free") +
  theme_bw() +
  stat_n_text(angle = 45, size = 2.5)+
  scale_x_discrete(limits=c("pmoC", "pmoB", "pmoA", "xoxF", "mxaF", "mxaI", "mmoB", "mmoC","mmoD","mmoX","mmoY","mmoZ", "hao", "pxmA", "pxmB", "pxmC")) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.7, color="black") +
  theme_classic2()+
  theme(axis.text.y  = element_text(size=12, colour = "black"),
        axis.text.x  = element_text(face="italic", angle = 55, hjust = 1, colour = "black"),
        axis.title.x = element_blank(),
        axis.title   = element_text(size=14, face="bold"),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab(label = "")+
  ylab(label = "Percentile rank of tAI")+
  NULL
ggsave(filename = "Fig3BCS1.pdf", plot = p_Fig3BCS1, device = "pdf", width = 10, height = 10, useDingbats=FALSE)
