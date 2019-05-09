#install.packages('ggplot2')
#install.packages('ggrepel')
#install.packages('gridExtra')
#install.packages('scales')
#install.packages("ggpubr")
library('ggplot2')
library('ggrepel')
library('gridExtra')
library('scales')
library("ggpubr")


df <- read.delim(file = '1_QC_CH4.txt', header = TRUE)
#show_col(hue_pal()(8))
manual_colors2 <- c("#CD9600", # Ia
                    "#7CAE00", # Ib
                    "#00BFC4", # IIa *
                    "#00A9FF", # IIb
                    "#C77CFF") # III


p_Fig2B <- ggplot(df[(df$gene_category=='pmo/amo' | df$gene_category=='mxa' | df$gene_category=='xox'), ], aes(x=pCAI_genome, y=pCAI_HiExp)) + 
  geom_polygon(data=data.frame(x=c(0, 0,     100/4), y=c(0, 100,   100)),   aes(x, y), fill="grey", alpha=0.3) + # slope = 4-top
  #geom_polygon(data=data.frame(x=c(0, 100,   100),   y=c(0, 0,     100/4)), aes(x, y), fill="grey", alpha=0.3) + # slope = 4-bottom
  #geom_point(aes(colour=group), alpha=0.3, size=1.5) + #geom_point(aes(colour=group, pch=sample), alpha=0.3, size=1.5) + 
  geom_point(aes(fill=group), alpha=0.7, size=2, colour="black", pch=21) + 
  scale_fill_manual(values = manual_colors2) +
  geom_abline(slope = 1,   intercept = c(0, 0), lty=2, colour='black') +
  geom_abline(slope = 2,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/2, intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 3,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/3, intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 4,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/4, intercept = c(0, 0), lty=2, colour='grey') +
  facet_wrap(~gene_category_subunit, scales = "free") +
  theme_bw() +
  facet_wrap(~gene_category_subunit, scales = "free") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.title=element_blank(),
        legend.position = "bottom",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.minor = element_blank())+
  xlab(label = "Percentile rank of CAIgenome")+
  ylab(label = "Percentile rank of CAIribosome")+
  NULL
ggsave(filename = "Fig2B.pdf", plot = p_Fig2B, device = "pdf", width = 5.5, height = 5, useDingbats=FALSE)



p_Fig2BS1 <- ggplot(df, aes(x=pCAI_genome, y=pCAI_HiExp)) + 
  geom_polygon(data=data.frame(x=c(0, 0,     100/4), y=c(0, 100,   100)),   aes(x, y), fill="grey", alpha=0.3) + # slope = 4-top
  #geom_polygon(data=data.frame(x=c(0, 100,   100),   y=c(0, 0,     100/4)), aes(x, y), fill="grey", alpha=0.3) + # slope = 4-bottom
  geom_point(aes(fill=group), alpha=0.7, size=2, colour="black", pch=21) + 
  scale_fill_manual(values = manual_colors2) +
  geom_abline(slope = 1,   intercept = c(0, 0), lty=2, colour='black') +
  geom_abline(slope = 2,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/2, intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 3,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/3, intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 4,   intercept = c(0, 0), lty=2, colour='grey') +
  geom_abline(slope = 1/4, intercept = c(0, 0), lty=2, colour='grey') +
  facet_wrap(~gene_category_subunit, scales = "free") +
  theme_bw() +
  facet_wrap(~gene_category_subunit, scales = "free") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.title=element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        panel.grid.minor = element_blank())+
  xlab(label = "Percentile rank of CAIgenome")+
  ylab(label = "Percentile rank of CAIribosome")+
  NULL
ggsave(filename = "Fig2BS1.pdf", plot = p_Fig2BS1, device = "pdf", width = 9, height = 9, useDingbats=FALSE)
