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
table(df$group)
#show_col(hue_pal()(8))
manual_colors2 <- c("#CD9600", # Ia
                    "#7CAE00", # Ib
                    #"#00BE67", # Ic
                    "#00BFC4", # IIa *
                    "#00A9FF", # IIb
                    "#C77CFF")#, # III
                    #"#FF61CC") # NC10)


p_Fig1 <- ggplot(df[(df$gene_category=='pmo/amo' | df$gene_category=='mmo' | df$gene_category=='mxa' | df$gene_category=='xox'), ], aes(x=GC, y=GC3)) +
  geom_polygon(data=data.frame(x=c(0.35, 0.35, 0.65, 0.65), y=c(0.2, 0.35, 0.65, 0.2)), aes(x, y), fill="grey", alpha=0.3) +
  geom_point(aes(fill=group, pch=sample), alpha=0.7, size=2, colour="black") + 
  scale_fill_manual(values = manual_colors2) +
  scale_shape_manual(values=c(21, 24))+
  geom_abline(slope = 1, intercept = c(0, 0), lty=2, colour='black') +
  theme_bw() +
  facet_wrap(~gene_category_subunit) + 
  theme(strip.text = element_text(face = "italic"))+
  guides(fill = guide_legend(title="Methanotroph type", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black"))) +
  NULL
ggsave(filename = "Fig1.pdf", plot = p_Fig1, device = "pdf", width = 9, height = 7, useDingbats=FALSE)


#p_Fig1_simpler <- ggplot(df[(df$gene_category=='pmo/amo' | df$gene_category=='mmo' | df$gene_category=='mxa' | df$gene_category=='xox'), ], aes(x=GC, y=GC3)) +
#  geom_polygon(data=data.frame(x=c(0.35, 0.35, 0.65, 0.65), y=c(0.2, 0.35, 0.65, 0.2)), aes(x, y), fill="grey", alpha=0.3) +
#  geom_point(aes(colour=group, pch=sample), alpha=1/2, size=2) + 
#  scale_color_manual(values = manual_colors2) +
#  geom_abline(slope = 1, intercept = c(0, 0), lty=2, colour='black') +
#  theme_classic2() +
#  facet_wrap(~gene_category_subunit, scales = "free") +
#  theme(strip.background = element_blank(),
#        panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),
#        strip.text = element_text(face = "italic"))+
#  NULL
#ggsave(filename = "Fig1_simpler.pdf", plot = p_Fig1_simpler, device = "pdf", width = 9, height = 7, useDingbats=FALSE)


p_Fig1_controls <- ggplot(df, aes(x=GC, y=GC3)) +
  geom_polygon(data=data.frame(x=c(0.35, 0.35, 0.65, 0.65), y=c(0.2, 0.35, 0.65, 0.2)), aes(x, y), fill="grey", alpha=0.3) +
  geom_point(aes(fill=group, pch=sample), alpha=0.7, size=2, colour="black") + 
  scale_fill_manual(values = manual_colors2) +
  scale_shape_manual(values=c(21, 24))+
  geom_abline(slope = 1, intercept = c(0, 0), lty=2, colour='black') +
  theme_bw() +
  facet_wrap(~gene_category_subunit) + 
  theme(strip.text = element_text(face = "italic"))+
  guides(fill = guide_legend(title="Methanotroph type", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black"))) +
  NULL
ggsave(filename = "Fig1_controls.pdf", plot = p_Fig1_controls, device = "pdf", width = 9, height = 8, useDingbats=FALSE)


#p_Fig1CS2 <- ggplot(df[df$group=='Ia' & df$gene_category=='pmo/amo', ], aes(x=GC, y=GC3)) +
#  geom_polygon(data=data.frame(x=c(0.35, 0.35, 0.6, 0.6), y=c(0.25, 0.35, 0.6, 0.25)), aes(x, y), fill="grey", alpha=0.3) +
#  geom_point(aes(colour=group, pch=sample), alpha=0.3, size=2.5) + 
#  geom_abline(lty=2) + 
#  scale_color_manual(values = c("#CD9600")) +
#  theme_bw() +
#  facet_wrap(~gene_category_subunit) + 
#  theme()
#ggsave(filename = "Fig1CS2.pdf", plot = p_Fig1CS2, device = "pdf", width = 10, height = 3.5, useDingbats=FALSE)
