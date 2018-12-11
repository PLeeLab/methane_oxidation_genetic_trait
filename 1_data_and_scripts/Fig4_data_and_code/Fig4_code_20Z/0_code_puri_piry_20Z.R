#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("EnvStats")
#source("https://bioconductor.org/biocLite.R"); biocLite("Biostrings")
#install.packages("seqinr")
#install.packages("reshape")
#install.packages("ggrepel")
#install.packages("dplyr")
#install.packages("scales")
#install.packages("gridExtra")
library('ggplot2')
library("ggpubr")
library("EnvStats")
library("Biostrings")
library("seqinr")
library("reshape")
library("ggrepel")
library("dplyr")
library("scales")
library("gridExtra")


fasta_seqinr <- read.fasta(file = "MALC.4-CDSs.fasta")
mage_code <- gsub(x = names(fasta_seqinr), pattern = "\\|.*", replacement = "")
mage_code <- gsub(x = mage_code, pattern = "MALC", replacement = "MEALZ")
rna_counts <- read.delim(file = "GSE51145_RPKM.tab")


get_A_comp <- function(x){table(x)[["a"]]} #a
get_G_comp <- function(x){table(x)[["g"]]} #g
get_T_comp <- function(x){table(x)[["t"]]} #t
get_C_comp <- function(x){table(x)[["c"]]} #c


mage_df <- data.frame(mage_code,
                      gc = sapply(X = fasta_seqinr, FUN = GC),
                      gc3 = sapply(X = fasta_seqinr, FUN = GC3),
                      a = sapply(X = fasta_seqinr, FUN = get_A_comp),
                      g = sapply(X = fasta_seqinr, FUN = get_G_comp),
                      t = sapply(X = fasta_seqinr, FUN = get_T_comp),
                      c = sapply(X = fasta_seqinr, FUN = get_C_comp),
                      length = sapply(X = fasta_seqinr, FUN = length))

integrated_df <- merge(x = mage_df, y = rna_counts, by.x="mage_code", by.y="locus.tag", sort=FALSE)
integrated_df_molten_basic <- melt(data = integrated_df, measure.vars = c(9:16))
colnames(integrated_df_molten_basic)[c(9,10)] <- c("experiment", "expression_level")
integrated_df_molten_basic$Purines <- (integrated_df_molten_basic$a * integrated_df_molten_basic$expression_level) + (integrated_df_molten_basic$g * integrated_df_molten_basic$expression_level)
integrated_df_molten_basic$Pyrimidines <- (integrated_df_molten_basic$t * integrated_df_molten_basic$expression_level) + (integrated_df_molten_basic$c * integrated_df_molten_basic$expression_level)

# Create a new column in the dataset to highligth the genes of the four operons "pmo", "smo","mxa", "xox" in the expression profile
integrated_df_molten_basic$operon <- ifelse(test = integrated_df_molten_basic$mage_code == "MEALZv4_0514" | integrated_df_molten_basic$mage_code == "MEALZv4_0515" | integrated_df_molten_basic$mage_code == "MEALZv4_0516", yes = "pmo",
                                            no = ifelse(test = integrated_df_molten_basic$mage_code == "MEALZv4_3497", yes = "xox",
                                                        no = ifelse(test = integrated_df_molten_basic$mage_code == "MEALZv4_3445" | integrated_df_molten_basic$mage_code == "MEALZv4_3448", yes = "mxa",
                                                                    no = "other")))
integrated_df_molten_basic$nPurines <-      (integrated_df_molten_basic$a + integrated_df_molten_basic$g) / integrated_df_molten_basic$length
integrated_df_molten_basic$nPyrimidines <-  (integrated_df_molten_basic$t + integrated_df_molten_basic$c) / integrated_df_molten_basic$length
integrated_df_molten_basic_length <- integrated_df_molten_basic[integrated_df_molten_basic$experiment == "br_1",]

p_puri_pyri_hist_genome_20Z <- ggplot(data = integrated_df_molten_basic_length)+
  geom_vline(xintercept = 0.5, lty=1, lwd=0.25)+
  geom_histogram(aes(x=nPurines),fill = "red", alpha = 0.25, bins = 50)+
  geom_histogram(aes(x=nPyrimidines),fill = "blue", alpha = 0.25, bins = 50)+
  theme_classic2()+
  theme(plot.margin = margin(5.5, 200, 5.5, 200, "pt"))+
  ggtitle(label = "a", subtitle = expression(paste(italic("Methylomicrobium alcaliphilum"), " 20Z"))) +
  xlab(label = "Content fraction in genome")+
  ylab(label = "Frequency")
t.test(integrated_df_molten_basic_length$nPurines, integrated_df_molten_basic_length$nPyrimidines, paired = FALSE)
wilcox.test(integrated_df_molten_basic_length$nPurines, integrated_df_molten_basic_length$nPyrimidines, paired = FALSE)

p_pyri_exp_20Z <- ggplot(data = integrated_df_molten_basic, mapping = aes(x = nPyrimidines, y = (expression_level+0.5))) +
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic$operon == "other", yes = 2, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic$operon == "other", yes = 0, no = 2),
             fill   = ifelse(test = integrated_df_molten_basic$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_vline(xintercept = 0.5, lty=2) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(sides = "l")+
  facet_wrap(~experiment, nrow = 2, ncol = 4, scales = "free")+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = "b")+
  xlab(label = "Pyrimidine fraction")+
  ylab(label = "mRNA abundance")


integrated_df_molten <- melt(data = integrated_df, measure.vars = c(4:7))
colnames(integrated_df_molten)[c(13,14)] <- c("base", "count_genome")
integrated_df_molten <- melt(data = integrated_df_molten, measure.vars = c(5:12))
colnames(integrated_df_molten)[c(7,8)] <- c("experiment", "expression_level")
integrated_df_molten <- integrated_df_molten[,c(1,5:8)]
integrated_df_molten$count_transcriptome <- integrated_df_molten$count_genome*integrated_df_molten$expression_level
integrated_df_molten$count_transcriptome_plus_genome <- integrated_df_molten$count_transcriptome + integrated_df_molten$count_genome # NEW
integrated_df_molten <- melt(data = integrated_df_molten, measure.vars = c(3,6,7))
colnames(integrated_df_molten)[c(5,6)] <- c("omic", "count")
integrated_df_molten$compound <- ifelse(test = integrated_df_molten$base == "a" | integrated_df_molten$base == "g", yes = "Purine", no = "Pyrimidine")

# removing mxa
integrated_df_molten_no_mxa <- integrated_df_molten[integrated_df_molten$mage_code != "MEALZv4_3445" & integrated_df_molten$mage_code != "MEALZv4_3448", ]

# removing xox
integrated_df_molten_no_xox <- integrated_df_molten[integrated_df_molten$mage_code != "MEALZv4_3497", ]

# removing pmo
integrated_df_molten_no_pmo <- integrated_df_molten[integrated_df_molten$mage_code != "MEALZv4_0514" & integrated_df_molten$mage_code != "MEALZv4_0515" & integrated_df_molten$mage_code != "MEALZv4_0516", ]

t_casted <- cast(data = integrated_df_molten[integrated_df_molten$omic=="count_transcriptome", ], formula = base~experiment, value = "count", fun.aggregate = sum )
for (u in 2:9) {t_casted[,u] <- t_casted[,u] / sum(t_casted[,u]) }; colSums(t_casted)
t_molten <- melt(data = t_casted, measure.vars = c(2:9))
t_molten$pmo <- "1T"

gtm_casted <- cast(data = integrated_df_molten_no_mxa[integrated_df_molten_no_mxa$omic=="count_transcriptome", ], formula = base~experiment, value = "count", fun.aggregate = sum )
for (u in 2:9) {gtm_casted[,u] <- gtm_casted[,u] / sum(gtm_casted[,u]) }; colSums(gtm_casted)
gtm_molten <- melt(data = gtm_casted, measure.vars = c(2:9))
gtm_molten$pmo <- "2T-mxa"

gtx_casted <- cast(data = integrated_df_molten_no_xox[integrated_df_molten_no_xox$omic=="count_transcriptome", ], formula = base~experiment, value = "count", fun.aggregate = sum )
for (u in 2:9) {gtx_casted[,u] <- gtx_casted[,u] / sum(gtx_casted[,u]) }; colSums(gtx_casted)
gtx_molten <- melt(data = gtx_casted, measure.vars = c(2:9))
gtx_molten$pmo <- "3T-xox"

gtp_casted <- cast(data = integrated_df_molten_no_pmo[integrated_df_molten_no_pmo$omic=="count_transcriptome", ], formula = base~experiment, value = "count", fun.aggregate = sum )
for (u in 2:9) {gtp_casted[,u] <- gtp_casted[,u] / sum(gtp_casted[,u]) }; colSums(gtp_casted)
gtp_molten <- melt(data = gtp_casted, measure.vars = c(2:9))
gtp_molten$pmo <- "4T-pmo"

z_dataframe <- rbind(t_molten, gtm_molten, gtx_molten, gtp_molten)
rownames(z_dataframe) <- NULL
z_dataframe$gc <- ifelse(test = z_dataframe$base == "g" | z_dataframe$base == "c", yes = "GC", no = "AT")
z_dataframe$compound <- ifelse(test = z_dataframe$base == "a" | z_dataframe$base == "g", yes = "Purine", no = "Pyrimidine")

p_base_20Z <- ggplot(data = z_dataframe, aes(x = pmo, y = value*100, colour=base, group=base))+
  geom_point()+
  geom_line()+
  geom_abline(slope = 0, intercept = 25, lty=2)+
  facet_wrap(~experiment, nrow = 2, ncol = 4, scales = "free")+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = "d")+
  ylab("Percentage")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")+
  guides(colour=guide_legend(title=""))

p_gc_20Z <- ggplot(data = z_dataframe, aes(x = pmo, y = value*100, colour=gc, group=gc))+
  scale_color_brewer(palette = "Set2")+
  stat_summary(fun.y = sum, fun.ymin = sum, fun.ymax = sum, geom = "point")+
  stat_summary(fun.y=sum, geom="line")+
  geom_abline(slope = 0, intercept = 50, lty=2)+
  facet_wrap(~experiment, nrow = 2, ncol = 4, scales = "free")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = "e")+
  ylab("Percentage")+
  xlab("")+
  guides(colour=guide_legend(title=""))

p_compound_20Z <- ggplot(data = z_dataframe, aes(x = pmo, y = value*100, colour=compound, group=compound))+
  scale_color_brewer(palette = "Set1")+
  stat_summary(fun.y = sum, fun.ymin = sum, fun.ymax = sum, geom = "point")+
  stat_summary(fun.y=sum, geom="line")+
  geom_abline(slope = 0, intercept = 50, lty=2)+
  facet_wrap(~experiment, nrow = 2, ncol = 4, scales = "free")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = "f")+
  ylab("Percentage")+
  xlab("")+
  guides(colour=guide_legend(title=""))

integrated_df_molten_basic$Carbon   <- ((integrated_df_molten_basic$a * 5) + (integrated_df_molten_basic$g * 5) + (integrated_df_molten_basic$c * 4) + (integrated_df_molten_basic$t * 4))/integrated_df_molten_basic$length   # U has 4 carbons and T has 5
integrated_df_molten_basic$Hydrogen <- ((integrated_df_molten_basic$a * 5) + (integrated_df_molten_basic$g * 5) + (integrated_df_molten_basic$c * 5) + (integrated_df_molten_basic$t * 4))/integrated_df_molten_basic$length  # U has 4 Hydrogens and T has 6
integrated_df_molten_basic$Nitrogen <- ((integrated_df_molten_basic$a * 5) + (integrated_df_molten_basic$g * 5) + (integrated_df_molten_basic$c * 3) + (integrated_df_molten_basic$t * 2))/integrated_df_molten_basic$length  # U and T have the same Nitrogen content = 2
integrated_df_molten_basic$Oxygen   <- ((integrated_df_molten_basic$a * 0) + (integrated_df_molten_basic$g * 1) + (integrated_df_molten_basic$c * 1) + (integrated_df_molten_basic$t * 2))/integrated_df_molten_basic$length  # U and T have the same Oxygen content = 2
integrated_df_molten_basic_atoms <- integrated_df_molten_basic[integrated_df_molten_basic$experiment == "br_1",] # only 1 rep is required to calculate atom content as mRNA abundance won't be accounted for

p_atoms_NC_20Z<- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = Nitrogen*3, y = Carbon*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 2, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_hline(yintercept = mean(integrated_df_molten_basic_atoms$Carbon)*3, lty=2) +
  geom_vline(xintercept = mean(integrated_df_molten_basic_atoms$Nitrogen)*3, lty=2) +
  geom_rug(alpha=1/20) +
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 200, "pt"))+
  ggtitle(label = "g")+
  xlab(label = "Nitrogen atoms per codon")+
  ylab(label = "Carbon atoms per codon")+
  NULL

p_atoms_OH_20Z <- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = Oxygen*3, y = Hydrogen*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 2, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_vline(xintercept = mean(integrated_df_molten_basic_atoms$Oxygen)*3, lty=2) + # x-axis mean
  geom_hline(yintercept = mean(integrated_df_molten_basic_atoms$Hydrogen)*3, lty=2) + # y-axis mean
  geom_rug(alpha=1/20) +
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 200, "pt"))+
  ggtitle(label = "h")+
  xlab(label = "Oxygen atoms per codon")+
  ylab(label = "Hydrogen atoms per codon")+
  NULL


Ncor <- cor.test(x = log10(integrated_df_molten_basic_atoms$expression_level+0.5), y = integrated_df_molten_basic_atoms$Nitrogen*3, method = "pearson", conf.level = 0.95)
p_Ncor <- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = expression_level+0.5, y = Nitrogen*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_smooth(method = "lm", colour="black", lty=2)+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 5.5, "pt"))+
  ggtitle(label = "i", subtitle = "br_1")+
  xlab(label = "mRNA abundance")+
  ylab(label = "Nitrogen atoms per codon")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  annotate("text", x = 4, y = 12.4, label = paste0("UCL = ", round(Ncor$conf.int[2], digits = 2)))+
  annotate("text", x = 4, y = 12.2, label = paste0("italic(R) == ", round(Ncor$estimate[[1]], digits = 2)), parse = TRUE)+
  annotate("text", x = 4, y = 12.0, label = paste0("LCL = ", round(Ncor$conf.int[1], digits = 2)))+
  annotation_logticks(sides = "b")+
  NULL

Ccor <- cor.test(x = log10(integrated_df_molten_basic_atoms$expression_level+0.5), y = integrated_df_molten_basic_atoms$Carbon*3, method = "pearson", conf.level = 0.95)
p_Ccor <- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = expression_level+0.5, y = Carbon*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_smooth(method = "lm", colour="black", lty=2)+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 5.5, "pt"))+
  ggtitle(label = "")+
  xlab(label = "mRNA abundance")+
  ylab(label = "Carbon atoms per codon")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  annotate("text", x = 10, y = 14.2, label = paste0("UCL = ", round(Ccor$conf.int[2], digits = 2)))+
  annotate("text", x = 10, y = 14.1, label = paste0("italic(R) == ", round(Ccor$estimate[[1]], digits = 2)), parse = TRUE)+
  annotate("text", x = 10, y = 14.0, label = paste0("LCL = ", round(Ccor$conf.int[1], digits = 2)))+
  annotation_logticks(sides = "b")+
  NULL


Ocor <- cor.test(x = log10(integrated_df_molten_basic_atoms$expression_level+0.5), y = integrated_df_molten_basic_atoms$Oxygen*3, method = "pearson", conf.level = 0.95)
p_Ocor <- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = expression_level+0.5, y = Oxygen*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_smooth(method = "lm", colour="black", lty=2)+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 5.5, "pt"))+
  ggtitle(label = "")+
  xlab(label = "mRNA abundance")+
  ylab(label = "Oxygen atoms per codon")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  annotate("text", x = 5, y = 3.8, label = paste0("UCL = ", round(Ocor$conf.int[2], digits = 2)))+
  annotate("text", x = 5, y = 3.7, label = paste0("italic(R) == ", round(Ocor$estimate[[1]], digits = 2)), parse = TRUE)+
  annotate("text", x = 5, y = 3.6, label = paste0("LCL = ", round(Ocor$conf.int[1], digits = 2)))+
  annotation_logticks(sides = "b")+
  NULL


Hcor <- cor.test(x = log10(integrated_df_molten_basic_atoms$expression_level+0.5), y = integrated_df_molten_basic_atoms$Hydrogen*3, method = "pearson", conf.level = 0.95)
p_Hcor <- ggplot(data = integrated_df_molten_basic_atoms, mapping = aes(x = expression_level+0.5, y = Hydrogen*3)) + # times 3 to calulate for codons
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 1, no = 0),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "purple", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "green", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "yellow", no = "grey"))),
             colour = "grey",
             pch    = 21)+
  geom_point(alpha  = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 1),
             size   = ifelse(test = integrated_df_molten_basic_atoms$operon == "other", yes = 0, no = 3),
             fill   = ifelse(test = integrated_df_molten_basic_atoms$operon == "pmo", yes = "#984EA3", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "mxa", yes = "#4DAF4A", no = ifelse(test = integrated_df_molten_basic_atoms$operon == "xox", yes = "#FF7F00", no = "grey"))),
             colour = "black",
             pch    = 21)+
  geom_smooth(method = "lm", colour="black", lty=2)+
  theme_classic2()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 200, 5.5, 5.5, "pt"))+
  ggtitle(label = "")+
  xlab(label = "mRNA abundance")+
  ylab(label = "Hydrogen atoms per codon")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  annotate("text", x = 5, y = 14.60, label = paste0("UCL = ", round(Hcor$conf.int[2], digits = 2)))+
  annotate("text", x = 5, y = 14.55, label = paste0("italic(R) == ", round(Hcor$estimate[[1]], digits = 2)), parse = TRUE)+
  annotate("text", x = 5, y = 14.50, label = paste0("LCL = ", round(Hcor$conf.int[1], digits = 2)))+
  annotation_logticks(sides = "b")+
  NULL

p_cors <- ggarrange(p_Ncor, p_Ccor, p_Ocor, p_Hcor, ncol=4)


atoms_pmo_rest <- integrated_df_molten_basic_atoms
atoms_pmo_rest$Carbon <- atoms_pmo_rest$Carbon * 3
atoms_pmo_rest$Hydrogen <- atoms_pmo_rest$Hydrogen * 3
atoms_pmo_rest$Oxygen <- atoms_pmo_rest$Oxygen * 3
atoms_pmo_rest$Nitrogen <- atoms_pmo_rest$Nitrogen * 3

atoms_pmo_rest_sampled <- atoms_pmo_rest[atoms_pmo_rest$operon=="pmo",]

set.seed(9301)
for (i in 1:1000) {
  atoms_rest_sampled <- atoms_pmo_rest[atoms_pmo_rest$operon!="pmo",][sample(x = nrow(atoms_pmo_rest[atoms_pmo_rest$operon!="pmo",]), size = 3), ]
  atoms_rest_sampled$operon <- paste0("rest", i)
  atoms_pmo_rest_sampled <- rbind(atoms_pmo_rest_sampled, atoms_rest_sampled)
}

N_of_random_samples <- 1000

p_stats_Pyri <- ggplot(data = atoms_pmo_rest_sampled, mapping = aes(x = operon, y = nPyrimidines, group=1)) +
  stat_summary(fun.y = "mean", fun.ymax = "mean", fun.ymin = "mean", colour = c("#984EA3", rep("grey50", N_of_random_samples)), size=c(0.5, rep(0.25, N_of_random_samples)))+
  geom_hline(yintercept = mean(atoms_pmo_rest_sampled$nPyrimidines[atoms_pmo_rest_sampled$operon=="pmo"]), colour="#984EA3", lty=2)+
  annotate(geom = "text", x = 80, y=0.54, label= expression(italic("pmoCAB")), colour="#984EA3")+
  theme_classic2()+
  xlab(label = "mean of pmoCAB and of 1000 random samples of n = 3")+
  ylab(label = "Pyrimidines fraction")+
  ggtitle(label = "c")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p_stats_C <- ggplot(data = atoms_pmo_rest_sampled, mapping = aes(x = operon, y = Carbon, group=1)) +
  stat_summary(fun.y = "mean", fun.ymax = "mean", fun.ymin = "mean", colour = c("#984EA3", rep("grey50", N_of_random_samples)), size=c(0.5, rep(0.25, N_of_random_samples)))+
  geom_hline(yintercept = mean(atoms_pmo_rest_sampled$Carbon[atoms_pmo_rest_sampled$operon=="pmo"]), colour="#984EA3", lty=2)+
  annotate(geom = "text", x = 80, y=13.4, label= expression(italic("pmoCAB")), colour="#984EA3")+
  theme_classic2()+
   xlab(label = "mean of pmoCAB and of 1000 random samples of n = 3")+
  ylab(label = "Carbon atom / codon")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p_stats_H <- ggplot(data = atoms_pmo_rest_sampled, mapping = aes(x = operon, y = Hydrogen, group=1)) +
  stat_summary(fun.y = "mean", fun.ymax = "mean", fun.ymin = "mean", colour = c("#984EA3", rep("grey50", N_of_random_samples)), size=c(0.5, rep(0.25, N_of_random_samples)))+
  geom_hline(yintercept = mean(atoms_pmo_rest_sampled$Hydrogen[atoms_pmo_rest_sampled$operon=="pmo"]), colour="#984EA3", lty=2)+
  annotate(geom = "text", x = 80, y=14.05, label= expression(italic("pmoCAB")), colour="#984EA3")+
  theme_classic2()+
   xlab(label = "mean of pmoCAB and of 1000 random samples of n = 3")+
  ylab(label = "Hydrogen atom / codon")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p_stats_O <- ggplot(data = atoms_pmo_rest_sampled, mapping = aes(x = operon, y = Oxygen, group=1)) +
  stat_summary(fun.y = "mean", fun.ymax = "mean", fun.ymin = "mean", colour = c("#984EA3", rep("grey50", N_of_random_samples)), size=c(0.5, rep(0.25, N_of_random_samples)))+
  geom_hline(yintercept = mean(atoms_pmo_rest_sampled$Oxygen[atoms_pmo_rest_sampled$operon=="pmo"]), colour="#984EA3", lty=2)+
  annotate(geom = "text", x = 80, y=3.35, label= expression(italic("pmoCAB")), colour="#984EA3")+
  theme_classic2()+
  xlab(label = "mean of pmoCAB and of 1000 random samples of n = 3")+
  ylab(label = "Oxygen atom / codon")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p_stats_N <- ggplot(data = atoms_pmo_rest_sampled, mapping = aes(x = operon, y = Nitrogen, group=1)) +
  stat_summary(fun.y = "mean", fun.ymax = "mean", fun.ymin = "mean", colour = c("#984EA3", rep("grey50", N_of_random_samples)), size=c(0.5, rep(0.25, N_of_random_samples)))+
  geom_hline(yintercept = mean(atoms_pmo_rest_sampled$Nitrogen[atoms_pmo_rest_sampled$operon=="pmo"]), colour="#984EA3", lty=2)+
  annotate(geom = "text", x = 80, y=10.85, label= expression(italic("pmoCAB")), colour="#984EA3")+
  theme_classic2()+
  xlab(label = "mean of pmoCAB and of 1000 random samples of n = 3")+
  ylab(label = "Nitrogen atom / codon")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# integrated figure
lay_matrix <- rbind(c(1, 2, NA),
                    c(1,  2, 8),
                    c(NA, 2, NA),
                    c(3, 4,  5),
                    c(3, 4,  5),
                    c(3, 4,  5),
                    c(6, 7,  9),
                    c(6, 7,  10),
                    c(11, 12, NA),
                    c(13, 13, 13),
                    c(13, 13, 13))
ggsave(filename = "0_puriPiry_20Z.pdf",
       plot   = grid.arrange(p_puri_pyri_hist_genome_20Z, p_pyri_exp_20Z, p_base_20Z, p_gc_20Z, p_compound_20Z, p_atoms_NC_20Z, p_atoms_OH_20Z, p_stats_Pyri, p_stats_C, p_stats_H, p_stats_O, p_stats_N, p_cors, layout_matrix = lay_matrix),
       device = "pdf",
       width  = 25,
       height = 20,
       useDingbats=FALSE)

## With one dot per operon and mean and SD per operon
integrated_df_molten_basic_atoms_only_module <- integrated_df_molten_basic_atoms[integrated_df_molten_basic_atoms$operon != "other", ]
integrated_df_molten_basic_atoms_only_module$Carbon <- integrated_df_molten_basic_atoms_only_module$Carbon*3
integrated_df_molten_basic_atoms_only_module$Oxygen <- integrated_df_molten_basic_atoms_only_module$Oxygen*3

atoms_C_per_codon_in_module_df <- data.frame(operon = "pmo",
                                             Mean = mean(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "pmo"]),
                                             SD = sd(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "pmo"]))

atoms_C_per_codon_in_module_df <- rbind(atoms_C_per_codon_in_module_df,
                                        data.frame(operon = "mxa",
                                                   Mean = mean(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "mxa"]),
                                                   SD = sd(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "mxa"])))
atoms_C_per_codon_in_module_df <- rbind(atoms_C_per_codon_in_module_df,
                                        data.frame(operon = "xox",
                                                   Mean = mean(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "xox"]),
                                                   SD = sd(integrated_df_molten_basic_atoms_only_module$Carbon[integrated_df_molten_basic_atoms_only_module$operon == "xox"])))
colnames(atoms_C_per_codon_in_module_df)[c(2,3)] <- c("Carbon_mean", "Carbon_sd")


atoms_O_per_codon_in_module_df <- data.frame(operon2 = "pmo",
                                             Mean = mean(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "pmo"]),
                                             SD = sd(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "pmo"]))

atoms_O_per_codon_in_module_df <- rbind(atoms_O_per_codon_in_module_df,
                                        data.frame(operon2 = "mxa",
                                                   Mean = mean(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "mxa"]),
                                                   SD = sd(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "mxa"])))
atoms_O_per_codon_in_module_df <- rbind(atoms_O_per_codon_in_module_df,
                                        data.frame(operon2 = "xox",
                                                   Mean = mean(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "xox"]),
                                                   SD = sd(integrated_df_molten_basic_atoms_only_module$Oxygen[integrated_df_molten_basic_atoms_only_module$operon == "xox"])))
colnames(atoms_O_per_codon_in_module_df)[c(2,3)] <- c("Oxygen_mean", "Oxygen_sd")

atoms_OC_per_codon_in_module_df <- cbind(atoms_O_per_codon_in_module_df, atoms_C_per_codon_in_module_df)

p_atoms_OC_20Z_mean_sd <- ggplot(data = atoms_OC_per_codon_in_module_df, mapping = aes(x = Carbon_mean, y = Oxygen_mean)) + # times 3 to calulate for codons
  geom_hline(yintercept = mean(integrated_df_molten_basic_atoms$Oxygen)*3, lty=2) +
  geom_vline(xintercept = mean(integrated_df_molten_basic_atoms$Carbon)*3, lty=2) +
  geom_errorbar(mapping=aes(x=Carbon_mean, ymin=Oxygen_mean-Oxygen_sd, ymax=Oxygen_mean+Oxygen_sd), width = 0.02, size=0.5) +
  geom_errorbarh(mapping=aes(y=Oxygen_mean, xmin=Carbon_mean-Carbon_sd, xmax=Carbon_mean+Carbon_sd), height = 0.02, size=0.5) +
  annotate("text", x = 13.67, y = 2.965, label = "Genome mean")+
  geom_point(fill   = c("#984EA3", "#4DAF4A", "#FF7F00"),
             colour = "black",
             pch    = 21,
             size = 5)+
  theme_linedraw()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = expression(paste(italic("Methylomicrobium alcaliphilum"), " 20Z")))+
  xlab(label = "Carbon atoms per codon")+
  ylab(label = "Oxygen atoms per codon")+
  NULL
ggsave(filename = "0_countergradient_mean_sd_20Z.pdf",
       plot = p_atoms_OC_20Z_mean_sd,
       device = "pdf",
       width  = 5,
       height = 5,
       useDingbats=FALSE)

p_delta_atoms_OC_20Z_mean_sd <- ggplot(data = atoms_OC_per_codon_in_module_df, mapping = aes(x = Carbon_mean-(mean(integrated_df_molten_basic_atoms$Carbon)*3), y = Oxygen_mean-(mean(integrated_df_molten_basic_atoms$Oxygen)*3))) + # times 3 to calulate for codons
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  geom_errorbar(mapping=aes(x=Carbon_mean-(mean(integrated_df_molten_basic_atoms$Carbon)*3), ymin=(Oxygen_mean-Oxygen_sd)-(mean(integrated_df_molten_basic_atoms$Oxygen)*3), ymax=Oxygen_mean+Oxygen_sd-(mean(integrated_df_molten_basic_atoms$Oxygen)*3)), width = 0.02, size=0.5) +
  geom_errorbarh(mapping=aes(y=Oxygen_mean-(mean(integrated_df_molten_basic_atoms$Oxygen)*3), xmin=(Carbon_mean-Carbon_sd)-(mean(integrated_df_molten_basic_atoms$Carbon)*3), xmax=Carbon_mean+Carbon_sd-(mean(integrated_df_molten_basic_atoms$Carbon)*3)), height = 0.02, size=0.5) +
  annotate("text", x = 0.15, y = 0.019, label = "Genome mean")+
  geom_point(fill   = c("#984EA3", "#4DAF4A", "#FF7F00"),
             colour = "black",
             pch    = 21,
             size = 5)+
  theme_linedraw()+
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggtitle(label = expression(paste(italic("Methylomicrobium alcaliphilum"), " 20Z")))+
  xlab(label =  expression(Delta*" Carbon atoms per codon relative to genome-scale mean"))+
  ylab(label = expression(Delta*" Oxygen atoms per codon relative to genome-scale mean"))+
  NULL
ggsave(filename = "0_delta_countergradient_mean_sd_20Z.pdf",
       plot = p_delta_atoms_OC_20Z_mean_sd,
       device = "pdf",
       width  = 5,
       height = 5,
       useDingbats=FALSE)