#install.packages('seqinr')
#install.packages('plotrix')
#install.packages("robustbase")
#install.packages("gplots")
library('seqinr')
library("plotrix")
library("robustbase")
library("gplots")

genes_metadata <- read.delim(file = "QC_data_CH4_IV.txt")
ALL_query_CDSs <- read.fasta(file = "query_ALL_clean.fasta", seqtype = 'DNA')
dim(genes_metadata)[1] == length(ALL_query_CDSs)

rscu_all <- sapply(ALL_query_CDSs, uco, index = "rscu")
rscu_all <- t(rscu_all)
rscu_all[is.na(rscu_all)] <- 0
rscu_all <- rscu_all[,-c(15, 49,51,57,59)] # removing Met, Trp and Stop codons


#----------------- FIGURE 3D_Suppl -----------------#
df_rscu_all <- as.data.frame(rscu_all)
df_rscu_all$CDS <- rownames(df_rscu_all)
df_rscu_integrated_all <- merge(x = genes_metadata, y = df_rscu_all, by="CDS")
s2c_n_translate_fun <- function(x) {translate(s2c(x))}
colnames(df_rscu_integrated_all)[25:83] <- paste(aaa(sapply(X = colnames(df_rscu_integrated_all)[25:83], FUN = s2c_n_translate_fun)),
                                                 "_",
                                                 toupper(colnames(df_rscu_integrated_all)[25:83]),
                                                 sep = "")

pdf(file = "Fig3DS1.pdf", width = 20, height = 5, useDingbats=FALSE)
{par(mfrow=c(1,4))
# RSCU plotrix::radialplot of group Ia-pmo
sub_medians         <- colMedians(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group=="Ia" & (df_rscu_integrated_all$gene_category == "pmo/amo"), c(25:83)]))
sub_medians_ordered <- sub_medians[order(names(sub_medians))]
sub_medians_ordered[sub_medians_ordered > 1.6] # frequent codons in pmo
plotrix::radial.plot(lengths = sub_medians_ordered,
                     labels = names(sub_medians_ordered),
                     rp.type = "sr",
                     radlab = TRUE,
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 20),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     line.col = ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                       no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     label.prop = 1.18,
                     lwd = ifelse(test = sub_medians_ordered >= 1.6, yes = 3, no = 3),
                     mar = c(4,3,8,3),
                     boxed.radial=TRUE,
                     grid.col = "grey",
                     radial.lim = c(0, 4),
                     main = "pmo")
plotrix::radial.plot(lengths = rep(4, length(sub_medians_ordered)),
                     rp.type = "s",
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 19),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     cex = abs(1-sub_medians_ordered),
                     radial.lim=c(0, 4),
                     add = TRUE)
plotrix::radial.plot(lengths = rep(1.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#d73027",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)
plotrix::radial.plot(lengths = rep(0.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#4575b4",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)

# RSCU plotrix::radialplot of group Ia-mxa
sub_medians         <- colMedians(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group=="Ia" & (df_rscu_integrated_all$gene_category == "mxa"), c(25:83)]))
sub_medians_ordered <- sub_medians[order(names(sub_medians))]
sub_medians_ordered[sub_medians_ordered > 1.6] # frequent codons in mxa
plotrix::radial.plot(lengths = sub_medians_ordered,
                     labels = names(sub_medians_ordered),
                     rp.type = "sr",
                     radlab = TRUE,
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 20),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     line.col = ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                       no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     label.prop = 1.18,
                     lwd = ifelse(test = sub_medians_ordered >= 1.6, yes = 3, no = 3),
                     mar = c(4,3,8,3),
                     boxed.radial=TRUE,
                     grid.col = "grey",
                     radial.lim = c(0, 4),
                     main = "mxa")
plotrix::radial.plot(lengths = rep(4, length(sub_medians_ordered)),
                     rp.type = "s",
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 19),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     cex = abs(1-sub_medians_ordered),
                     radial.lim=c(0,4),
                     add = TRUE)
plotrix::radial.plot(lengths = rep(1.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#d73027",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0,4), 
                     add = TRUE)
plotrix::radial.plot(lengths = rep(0.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#4575b4",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0,4), 
                     add = TRUE)

# RSCU plotrix::radialplot of group Ia-(s)mmo
sub_medians         <- colMedians(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group=="Ia" & (df_rscu_integrated_all$gene_category == "mmo"), c(25:83)]))
sub_medians_ordered <- sub_medians[order(names(sub_medians))]
sub_medians_ordered[sub_medians_ordered > 1.6] # frequent codons in (s)mmo
plotrix::radial.plot(lengths = sub_medians_ordered,
                     labels = names(sub_medians_ordered),
                     rp.type = "sr",
                     radlab = TRUE,
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 20),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     line.col = ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                       no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     label.prop = 1.18,
                     lwd = ifelse(test = sub_medians_ordered >= 1.6, yes = 3, no = 3),
                     mar = c(4,3,8,3),
                     boxed.radial=TRUE,
                     grid.col = "grey",
                     radial.lim = c(0, 4),
                     main = "mmo")
plotrix::radial.plot(lengths = rep(4, length(sub_medians_ordered)),
                     rp.type = "s",
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 19),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     cex = abs(1-sub_medians_ordered),
                     radial.lim=c(0, 4),
                     add = TRUE)
plotrix::radial.plot(lengths = rep(1.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#d73027",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)
plotrix::radial.plot(lengths = rep(0.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#4575b4",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)

# RSCU plotrix::radialplot of group Ia-xox
sub_medians         <- colMedians(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group=="Ia" & (df_rscu_integrated_all$gene_category == "xox"), c(25:83)]))
sub_medians_ordered <- sub_medians[order(names(sub_medians))]
sub_medians_ordered[sub_medians_ordered > 1.6] # frequent codons in xox
plotrix::radial.plot(lengths = sub_medians_ordered,
                     labels = names(sub_medians_ordered),
                     rp.type = "sr",
                     radlab = TRUE,
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 20),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     line.col = ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                       no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     label.prop = 1.18,
                     lwd = ifelse(test = sub_medians_ordered >= 1.6, yes = 3, no = 3),
                     mar = c(4,3,8,3),
                     boxed.radial=TRUE,
                     grid.col = "grey",
                     radial.lim = c(0, 4),
                     main = "xox")
plotrix::radial.plot(lengths = rep(4, length(sub_medians_ordered)),
                     rp.type = "s",
                     point.symbols = ifelse(test = sub_medians_ordered >= 1.6, yes = 19, no = 19),
                     point.col =  ifelse(test = sub_medians_ordered >= 1.6, yes = "#d73027",
                                         no = ifelse(test = sub_medians_ordered <= 0.6, yes = "#4575b4", no = "grey")),
                     cex = abs(1-sub_medians_ordered),
                     radial.lim=c(0, 4),
                     add = TRUE)
plotrix::radial.plot(lengths = rep(1.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#d73027",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)
plotrix::radial.plot(lengths = rep(0.6, length(sub_medians_ordered)),
                     rp.type = "p",
                     line.col = "#4575b4",
                     lty=2,
                     lwd=2,
                     radial.lim=c(0, 4), 
                     add = TRUE)
dev.off()
}

#----------------- FIGURE Fig3D Suppl 2 -----------------#
{my_palette <- colorRampPalette(c("grey80", "#f7f7f7", "#ca0020"))(n = 23)
col_breaks = c(seq(0.00, 0.60, length=6),
               seq(0.61, 1.60, length=6),
               seq(1.61, 6.00, length=12))
gr_row        <- df_rscu_integrated_all$gene_category[df_rscu_integrated_all$group=="Ia"]
gr_row_labels <- df_rscu_integrated_all$gene_category_subunit[df_rscu_integrated_all$group=="Ia"]
col1 <- ifelse(test = gr_row=="mmo", yes = "#d9ef8b", 
               ifelse(test = gr_row=="hao", yes = "orange", 
                      ifelse(test = gr_row=="mxa", yes = "#4575b4", 
                             ifelse(test = gr_row=="xox", yes = "#b2abd2", 
                                    ifelse(test = gr_row=="pxm", yes = "#c51b8a", 
                                           ifelse(test = gr_row=="pmo/amo", yes = "#1a9850", "red"))))))
pdf(file = "Fig3DS2.pdf", width = 10, height = 10)
heatmap.2(t(as.matrix(df_rscu_integrated_all[df_rscu_integrated_all$group=="Ia", c(25:83)])),
          main = "Type Ia",
          notecol = "black",
          density.info = "none",
          trace = "none",
          margins =c(15,5),
          col=my_palette,
          breaks=col_breaks,
          rowsep = c(1:60),
          labCol = FALSE,
          ColSideColors = col1,
          keysize=0.6, key.title = "", key.xlab = "RSCU",
          key.par = list(cex=0.6),
          xlab = NULL,
          hclustfun = function(x) hclust(x,method = 'complete'))
legend("bottom",      
       legend = c("pmo",     "mmo",     "mxa",     "xox",     "hao",   "pxm"),
       col    = c("#1a9850", "#d9ef8b", "#4575b4", "#b2abd2", "orange", "#c51b8a"), 
       lty    = 1,             
       lwd    = 8,           
       cex    = 1,
       ncol   = 6)
dev.off()