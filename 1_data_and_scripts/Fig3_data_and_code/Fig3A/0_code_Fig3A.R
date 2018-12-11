#install.packages("gplots")
#install.packages("RColorBrewer")
library("gplots")
library("RColorBrewer")

trna_matrix <- as.matrix(read.delim("1_tRNA_anticodon_CH4.txt", header = TRUE, row.names = 1))

my_palette <- colorRampPalette(c("white", "#b3cde3", "#8c96c6", "#8856a7", "#810f7c"))

col_breaks = c(0, 0.5, 1.5, 2.5, 4.5, 7)

pdf(file = "Fig3A.pdf", width = 8, height = 10)
heatmap.2(trna_matrix,
          main = "tRNA anticodon", # heat map title
          notecol = "black",      # change font color of cell labels to black
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins =c(15,5),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,
          colsep = c(0:dim(trna_matrix)[2]),
          rowsep = c(0:dim(trna_matrix)[1]),
          #sepcolor = "black",
          #dendrogram= "none",
          Colv="NA",
          #Rowv = "NA", 
          keysize=0.6,
          key.par = list(cex=0.6))
dev.off()

pdf(file = "Fig3AS1.pdf", width = 8, height = 10)
heatmap.2(trna_matrix,
          main = "tRNA anticodon", # heat map title
          cellnote=trna_matrix,
          notecex=0.6,
          notecol = "black",      # change font color of cell labels to black
          density.info = "none",  # turns off density plot inside color legend
          trace = "none",         # turns off trace lines inside the heat map
          margins =c(15,5),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,
          colsep = c(0,23,29, 30, 41, 46, 49, 50, 68), #colsep = c(1:dim(mm)[2]),
          rowsep = c(0, 4, 10, 12, 14, 16, 18, 20, 24, 26, 29, 35, 37, 38, 40, 44, 48, 51, 57, 58, 60, 64), #rowsep = c(1:dim(mm)[1]),
          sepcolor = "black",
          dendrogram= "none",
          Colv="NA",
          Rowv = "NA", 
          keysize=0.6,
          key.par = list(cex=0.6))
dev.off()