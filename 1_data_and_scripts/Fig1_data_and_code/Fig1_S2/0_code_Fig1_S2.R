#install.packages("BiocManager")
#BiocManager::install("ape")
#BiocManager::install("ggtree")
#install.packages('ggridges')
#install.packages("ggrepel")
#install.packages("ggstance")
library('ape')
library('ggtree')
library('ggridges')
library("ggrepel")
library('ggstance')

tree_raw <- read.tree(file = "1_proteomes_tree.nwk")
tree_rooted <- root(tree_raw, outgroup = "Bacteroides_ovatus_ATCC_8483", edgelabel = TRUE, resolve.root=TRUE)

genomes_list <- read.delim(file = "2_genomes_list.txt", header = T, sep = '\t')[,1:2]
rownames(genomes_list) <- genomes_list[,1]

groupInfo <- list(Ia=rownames(genomes_list)[genomes_list$group=='Ia'],
                  Ib=rownames(genomes_list)[genomes_list$group=='Ib'],
                  #Ic=rownames(genomes_list)[genomes_list$group=='Ic'],
                  IIa=rownames(genomes_list)[genomes_list$group=='IIa'],
                  IIb=rownames(genomes_list)[genomes_list$group=='IIb'],
                  III=rownames(genomes_list)[genomes_list$group=='III'])

tree_grouped <- groupOTU(.data = tree_rooted, .node = groupInfo)

all_df <- data.frame(read.table(file = "3_all_dataframe.txt", header = TRUE, stringsAsFactors = FALSE))
all_df <- all_df[,c(2,1,3,4,5)]

p_tree <- ggtree(tree_grouped) + geom_tiplab(size=3, offset = 0.03, align=TRUE, linesize=0.25)
length_df <- as.data.frame(table(all_df$BIN))
colnames(length_df) <- c('strain', 'Freq')
borrow_QC <- read.delim(file = "4_QC_data_IV.txt")[,c(1,3)]
borrow_QC <- unique(borrow_QC)
outgroup_sample <- data.frame(sample = "Isolate", strain="Bacteroides_ovatus_ATCC_8483")
borrow_QC <- rbind(borrow_QC, outgroup_sample)
length_df <- merge(x = length_df, y = borrow_QC, by = 'strain', all.x = TRUE)

#library('scales')      #to pick manual colors
#show_col(hue_pal()(8)) #to pick manual colors

#--------------- Fig1B --------------- #
p_GC <- facet_plot(p = p_tree, 
                   panel="GC/GC3 distribution (fraction)", 
                   data = all_df,
                   geom_density_ridges,
                   mapping = aes(x = GC, group=label, fill=group, alpha=0.8), 
                   colour="black", 
                   lwd = 0.3) 

p_GC3 <- facet_plot(p_GC, 
                    panel="GC/GC3 distribution (fraction)", 
                    data = all_df, 
                    geom_density_ridges, 
                    mapping = aes(x=GC3, group=label, fill=group, alpha=0.6),
                    colour="black",
                    lwd=0.3)

manual_colors2 <- c("grey",    # outgroup,
                    "#CD9600", # Ia
                    "#7CAE00", # Ib
                    "#00BFC4", # IIa
                    "#00A9FF", # IIb
                    "#C77CFF") # III

pdf(file = "Fig1_S2.pdf", width = 7, height = 10, useDingbats=FALSE)
p_GC3 + scale_fill_manual(values = manual_colors2) + scale_color_manual(values = manual_colors2) + theme_tree2(legend.position="bottom")
dev.off()
#--------------- end --------------- #




#--------------- Fig1BS1 --------------- #
p_GS <- facet_plot(p_tree, 
                   panel="Number of CDSs / 1000", 
                   data = length_df,
                   geom_point, 
                   mapping = aes(x=Freq/1E3, size=0.81, fill=group, alpha=0.7, pch=sample))

p_GC <- facet_plot(p = p_GS, 
                   panel="GC/GC3 distribution (fraction)", 
                   data = all_df,
                   geom_density_ridges,
                   mapping = aes(x = GC, group=label, fill=group, alpha=0.8), 
                   colour="black", 
                   lwd = 0.3) 

p_GC3 <- facet_plot(p_GC, 
                    panel="GC/GC3 distribution (fraction)", 
                    data = all_df, 
                    geom_density_ridges, 
                    mapping = aes(x=GC3, group=label, fill=group, alpha=0.6),
                    colour="black",
                    lwd=0.3)

p_N  <- facet_plot(p_GC3, 
                    panel="n", 
                    data = length_df, 
                    geom_text, 
                    mapping = aes(x=0, label=paste0("n=",Freq)),
                    colour="black")
pdf(file = "Fig1_S2b.pdf", width = 7, height = 10, useDingbats=FALSE)
p_N + scale_fill_manual(values = manual_colors2) + scale_color_manual(values = manual_colors2) + theme_tree2(legend.position="bottom")
dev.off()
#--------------- end --------------- #

