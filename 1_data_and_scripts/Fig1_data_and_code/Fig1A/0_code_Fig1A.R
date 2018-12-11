#install.packages("maps")
#install.packages("ggplot2")
library("maps")
library("ggplot2")


map <- map_data("world")
GEO_df <- read.delim(file = "1_Table_S1_genomes_CH4_GEO.txt")
manual_colors2 <- c("#CD9600", # Ia
                    "#7CAE00", # Ib
                    "#00BE67", # Ic
                    "#00BFC4", # IIa *
                    "#00A9FF", # IIb
                    "#C77CFF", # III
                    "#FF61CC") # NC10)
pdf(file = "Fig1A.pdf", width = 9, height = 9)
ggplot() + 
  geom_polygon(data=map, aes(x=long, y=lat, group=group), colour="grey", fill="grey") +
  coord_fixed() + 
  xlab("Longitude") + 
  ylab("Latitude") +
  geom_point(data=GEO_df, aes(x=Longitude_country, y=Latitude_country, fill=Group, shape=Sample),
             colour="black", size=3, alpha=0.7, position = position_jitter(w = 3.1, h = 3.1)) +
  scale_shape_manual(values=c(21, 24))+
  scale_fill_manual(values = manual_colors2) +
  guides(fill = guide_legend(title="Type", override.aes=list(shape=21)), shape = guide_legend(title="", override.aes=list(color="black"))) +
  theme_bw()
dev.off()
