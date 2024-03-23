library(Seurat)
library(ggplot2)
library(cowplot)
setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/1) ADAMTS12 expression in other Datasets/")

# Read dataset
sample<-"Kirita"
Samples.combined <- LoadSeuratRds("Kirita.rds")

# Supplemental Figure 1K
FeaturePlot(Samples.combined, features ="Adamts12", pt.size = 1.3, order=TRUE, label = TRUE, label.size = 2)+ scale_colour_gradient2(mid = 'lightgrey', high = 'red')+
  NoLegend()+ 
  theme_void()+
  theme(legend.position = "none",
        title = element_blank())
ggsave(filename = paste0("LK for publication/", sample,"_ADAMTS12_red_l.svg"), width= 1.7 , height = 1.7)

# Split time
Samples.combined$time <- paste(Samples.combined$name, Samples.combined$time, sep = "_")
Idents(Samples.combined) <- "time"
ADAMTS12_Data <- DotPlot(Samples.combined, features = "Adamts12")
df <- ADAMTS12_Data$data %>% separate("id", into = c("name", "time"),sep = "_" )

# Supplemental Figure 1L
ggplot(data=df)+ 
  geom_point(mapping = aes(x= name, y=time, color=avg.exp, size=pct.exp))+
  scale_y_discrete(limits = rev(c("Control", "4hours", "12hours", "2days", "14days", "6weeks")))+
  bgcolor("white")+
  scale_colour_gradient2(mid = "lightgrey" , high = "blue")+
  scale_size_continuous(range = c(0.1,3))+
  theme_cowplot()+
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size = 8, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"))
ggsave(filename = paste0("LK for publication/", sample,"_ADAMTS12_L2_MarkerGenes_Dotplot_l_time.svg"), width= 5, height = 1.5)









