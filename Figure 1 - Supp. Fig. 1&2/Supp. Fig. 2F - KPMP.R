library(Seurat)
library(SeuratDisk)
library(ggplot2)
setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/1) ADAMTS12 expression in other Datasets/LK for publication")

# Read dataset
sample="KPMP_"
Samples.combined <- LoadH5Seurat("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/KPMP/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")

# Annotation level
Idents(Samples.combined)<-"subclass.l2"

# Supplemental Figure 2F
DotPlot(Samples.combined, features = c("POSTN", "COL1A1", "PDGFRB", "ACTA2"),dot.scale = 3)+ 
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"))+
  coord_flip()
ggsave(filename = paste0(sample,"_ADAMTS12_L2_MarkerGenes_Dotplot_l.svg"), width= 6.5, height = 1.5)

# Figure 1I - KPMP ADAMTS12 by disease
Samples.combined$subclass.l2_diseasetype <- paste(Samples.combined$subclass.l2, Samples.combined$diseasetype, sep = "_")
Idents(Samples.combined) <- "subclass.l2_diseasetype"
ADAMTS12_Data <- DotPlot(Samples.combined, features = "ADAMTS12")
df <- ADAMTS12_Data$data %>% separate("id", into = c("Annotation", "diseasetype"),sep = "_" )
ggplot(data=df)+ 
  geom_point(mapping = aes(x= diseasetype, y=Annotation, color=avg.exp, size=pct.exp))+
  scale_x_discrete(limits = c("LivingDonor", "AKI", "CKD"))+
  bgcolor("white")+
  scale_colour_gradient2(mid = "lightgrey" , high = "blue")+
  theme_cowplot()+
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), 
        axis.title = element_blank(),
       axis.text.x = element_text(angle = 90, hjust = 1, size =6, face= "plain", family= "Arial"),
       axis.text.y = element_text(angle = 45, size = 8, face= "plain", family = "Arial"),
       legend.title = element_text(size =8, face= "plain", family = "Arial"))+
  coord_flip()
ggsave(filename = paste0(sample,"_ADAMTS12_L2_Dotplot1.svg"),width= 6.5, height = 1.0, bg = "white")

