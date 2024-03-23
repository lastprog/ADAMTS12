library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(mgsub)
setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/2) ADAMTS12 Expression in PDGFRB and CD10 neg scRNA Datasets")

# Load PDGFRb dataset
sample="huCKD_PDGFRb_"
Samples.combined <- readRDS("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/alt/huCKD_Nature_Decode_Myofibroblast_Dataset/Human_PDGFRBpositive.rds")

# Basic processing
Samples.combined <- NormalizeData(Samples.combined, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE, features = rownames(Samples.combined))
Samples.combined <- RunPCA(Samples.combined, verbose = FALSE)

# Introduce Abbreviations
Samples.combined@meta.data[["Annotation.Level.2"]] <- 
  mgsub(Samples.combined@meta.data[["Annotation.Level.2"]],
        c("Glomerular Capillaries",  "Venular Endothelium", "Thick Ascending Limb", "Injured Proximal Tubule", "Prinicipal Cells", "Intercalated Cells",
          "Proximal Tubule", "Vascular Smooth Muscle Cells"),
        c("Glom Cap", "Venular Endo", "TAL", "Injured PT", "PC", "IC", "PT", "VSMC"))

Samples.combined@meta.data[["Annotation.Level.3"]] <- 
  mgsub(Samples.combined@meta.data[["Annotation.Level.3"]],
        c("Glomerular Capillaries 1","Vasa Recta","Glomerular Capillaries 2","Arterioles 1","Venular Endothelium",
          "Arterioles 2","Thick Ascending Limb","Injured Proximal Tubule","Collecting Duct Principal Cells","Intercalated Cells",
          "Proximal Tubule","Fibroblast 6","Injured Monocytes","Pericytes","Vascular Smooth Muscle Cells 2", 
          "Myofibroblasts 2b","Myofibroblasts 1","Myofibroblasts 3a","Fibroblasts 4a","Fibroblasts 4b",
          "Myofibroblasts 2a","Mesangial Cells 1","Fibroblast 2","Myofibroblasts 3b","Dendritic Cells",
          "Monocytes"),
        c("Glom Cap 1", "VR", "Glom Cap 2", "Art 1", "Ven Endo", "Art 2", "TAL", "IPT", 
          "CDPC", "IC", "PT", "Fib 6", "Inj Mono", "Peri", "VSMC 2" ,"MyoF 2b", 
          "MyoF 1", "MyoF 3a", "Fib 4a", "Fib 4b", "MyoF 2a", "Mesangial 1", 
          "Fib 2", "MyoF 2b", "Dendritic", "Mono"))

# Annotation level 2
Idents(Samples.combined)<-"Annotation.Level.2"

# Figure 1E Featureplot
FeaturePlot(Samples.combined, features ="ADAMTS12", pt.size = 0.8, order=TRUE, label = TRUE)+ scale_colour_gradient2(mid = 'lightgrey', high = 'red')+
  theme_void()+
  theme(title = element_blank())
ggsave(filename = paste0(sample,"_ADAMTS12.jpeg"), width=7 , height = 7)

# Figure 1E Dotplot
DotPlot(Samples.combined, features = "ADAMTS12",dot.scale = 3)+  
  theme(plot.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"_ADAMTS12_Level2_Dotplot.svg"),width=4.2, height = 3)


Samples.combined@active.ident <- factor(Samples.combined@active.ident, levels = rev(c("TAL", "PC", "IC", "Glom Cap", "Mesangial Cells", "Dendritic Cells", "VSMC",
                                                                                  "Vasa Recta", "Arterioles", "Venular Endo", "Monocytes", "Pericytes", "Injured PT",
                                                                                  "PT", "Fibroblasts", "Myofibroblasts")))

# Supplemental Figure 2 A
DotPlot(Samples.combined, features = c("POSTN", "COL1A1", "PDGFRB", "ACTA2"), dot.scale = 3)+  
  theme(plot.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2))+
  coord_flip()
ggsave(filename = paste0(sample,"_ADAMTS12_Level2_MarkerGenes_Dotplot.jpg"),width=4.6, height = 3, bg = "white")

# Annotation level 3
Idents(Samples.combined) <- "Annotation.Level.3"

# Supplemental Figure 2B
FeaturePlot(Samples.combined, features ="ADAMTS12", pt.size = 0.8, order=TRUE, label = TRUE)+ scale_colour_gradient2(mid = 'lightgrey', high = 'red')+
  theme_void()+
  theme(title = element_blank())
ggsave(filename = paste0(sample,"_ADAMTS12_ann_lvl3.jpeg"), width=7 , height = 7)

Samples.combined@active.ident <- factor(Samples.combined@active.ident, levels = c("MyoF 1", "MyoF 2a", "MyoF 2b", "MyoF 3a", "MyoF 3b", "Fib 2", "Fib 4a",
                                                                                  "Fib 4b", "Fib 6", "VSMC 2", "Pericytes", "Mesangial 1", "Mono",
                                                                                  "Dendritic", "Inj Mono", "PT", "IC", "CDPC", "IPT", "TAL", "Ven Endo",
                                                                                  "Art 2", "Art 1", "VR", "Glom Cap 2", "Glom Cap 1"))
# Supplemental Figure 2C
DotPlot(Samples.combined, features = "ADAMTS12",dot.scale = 3)+  
  theme(plot.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, size =8, face= "plain", family= "Arial", hjust =1),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))+
  coord_flip()
ggsave(filename = paste0(sample,"_ADAMTS12_Level3_Dotplot.svg"),width=4.2, height = 3)

DotPlot(Samples.combined, features = c("POSTN", "COL1A1", "PDGFRB", "ACTA2"), dot.scale = 3)+  
  theme(plot.title = element_blank(), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2))+
  coord_flip()
ggsave(filename = paste0(sample,"_Level3_MarkerGenes_Dotplot.svg"),width=4.6, height = 3, bg = "white")



