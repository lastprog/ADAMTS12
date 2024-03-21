library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(mgsub)
library(svglite)
library(tidyverse)

setwd("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/0) Figures, Tables and Western Blots/Figure 1 Daten und Abbildungen")
# Read dataset
Samples.combined.CD10neg <- readRDS("/Users/larskoch/Library/CloudStorage/Dropbox/3. Adamts12/Revision/2) ADAMTS12 Expression in PDGFRB and CD10 neg scRNA Datasets/dataset/Human_CD10negative.rds")
sample <- "CD10neg"

# Basic processing steps
Samples.combined.CD10neg <- NormalizeData(Samples.combined.CD10neg, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
Samples.combined.CD10neg <- FindVariableFeatures(Samples.combined.CD10neg, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined.CD10neg <- ScaleData(Samples.combined.CD10neg, verbose = FALSE, features = rownames(Samples.combined.CD10neg))
Samples.combined.CD10neg <- RunPCA(Samples.combined.CD10neg, verbose = FALSE)

# Introduce Abbreviations
Samples.combined.CD10neg@meta.data[["Annotation.Level.2"]] <- 
  mgsub(Samples.combined.CD10neg@meta.data[["Annotation.Level.2"]],
        c("Vasa Recta", "Glomerular Capillaries","Arterioral Endothelium","Venular Endothelium","Lymph Endothelium",
          "Injured Endothelium","Distal Convoluted Tubule","Injured Tubule","Proximal Tubule","Thick Ascending Limb",
          "Collecting Duct Principal Cells","Descending Thin Limb","Intercalated Cells","Connecting Tubule",
          "Uroethlial Cells","Smooth Muscle Cells","Fibroblast","Pericytes","Myofibroblast",
          "Podocytes","Schwann Cells","Plasma Cells","Dendritic Cells","Monocytes","Macrophages",
          "Basophils","T Cells","B Cells","Natural Killer Cells"),
        c("VR","Glom Cap", "Art Endo", "Venular Endo", "Lymph Endo", "Inj Endo", "DCT",
          "IT", "PT","TAL", "CDPC", "Des Thin Limb", "IC", "CT", "Uro", "SMC", "Fib", "Pericytes",
          "MyoF", "Podo", "Schwann", "Plasma", "Dendritic", "Monocytes", "Macrophages", "Baso", "T", "B",
          "Nat Killer"))

Idents(Samples.combined.CD10neg)<-"Annotation.Level.2"

# Supplemental Figure 2B
FeaturePlot(Samples.combined.CD10neg, features ="ADAMTS12", pt.size = 0.01, order=TRUE, label = T, label.size = 2)+ 
  scale_colour_gradient2(mid = 'lightgrey', high = 'red')+
  NoLegend()+ 
  theme_void()+
  theme(legend.position = "none",
        title = element_blank())
ggsave(filename = paste0(sample,"_ADAMTS12_red_wl.tiff"), width= 2 , height = 2)

# Supplemental Figure 2C ADAMTS12
DotPlot(Samples.combined.CD10neg, features = c("ADAMTS12"),dot.scale = 3)+ 
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"))+
  coord_flip()
ggsave(filename = paste0(sample, "_Supp Fig 1b ADAMTS12_l.svg"), width= 4, height = 1.2)

# Supplemental Figure 2C POSTN COL1A1 PDGFRB ACTA2
DotPlot(Samples.combined.CD10neg, features = c("POSTN", "COL1A1", "PDGFRB", "ACTA2"),dot.scale = 3)+ 
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), 
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),axis.ticks = element_line(size = 0.2),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"))+
  coord_flip()
ggsave(filename = paste0(sample, "_Supp Fig 1b A_wl.svg"), width= 3.9, height = 1.5)




