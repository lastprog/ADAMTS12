# Fig1e plotting ADAMTS12 expression in huCKD
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scico)
library(reticulate)
library(tibble)
library(tidyr)
windowsFonts("Arial" = windowsFont("Arial"))

setwd("xxx") # Set your working directory here

sample="huCKD_PDGFRb_"
Samples.combined <- readRDS("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/alt/huCKD_Nature_Decode_Myofibroblast_Dataset/Human_PDGFRBpositive.rds")
Idents(Samples.combined)<-"Annotation.Level.2"

Samples.combined <- NormalizeData(Samples.combined, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
Samples.combined <- FindVariableFeatures(Samples.combined, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined <- ScaleData(Samples.combined, verbose = FALSE, features = rownames(Samples.combined))
Samples.combined <- RunPCA(Samples.combined, verbose = FALSE)

FeaturePlot(Samples.combined, features ="ADAMTS12", pt.size = 0.8, order=TRUE, label = FALSE)+ scale_colour_gradient2(mid = 'lightgrey', high = 'red')
ggsave(filename = paste0(sample,"_ADAMTS12_red.jpeg"), width=7 , height = 7)
FeaturePlot(Samples.combined, features ="ADAMTS12", pt.size = 0.8, order=TRUE, label = TRUE)+ scale_colour_gradient2(mid = 'lightgrey', high = 'red')
ggsave(filename = paste0(sample,"_ADAMTS12_red_labeled.jpeg"), width=7 , height = 7)

VlnPlot(Samples.combined, features = "ADAMTS12",  group.by="Annotation.Level.2", flip = TRUE, pt.size = 0.01, sort = "Decreasing") &NoLegend() & coord_flip() &
  theme(plot.title = element_text(face="plain", size = 10, family="Arial", hjust = 0.5), axis.title = element_text(face = "plain", size = 8, family="Arial"),
        axis.text.x = element_text(angle = 45, hjust = 1, size =8, face= "plain", family= "Arial"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"_ADAMTS12_Vln.jpeg"), width=2.8 , height = 2.8)
ggsave(filename = paste0(sample,"_ADAMTS12_Vln.jpeg"), width=2.8 , height = 2.8)
