# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(Seurat)
library(writexl)
library(dplyr)
library(readxl)
library(mgsub)
library(ggpubr)
setwd("xxx") # Set your working directory here

# Read dataset
Samples.combined_PDGFRb <- readRDS("dataset/Human_PDGFRBpositive.rds")

# Abbreviate long names
Samples.combined_PDGFRb@meta.data[["Annotation.Level.2"]] <- 
  mgsub(Samples.combined_PDGFRb@meta.data[["Annotation.Level.2"]],
        c("Venular Endothelium", "Vascular Smooth Muscle Cells", "Thick Ascending Limb", "Proximal Tubule", "Prinicipal Cells","Intercalated Cells", "Injured Proximal Tubule","Glomerular Capillaries", "Arterioles"), 
        c("Venular Endo", "VSMC", "TAL", "PT", "PC", "IC", "Injured PT", "Glom Cap", "Arterioles"))
Samples.combined_PDGFRb$Annotation.Level.2 <- factor(Samples.combined_PDGFRb$Annotation.Level.2, 
                                                     levels = c("TAL", "PC", "IC", "Glom Cap", 
                                                                "Mesangial Cells", "Dendritic Cells", 
                                                                "VSMC", "Vasa Recta", "Arterioles", 
                                                                "Venular Endo", "Monocytes", "Pericytes", 
                                                                "Injured PT", "PT", "Fibroblasts", "Myofibroblasts"))
Idents(Samples.combined_PDGFRb)<-"Annotation.Level.2"

# Basic processing
Samples.combined_PDGFRb <- NormalizeData(Samples.combined_PDGFRb, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
Samples.combined_PDGFRb <- FindVariableFeatures(Samples.combined_PDGFRb, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
Samples.combined_PDGFRb <- ScaleData(Samples.combined_PDGFRb, verbose = FALSE, features = rownames(Samples.combined_PDGFRb))
Samples.combined_PDGFRb <- RunPCA(Samples.combined_PDGFRb, verbose = FALSE)

###################### Ext data Fig 1G - UMAP huPDGFRb+ ########################

for(i in  c(T, F)){
  if(i == T){text <- "ann"} 
  else{text <- "no_ann"}
  DimPlot(Samples.combined_PDGFRb, group.by = "Annotation.Level.2", label = i, label.size = 2)+ 
    NoLegend()+ 
    theme_void()+
    theme(legend.position = "none",
          title = element_blank())
  ggsave(filename = paste0("Output/Ext data Fig 1G - UMAP huPDGFRb+_", text,".tiff"), 
         width= 2.4, height = 2.3, dpi = 1000)
}
# Plot to extract legend
DimPlot(Samples.combined_PDGFRb, group.by = "Annotation.Level.2", label = F)+
  theme_void()+
  theme(title = element_blank(),
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.key.height = unit(0.1, "inches"))
ggsave(filename = "Output/Ext data Fig 1G - UMAP Legend.svg", 
       width= 2.3, height = 2.3, dpi = 1000)


################### Ext data Fig 1H - VlnPlots_combined ########################

goi <- c("PDGFRB", "PDGFRA", "GLI1")
myplots <- vector('list')
for(i in goi){
  myplots[[i]] <- 
    VlnPlot(Samples.combined_PDGFRb, flip=F, features = i,  group.by= "Annotation.Level.2",
            pt.size = 0.01, sort = )+ 
    ggtitle(i)+
    xlab("")+
    ylab("Expression Level")+
    theme_classic()+
    NoLegend()+
    theme(title = element_blank(),
          axis.text.x = element_text(size = 6, 
                                     family = "Arial", 
                                     angle = 0, 
                                     hjust =1,
                                     colour = "black"),
          axis.text.y = element_text(size = 8, 
                                     family = "Arial", 
                                     colour = "black"),
          axis.title = element_text(size =8, face= "plain", family = "Arial"),
          plot.title = element_text(face="plain", size = 10, hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    coord_flip()
}

# Plot
ggarrange(myplots[["PDGFRB"]]+ theme(axis.text.y = element_blank(),
                                     axis.text.x = element_blank(),
                                     plot.margin = margin(r = 1, l = 1)),
          myplots[["PDGFRA"]]+ theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.text.x = element_blank(),
                                     plot.margin = margin(r = 1, l = 1)),
          myplots[["GLI1"]]+ theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.text.x = element_blank(),
                                   plot.margin = margin(r = 1, l = 1)),
          nrow = 1)
ggsave(filename = "Output/Ext data Fig 1H - VlnPlots_combined.tiff", 
       width=6, height = 4, dpi = 1000)
