# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(Seurat)
library(ggplot2)
library(dplyr)
library(mgsub)
library(tidyr)
setwd("xxx") # Set your working directory here

# Specifiy samples and zones of input files
samples <- list.files("Input/",pattern = "*.rds",recursive = TRUE, full.names = TRUE)
n <- c("control", "fibrotic", "ischemic", "border", "remote")

# Read human MI dataset and plot SpatialDimPlot 
for(i in 1:length(samples)){
  # Load dataset
  vis_data <- readRDS(samples[i])
  Idents(object = vis_data) <- vis_data@meta.data$CellTypes
  
  # Plot Cell type annotation
  SpatialDimPlot(vis_data)+
    theme(plot.title = element_blank(),
          legend.position = "none")+
    scale_fill_manual(
      values = c(
        "Fibroblast" = "#F8766D",
        "Cycling.cells" = "#DB8E00",
        "Myeloid" = "#AEA200",
        "Cardiomyocyte" = "#64B200",
        "Endothelial" = "#00BD5C",
        "vSMCs" = "#00C1A7",
        "Pericyte" = "#00BADE",
        "Mast" = "#00A6FF",
        "Adipocyte" = "#B385FF",
        "Neuronal" ="#EF67EB",
        "Lymphoid" ="#FF63B6"))
  ggsave(filename =  paste0("Output/SpatialDimPlot_celltypes_", n[i], ".svg"), width = 1, height = 1)
  
  # Plot Spatial Feature Plot
  SpatialFeaturePlot(vis_data, features = "ADAMTS12",min.cutoff = 0,max.cutoff = 1.5)+
    theme(legend.position = "none")
  ggsave(filename =  paste0("Output/SpatialFeaturePlot_ADAMTS12_", n[i], ".svg"), width = 1, height = 1)
}

# Plot to extract Legend
SpatialFeaturePlot(vis_data, features = "ADAMTS12",min.cutoff = 0,max.cutoff = 1.5)+
  theme(legend.text = element_text(size =6, face= "plain", family = "Arial", angle = 270),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))
ggsave("Output/huMI_adamts12_legend.svg")





