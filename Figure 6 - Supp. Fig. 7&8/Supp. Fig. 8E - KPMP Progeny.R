# Fig1e plotting ADAMTS12 expression in huCKD
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(scico)
library(reticulate)
library(tibble)
library(tidyr)
library(progeny)
library(pheatmap)
library(readr)
windowsFonts("Arial" = windowsFont("Arial"))

setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/Revision/PROGENyDorothea/")

# Read dataset
sample="KPMP_"
Samples.combined <- LoadH5Seurat("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/KPMP/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")
Idents(Samples.combined)<-"subclass.l2"
Samples.combined$subclass.l2.disease<-paste(Samples.combined$subclass.l2, Samples.combined$diseasetype, sep = "_")
Idents(Samples.combined)<-"subclass.l2.disease"

#PROGENy Analysis
Samples.combined <- progeny(Samples.combined, scale=FALSE, organism="Human", top=500, perm=1, return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
Samples.combined <- Seurat::ScaleData(Samples.combined, assay = "progeny") 

#Subset interstitial cells for better visualization
Idents(Samples.combined)<-"subclass.l1"
Samples.combined.Interstitial<-subset(Samples.combined, idents = "Interstitial")
Idents(Samples.combined.Interstitial)<-"subclass.l2.disease"

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <- as.data.frame(t(GetAssayData(Samples.combined.Interstitial, slot = "scale.data", assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## We match Progeny scores with the cell clusters.
CellsClusters <- data.frame(Cell = names(Idents(Samples.combined.Interstitial)), 
                            CellType = as.character(Idents(Samples.combined.Interstitial)),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% group_by(Pathway, CellType) %>%  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

write_csv(summarized_progeny_scores_df, file = paste0(sample, "03_PROGENy_Interstitial_disease.csv"))
progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

# Supp. Fig. 8E
pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14,fontsize_row = 10, 
         color=myColor, breaks = progenyBreaks, main = "PROGENy (500)", angle_col = 45,
         treeheight_col = 0,  border_color = NA, filename = paste0(sample, "04_PROGENY_L2_disease.pdf"), width = 10, height=6)

