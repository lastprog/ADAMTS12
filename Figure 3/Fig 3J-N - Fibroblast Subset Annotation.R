library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)
library(readxl)
library(harmony)
windowsFonts("Arial" = windowsFont("Arial"))

#set working directory to the output folder of your data
setwd("xxx") # Set your working directory here
sample="ADAMTS12_Visium_MI"

#Load integrated Visium Seurat Object
Samples.combined <- readRDS(file = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/MI/Visium+Single WT and ADAMTS12KO 7 days after MI/Visium dataset/huMI ann/Scar cell type composition/Input/all_samples_ann_huMI.rds") #if its an rds file

#Annotate Genotype
Idents(Samples.combined)<-"Sample"
current.cluster.ids = levels(Idents(object = Samples.combined))
new.cluster.ids = c("ADAMTS12 KO", "ADAMTS12 KO", "WT", "WT")
Samples.combined$genotype<- plyr::mapvalues(x=Samples.combined$Sample, from = current.cluster.ids, to=new.cluster.ids)


#Prep for annotation according Seurat Tutorial: https://satijalab.org/seurat/articles/spatial_vignette.html
Samples.combined <- SCTransform(Samples.combined, assay = "Spatial", verbose = FALSE)
DefaultAssay(Samples.combined) <- "SCT"
Samples.combined <- RunPCA(Samples.combined, verbose = FALSE)
Samples.combined <- FindNeighbors(Samples.combined, dims = 1:30)
Samples.combined <- FindClusters(Samples.combined, verbose = FALSE)
Samples.combined <- RunUMAP(Samples.combined, dims = 1:30)


#Load Reference scRNA Dataset (Peisker et al. Nat. com) and prep according to seurat tutorial
ref <- readRDS(file = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 ?berexpression/2021 ADAMTS12 ?berexpression/Experiments/Bulk-Seq/Analysis/Fibroblast_integrated_filtered_processed_annotated.rds")
Idents(ref)<-"Annotated_Subsets"
DefaultAssay(ref)="RNA"

# SCTransform and then perform Fibroblast subset annotation
ref <- SCTransform(ref, ncells = 3000, verbose = FALSE)
anchors <- FindTransferAnchors(reference = ref, query = Samples.combined, normalization.method = "SCT", npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$Annotated_Subsets, prediction.assay = TRUE, weight.reduction = Samples.combined[["pca"]], dims = 1:50)
Samples.combined[["predictions"]] <- predictions.assay
DefaultAssay(Samples.combined) <- "predictions"


#extract Fibroblast subset probabilities and correct them for the overall likelihood of a fibroblast being called for each spot
Fib_Probability <- data.frame(Samples.combined$Fibroblast) # grab the intially called fibroblast probability for each spot
predictions.data <- t(data.frame(GetAssayData(Samples.combined[["predictions"]], slot = "data"))) # grab fibroblast subset prediction data
predictions.data <- as.data.frame(predictions.data)
Fib_df <- cbind(Fib_Probability, predictions.data)
Fib_df_Adj <- Fib_df*Fib_df$Samples.combined.Fibroblast
Fib_df_Adj <-t(Fib_df_Adj[,-1])
Fib_df_Adj_assay <- CreateAssayObject(counts = Fib_df_Adj)
Samples.combined[["predictionsadjusted"]]<-Fib_df_Adj_assay
DefaultAssay(Samples.combined) <- "predictionsadjusted"


# Create a meta-data column, where for each spot the most likely Fibroblast-Subset within the respective spot is the label (Max_Fibroblast Probability to Metadata to Plot Fibroblast Subset Annotation)
Fib_PA<-as.data.frame(t(as.data.frame(GetAssayData(Samples.combined, assay = "predictionsadjusted", slot = "counts")))) 
# Assign Fibroblast-Subset with Max Probability to each barcode
Fib_PA<-mutate(Fib_PA, Max_Fib = if_else( max == `ECM-Fibroblast`, "ECM-Fibroblast",if_else(
                                        max == `Atf3 IR Fibroblast`, "Atf3 IR Fibroblast", if_else(
                                        max == `Fibroblast 1`, "Fibroblast 1", if_else(
                                        max == `Fibroblast 2`, "Fibroblast 2", if_else(
                                        max == `Fibroblast 3`, "Fibroblast 3", if_else(
                                        max == `Interferon Fibroblast`, "Interferon Fibroblast", "Other"
                                        )))))))
Fib_PA$Max_Fib <- as.factor(Fib_PA$Max_Fib)
Fib_PA<-Fib_PA$Max_Fib
Samples.combined<-AddMetaData(Samples.combined, metadata = Fib_PA, col.name = "Fib_PA")

# Fig. 3j: Spatial Dimplot of most likely Fibroblast subset for each spot
Idents(Samples.combined)<-"Fib_PA"
SpatialDimPlot(Samples.combined, label = FALSE)&   theme(legend.position = "right",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))&
  ggtitle(" ")&
  scale_fill_manual(values = c( "Fibroblast 1" = "#AEA200",
                "Fibroblast 2" = "#64B200",
                "Fibroblast 3" = "#00C1A7",
                "Interferon Fibroblast" = "#00BADE",
                "ECM-Fibroblast" = "#F8766D",
                "Atf3 IR Fibroblast" = "#EF67EB"))&  
  guides(color = guide_legend(override.aes = list(size = 100)))
ggsave(filename=paste0(sample, "Fig3j_Dimplot_Fib_Annotation.jpeg"), width = 12, height = 3)
ggsave(filename=paste0(sample, "Fig3j_Dimplot_Fib_Annotation.svg"), width = 12, height = 3)


# Figure 3k: Barplot of Fibroblast subsets for each zone
library(tibble)
library(tidyr)
Idents(Samples.combined) <- "concat_leiden_0.5"
avgprob = as.data.frame(AverageExpression(Samples.combined, return.seurat = F, assays = "predictionsadjusted", slot = "counts",group.by = "concat_leiden_0.5"))
rownames(avgprob) <- c("2 Atf3 IR Fibroblast", "4 Fibroblast 3", "5 Fibroblast 2", "6 Fibroblast 1", "3 Interferon Fibroblast", "1 ECM-Fibroblast", "max")
avgprob <- avgprob[-7,] %>% rownames_to_column(var="Annotated_Subsets")
colnames(avgprob)<- c("celltype","HV #1", "HV #2", "IZ", "BZ", "RBZ", "Epi")
zone_order <- c("Epi", "HV #2","HV #1","RBZ", "BZ", "IZ")
cell_count_long <- gather(avgprob, key = "region", value = "AU")
cell_count_long <- cell_count_long[cell_count_long$region != "celltype",]
cell_count_long$celltype <- rep(avgprob$celltype, nrow(avgprob)/length(avgprob$celltype))
cell_count_long$AU<-as.numeric(cell_count_long$AU)

ggplot(cell_count_long, aes(x = factor(region, levels =  zone_order), 
        y =AU, fill=celltype)) + geom_bar(stat = "identity") + coord_flip()+
  xlab("Zone")+  theme_minimal()+
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5), 
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size =6, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))+
  scale_fill_manual(values = c( "6 Fibroblast 1" = "#AEA200",
                                "5 Fibroblast 2" = "#64B200",
                                "4 Fibroblast 3" = "#00C1A7",
                                "3 Interferon Fibroblast" = "#00BADE",
                                "1 ECM-Fibroblast" = "#F8766D",
                                "2 Atf3 IR Fibroblast" = "#EF67EB"))
ggsave(filename = paste0(sample,"Fig3k_CompositionperZone.svg"), width=3, height = 2)



# Plot Fibroblast Subset Composition for Ischemic Zone
Idents(Samples.combined) <- "concat_leiden_0.5"
VlnPlot(subset(Samples.combined,idents = 2), features=c(levels(ref$Annotated_Subsets)),ncol = 6,combine = TRUE, 
        split.by = "genotype",split.plot= T,cols = c('#00C5C0','#FA756C'), pt.size = 0) &
  geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0) &
  theme(
    axis.title = element_blank(),
    text = element_text(size=8, family = "Arial"),
    axis.text.x = element_text(size=6, family = "Arial", angle=45, hjust=1),
    axis.text.y = element_text(size=6, family = "Arial"),
    plot.title = element_text(size=8, family="Arial", face = "plain"),
    axis.line = element_line(size = 0.2), 
    axis.ticks = element_line(size = 0.2))
ggsave(filename = paste0(sample,"Fig3l_Overall_Vln.jpeg"), width=4 , height = 1.4)
ggsave(filename = paste0(sample,"Fig3l_Overall_Vln.svg"), width=4 , height = 1.4)

#SpatialFeaturePlots of Fib3 and IR Fib
SpatialFeaturePlot(Samples.combined, features="Fibroblast 3", ncol = 1) & scale_fill_gradient2(low="darkblue", mid="yellow", high="darkred", midpoint = 40, limits=c(0,80))
ggsave(filename = paste0(sample,"Fig3m_And_ExtendedData3j_Fibroblast 3_f.jpeg"), width=3 , height = 12)
ggsave(filename = paste0(sample,"Fig3m_And_ExtendedData3j_Fibroblast 3_f.svg"), width=3 , height = 12)

SpatialFeaturePlot(Samples.combined, features="Atf3 IR Fibroblast", ncol = 1) & scale_fill_gradient2(low="darkblue", mid = "yellow", high="indianred", midpoint = 20, limits=c(0,40))
ggsave(filename = paste0(sample,"Fig3n_And_ExtendedData3k_Atf3_IR_Fibroblast_f.jpeg"), width=3 , height = 12)
ggsave(filename = paste0(sample,"Fig3n_And_ExtendedData3k_Atf3_IR_Fibroblast_f.svg"), width=3 , height = 12)

#Save Annotated RDS
saveRDS(Samples.combined, file = paste0(sample, "_11_SeuratObject_Annotated.rds"))
