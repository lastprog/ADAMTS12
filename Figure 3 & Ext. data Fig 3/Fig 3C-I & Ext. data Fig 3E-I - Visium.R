# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(ggplot2)
library(tidyr)
library(ggpubr)
library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(progeny)
library(dorothea)
library(ggrepel)
library(gprofiler2)
library(mgsub)
library(tibble)
setwd("xxx") # Set your working directory here
source("Source/my_vulcano_nice.R")
source("Source/colorado_bold text ggplot2.R")

data <- readRDS("Input/all_samples.rds")

################# Preprocessing, Normalization and Integration ################# 

data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, verbose = F)
data <- RunPCA(data, pc.genes = data@var.genes, npcs = 20, verbose = FALSE)

# Harmony Integration
options(repr.plot.height = 2.5, repr.plot.width = 6)
data <- data %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

# Removing malat1 and mt Genes
data <- data[-which(grepl("^mt-",rownames(data))),]
data <- data[-which(grepl("^Malat1",rownames(data))),]

# Assign Genotype
data@meta.data$genotype <- ifelse(str_starts(data@meta.data$orig.ident, "KO"), "KO", "WT")

# Name clusters
data@meta.data['concat_leiden_0.5'] <- mgsub(data@meta.data['concat_leiden_0.5'],
                                             pattern = c("0","1","2","3","4","5"), 
                                             replacement = c("HM #1", "HM #2", "IZ", "BZ", "RBZ", "Epi"))
data@meta.data['concat_leiden_0.5'] <- factor(data@meta.data['concat_leiden_0.5'][,1],
                                              levels = c("IZ", "BZ", "RBZ", "HM #1", "HM #2", "Epi"))
Idents(data) <- data@meta.data['concat_leiden_0.5']

# Split object by sample
sam_split <- SplitObject(data, split.by = "orig.ident")
samples <- unique(data@meta.data[["orig.ident"]])
images <- names(sam_split[["KO.1"]]@images)

###############################################################################

################# Ext. Data Fig 3E - Dotplot Marker Genes ######################
# Find Marker Genes
marker_genes <- FindAllMarkers(data)

# Extract top marker genes to plot
top5 <- marker_genes %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% pull(gene)
top5 <- unique(top5)

t <- marker_genes %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Plot Marker Genes
DotPlot(data, features = top5, dot.scale = 3)+
  scale_colour_gradient2(
    low = "darkblue",
    mid = "whitesmoke",
    high = "indianred",
    midpoint = 0)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size =8, angle = 90, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.2, 'cm'),
        legend.position = "bottom")
ggsave("Output/Ext. data Fig 3D - Dotplot Markergenes per cluster_on_logFC.svg", width = 4, height = 2)

################### Ext. data Fig 3F - cell type annotation ###################

for(i in seq(1, length(images), 1)){
  sam <- samples[i]
  img <- images[i]
  Idents(sam_split[[i]]) <- sam_split[[i]]@meta.data['max_like_cell']
  SpatialDimPlot(sam_split[[i]], label = FALSE, images = img) + 
    theme(legend.position = "none",
          plot.title = element_blank())+
    scale_fill_manual(
      values = c("Fibroblast" = "#F8766D",
                 "Cycling.cells" = "#DB8E00",
                 "Myeloid" = "#AEA200",
                 "Cardiomyocyte" = "#64B200",
                 "Endothelial" = "#00BD5C",
                 "vSCMs" = "#00C1A7",
                 "Pericyte" = "#00BADE",
                 "Mast" = "#00A6FF",
                 "Adipocyte" = "#B385FF",
                 "Neuronal" ="#EF67EB",
                 "Lymphoid" ="#FF63B6"))
  ggsave(filename = paste0("Output/Ext. data Fig 3_celltypes_",i,".svg"))
}

# Plot for extracting legend
SpatialDimPlot(sam_split[["KO.1"]], label = FALSE, images = images[1])+ 
  theme_minimal()+
  theme(legend.position = "right",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))+
  scale_fill_manual(
    values = c("Fibroblast" = "#F8766D",
               "Cycling.cells" = "#DB8E00",
               "Myeloid" = "#AEA200",
               "Cardiomyocyte" = "#64B200",
               "Endothelial" = "#00BD5C",
               "vSCMs" = "#00C1A7",
               "Pericyte" = "#00BADE",
               "Mast" = "#00A6FF",
               "Adipocyte" = "#B385FF",
               "Neuronal" ="#EF67EB",
               "Lymphoid" ="#FF63B6"))
ggsave(filename = "Output/Ext. data Fig 3 - celltypes_legend.svg")

######################## Fig 3C - Spatial niches ###############################
for(i in seq(1, length(images), 1)){
  sam <- samples[i]
  img <- images[i]
  SpatialDimPlot(sam_split[[i]], label = FALSE, images = img) + 
    theme(legend.position = "none")+
    ggtitle(" ")+
    scale_fill_manual(values = c("HM #1" = "#00BADE",
                                 "HM #2" = "#64B200",
                                 "IZ" = "#F8766D",
                                 "BZ" = "#AEA200",
                                 "RBZ" = "#00C1A7",
                                 "Epi" = "#EF67EB"))
  ggsave(filename = paste0("Output/Fig 3C - spatial niches_", sam,".svg"))
}

# Plot for extracting legend
SpatialDimPlot(sam_split[[1]], label = FALSE)+ 
  ggtitle(" ")+
  scale_fill_manual(
    values = c("HV #1" = "#00BADE",
               "HV #2" = "#64B200",
               "IZ" = "#F8766D",
               "BZ" = "#AEA200",
               "RBZ" = "#00C1A7",
               "Epi" = "#EF67EB"))+
  guides(color = guide_legend(override.aes = list(size = 100)))+
  theme(legend.position = "right",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))
ggsave(filename = "Output/Fig 3C - spatial niches_legend.svg")

################### Fig 3 D - Spatial ADAMTS12 expression ######################

for(i in seq(1,length(images), 1)){
  sam <- samples[i]
  img <- images[i]
  SpatialFeaturePlot(sam_split[[i]], features = "Adamts12", max.cutoff = 1.5, images = img)+
    theme(legend.position = "right",
          legend.key.height= unit(0.3, 'cm'),
          legend.key.width= unit(0.2, 'cm'),
          legend.title = element_text(size=8),
          legend.text = element_text(size=6))
  ggsave(filename = paste0("Output/Fig 3D - ADAMTS12_exp_",sam,".svg"))
}

################# Fig 3 E - Stacked barplot of spatial ADAMTS12 ################

# Extract only WT samples
wt_comb <- merge(sam_split[["WT.1"]], sam_split[["WT.2"]])
wt_comb@meta.data[["genotype"]] <- "WT"

# Extract ADAMTS12 expression, Celltype and Cluster from Seurat object
adamts12_wt <- FetchData(wt_comb, vars = "Adamts12")
adamts12_wt$celltype <- wt_comb@meta.data[["max_like_cell"]]
adamts12_wt$cluster <- wt_comb@meta.data[["concat_leiden_0.5"]]

# Sum ADAMTS12 Expression by Celltype and Cluster
adamts12_wt <- aggregate(adamts12_wt$Adamts12  ~ adamts12_wt$cluster + adamts12_wt$celltype,data = adamts12_wt, sum)

# Devide expression by 2 to get the average expression per sample
adamts12_wt$per_sample <- adamts12_wt$`adamts12_wt$Adamts12`/2

# Plot stacked barplot
adamts12_wt$`adamts12_wt$cluster` <- factor(adamts12_wt$`adamts12_wt$cluster`,
       levels = c("IZ", "BZ", "RBZ", "HM #1", "HM #2", "Epi"))

ggplot(adamts12_wt, aes(x = adamts12_wt$`adamts12_wt$cluster`,
                        y = adamts12_wt$`adamts12_wt$Adamts12`, 
                        fill = reorder(adamts12_wt$`adamts12_wt$celltype`, adamts12_wt$`adamts12_wt$Adamts12`)))+
  geom_bar(stat = "identity")+
  ylab("total ADAMTS12\n per sample")+
  scale_fill_manual(
    values = c(
      "Fibroblast" = "#F8766D",
      "Cycling.cells" = "#DB8E00",
      "Myeloid" = "#AEA200",
      "Cardiomyocyte" = "#64B200",
      "Endothelial" = "#00BD5C",
      "vSCMs" = "#00C1A7",
      "Pericyte" = "#00BADE",
      "Mast" = "#00A6FF",
      "Adipocyte" = "#B385FF",
      "Neuronal" ="#EF67EB",
      "Lymphoid" ="#FF63B6"))+
  theme_minimal()+
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5), 
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size =8, angle = 45, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.text = element_text(size =8, face= "plain", family = "Arial"), 
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))
ggsave("Output/Fig. 3E stacked barplot_ADAMTS12 Expression_sum.svg", height = 1.3, width = 2.9)

############## Fig. 3F & Ext. data Fig 3G - cell type composition ##############

# Celltypes and samples in the dataset
samples <- names(sam_split)
zones <- unique(data@meta.data[["concat_leiden_0.5"]])
ct <- c("Fibroblast", "Cardiomyocyte", "Endothelial", "Myeloid", "Neuronal", 
        "Pericyte", "Mast", "Lymphoid", "Cycling.cells", "vSMCs", "Adipocyte")

# Assigning spatial niches as Idents
for(i in samples){
  Idents(sam_split[[i]]) <- sam_split[[i]]@meta.data[["concat_leiden_0.5"]]
}

# Extracting the cell type probabilities for each sample for each zone
ls <- list()
for(z in zones){
  ls[[z]] <- list()
  for(g in samples){
    ls[[z]][[g]] <- list()
    sd_subs <- subset(sam_split[[g]], idents = z) # take only one zone and 1 genotype
    for(i in ct){
      tmp <- sd_subs@meta.data[[i]]
      ls[[z]][[g]][[i]] <- list()
      ls[[z]][[g]][[i]] <- mean(tmp) # Calculate already the mean for all Probabilities
    }
  }
}

# Scale probability scores so that they sum to 100%
zones_prob <- list()
for(z in zones){
  tmp <- lapply(ls[[z]], as.numeric)
  tmp <- as.data.frame(do.call(rbind, tmp))
  
  # Scale all probabilites to 100%
  tmp$sum <- rowSums(tmp)
  tmp$scale <- 100/tmp$sum
  tmp <- tmp*tmp$scale
  colnames(tmp) <- ct
  zones_prob[[z]] <- t(tmp[,-c(12, 13)])
}

########################### Fig 3F - log2fc of IZ ##############################

zone_df <- as.data.frame(zones_prob[["IZ"]])
means_geno <- data.frame(KO = (zone_df$KO.1+zone_df$KO.2)/2,
                         WT = (zone_df$WT.2+zone_df$WT.1)/2)
rownames(means_geno) <- rownames(zone_df)

# Calculating log2 fc
means_geno$fc <- (means_geno$KO-means_geno$WT)/means_geno$WT
means_geno$log2fc <- log2(means_geno$fc+1)
means_geno <- rownames_to_column(means_geno)

# Plotting Fig 3F
ggplot(means_geno,aes(x = reorder(means_geno$rowname, means_geno$log2fc), y =means_geno$log2fc)) + 
  geom_bar(aes(fill = means_geno$log2fc), stat = "identity") + 
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  ggtitle("cell type composition\nischemic zone")+
  ylab("log2fc")+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0)+
  scale_y_continuous(breaks = c(-0.5, 0, 0.5,1))+
  theme_minimal() +
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5), 
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.text.x = element_text(size =6, angle = 0, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.position = "none",
        legend.key.width= unit(0.4, 'cm'),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank())+
  xlab("")+
  coord_flip()
ggsave(filename = "Output/Fig 3F - celltypes_log2fc_barplot.svg", width = 1.5, height = 2.1, units = "in")


################# Ext. Data Fig 3G - cell type composition by zone ##############

p <- list()
for(z in zones){
  zone_df <- as.data.frame(zones_prob[[z]])
  means_geno <- data.frame(KO = (zone_df$KO.1+zone_df$KO.2)/2,
                           WT = (zone_df$WT.2+zone_df$WT.1)/2)
  means_geno_long <- gather(means_geno, key = "sample", value = "percent")
  means_geno_long$celltype <- rep(ct, 2)
  means_geno_long$sample <- factor(means_geno_long$sample, levels = c("WT", "KO"))
  
  p[[z]] <- ggplot(means_geno_long, aes_string(y = means_geno_long$percent,
                                               x = means_geno_long$sample,
                                               fill = reorder(means_geno_long$celltype, means_geno_long$percent)))+
    geom_bar(stat = "identity", position = "stack")+
    ggtitle(z)+
    scale_fill_manual(
      values = c("Fibroblast" = "#F8766D",
                 "Cycling.cells" = "#DB8E00",
                 "Myeloid" = "#AEA200",
                 "Cardiomyocyte" = "#64B200",
                 "Endothelial" = "#00BD5C",
                 "vSMCs" = "#00C1A7",
                 "Pericyte" = "#00BADE",
                 "Mast" = "#00A6FF",
                 "Adipocyte" = "#B385FF",
                 "Neuronal" ="#EF67EB",
                 "Lymphoid" ="#FF63B6"))+
    theme_minimal()+
    theme(plot.title = element_text(size = 8, face="plain", hjust = 0.5),
          axis.title = element_blank(),
          axis.text.x = element_text(size =8, angle = 0, hjust = 1, face= "plain", family= "Arial", color = "black"),
          axis.text.y = element_text(size =6, face= "plain", family = "Arial", color = "black"),
          legend.text = element_text(size =8, face= "plain", family = "Arial"),
          legend.title = element_blank(),
          legend.key.width= unit(0.2, 'cm'),
          legend.key.height= unit(0.1, 'cm'))
}

# Plot
ggarrange(p[["IZ"]], p[["BZ"]],p[["RBZ"]],  p[["HM #1"]], p[["HM #2"]], p[["Epi"]],
          common.legend = T, legend = "right")
ggsave("Output/Ext. data Fig 3G - Cell_composition different zones.svg",  height = 2.5, width = 3.5)

############## Ext. Data Fig 3G - cell type composition all niches ############# 
means <- list()
for(g in samples){
  ls <- list()
  for(i in ct){
    tmp <- sam_split[[g]]@meta.data[[i]]
    ls[[i]] <- tmp
  }
  tmp <- as.data.frame(ls)
  # Scale all probabilites to 100%
  tmp$sum <- rowSums(tmp)
  tmp$scale <- 100/tmp$sum
  tmp <- tmp*tmp[,"scale"]
  tmp <- tmp[, -c(13,12)]
  
  # Mean of cell type probability
  means[[g]] <- colMeans(tmp)
}
means <- as.data.frame(means)

# Plot
means_sample_long <- gather(means, key = "sample", value = "percent")
means_sample_long$celltype <- rep(rownames(means), 4)
means_sample_long$sample <- factor(means_sample_long$sample, levels = c("WT.1", "WT.2", "KO.1", "KO.2"))

ggplot(means_sample_long, aes(y = means_sample_long$percent, 
                              x = means_sample_long$sample, 
                              fill = reorder(means_sample_long$celltype,means_sample_long$percent)))+
  geom_bar(stat = "identity", position = "stack")+
  ggtitle("Celltype composition\n(total)")+
  ylab("Celltype composition [%]")+
  xlab("")+
  scale_fill_manual(
    values = c("Fibroblast" = "#F8766D",
               "Cycling.cells" = "#DB8E00",
               "Myeloid" = "#AEA200",
               "Cardiomyocyte" = "#64B200",
               "Endothelial" = "#00BD5C",
               "vSMCs" = "#00C1A7",
               "Pericyte" = "#00BADE",
               "Mast" = "#00A6FF",
               "Adipocyte" = "#B385FF",
               "Neuronal" ="#EF67EB",
               "Lymphoid" ="#FF63B6"))+
  theme_minimal()+
  labs(fill = "Cell Types")+
  theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5), 
        axis.title = element_text(size = 8, face = "plain", family="Arial"),
        axis.text.x = element_text(size =8, angle = 0, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))
ggsave("Output/Supp Fig 3G -  cell_composition_whole heart.svg", height = 2, width = 2.8)

############################## DEGs for ischemic zone ##########################
# Subsetting only the ischemic zone
Idents(data)<- data@meta.data$genotype
iz_only <- data[,data@meta.data[["concat_leiden_0.5"]] == "IZ"]

# Calculate DEGs
degs_iz <- FindAllMarkers(iz_only, only.pos = T, min.pct = 0.25, logfc.threshold = 0.0, test.use = "MAST")

######################### Supp table 2 - DEGs of IZ ############################
write.csv(degs_iz,"Output/Supp. table 2 - Visium DEGs IZ.csv")

###################### Fig 3G - Vulcano Plot of DEGs in IZ #####################
degs_iz <- read.csv("Output/DEGs IZ.csv")
degs_iz$avg_log2FC <- ifelse(degs_iz$cluster == "WT", degs_iz$avg_log2FC*-1, degs_iz$avg_log2FC)

# Prepare labels
top10 <- degs_iz %>% top_n(n = 10, wt = p_val_adj) %>% pull(`gene`)
bottom10 <- degs_iz  %>% top_n(n = -10, wt = p_val_adj) %>% pull(`gene`)
labels <- c(top10,bottom10)

# Vulcano Plot for annotation and for plotting
for(i in c(T, F)){
  if(i == T){
    txt <- "with_label"
  }else{txt <- ""}
  volcano_nice(df = degs_iz, hAss = 0.05, FCIndex = 3,pValIndex = 6, 
               IDIndex = 8, nlabels = 10, label = i, manual_labels = c(labels, "Adamts12"))+
    ggtitle("DEGs")+
    labs(x="logFC", y="-log(P-Val.)")+
    scale_x_continuous(breaks = c(-1.5, 0, 1.5))+
    theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
          axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
          axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
          axis.title = element_text(size =8, face= "plain", family = "Arial"),
          panel.grid.major = element_line(colour = "grey", size = 0.05),
          panel.grid.minor = element_line(colour = "grey", size = 0.05))
  ggsave(filename = paste0("Output/Fig 3G - Vulcano DEGs_IZ_", txt,".svg"), width = 2.1, height =2, units = "in", dpi = 1000)
}

############################### GO Terms ######################################
# Extracting only genes which are significantly upregulated in KO
ko_deg <- degs_iz[degs_iz$p_val_adj < 0.05 & grepl("KO", degs_iz$cluster),]

# Calculate GO Terms for sig. upregulated genes in KO (ischemic zone)
ko_up <- gost(query = ko_deg$gene, 
              organism = "mmusculus", ordered_query = T, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              domain_scope = "annotated", custom_bg = NULL, 
              numeric_ns = "", sources = NULL, as_short_link = FALSE) 
res_up <- ko_up$result

######################### Fig 3H - Reactome Terms #############################

# Select only Reactome Terms with termsize > 10
go_reac_up <- res_up[res_up$source == "REAC" & res_up$term_size > 10, ]

# Prepare for plotting
go_reac_up <- go_reac_up[order(go_reac_up$p_value),]
top10_up <- head(go_reac_up, 10)
top10_up$mlog10_pvalue <- -log10(top10_up$p_value)
top10_up$term_name <- mgsub(top10_up$term_name, c("Striated Muscle Contraction",
                                                  "Cell-extracellular matrix interactions"),
                            c("Str Muscle Contraction",
                              "Cell-EZM interactions"))

# Plot Top 10 Reactome Terms upregulated in KO(ischemic zone)
ggplot(top10_up,
       aes(x = reorder(factor(term_name), top10_up$mlog10_pvalue), y = top10_up$mlog10_pvalue))+ 
  geom_bar(aes(fill = top10_up$mlog10_pvalue), stat = "identity")+
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = median(top10_up$precision))+
  scale_y_continuous(breaks = c(0, 1, 2))+
  labs(title = "Reactome Terms")+
  xlab("")+
  ylab("-log10(p-value)")+
  theme_minimal()+
  theme(axis.title = element_text(size =8, face= "plain", family = "Arial"),
        plot.title = element_text(size = 9, face="plain",  hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial",
                                   color = "black"),
        axis.text.y = element_text(size =8,
                                   family = "Arial",
                                   color = "black"),
        legend.position = "none",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"), 
        panel.grid.minor = element_blank(),  
        legend.key.width= unit(0.4, 'cm')) +
  coord_flip()
ggsave(filename = "Output/Fig 3H - Reactome_Terms_barplot.svg", width = 1.8, height = 2, units = "in")

####################### Ext. data. Fig 3H - GO:BP Terms ########################

# Select only GO:BP Terms with term size <2000
go_bp_up <- res_up[res_up$source == "GO:BP" & res_up$term_size < 2000, ]

# Prepare for plotting
go_bp_up <- go_bp_up[order(go_bp_up$p_value),]
top10_up <- head(go_bp_up, 10)
top10_up$mlog10_pvalue <- -log10(top10_up$p_value)
top10_up$term_name <- mgsub(top10_up$term_name, c("circulatory system development",
                                                  "establishment of protein localization",
                                                  "anatomical structure formation involved in morphogenesis",
                                                  "muscle structure development",
                                                  "cellular catabolic process",
                                                  "regulation of protein modification process"),
                            c("circulatory sys dev",
                              "est protein localization",
                              "anatomical struc from",
                              "muscle structure dev",
                              "cell catabolic proc",
                              "reg of prot mod process"))

# Plot Top 10 GO:BP Terms upregulated in KO vs WT (ischemic zone)
ggplot(top10_up,
       aes(x = reorder(factor(term_name), top10_up$mlog10_pvalue), y = top10_up$mlog10_pvalue))+ 
  geom_bar(aes(fill = top10_up$mlog10_pvalue), stat = "identity")+
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = median(top10_up$precision))+ 
  scale_y_continuous(breaks = c(0, 20, 40))+
  labs(title = "Top 10 GO:BP")+
  xlab("")+
  ylab("-log10(p-value)")+
  theme_minimal()+
  theme(axis.title = element_text(size =8, face= "plain", family = "Arial"),
        plot.title = element_text(size = 9, face="plain",  hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial",
                                   color = "black"),
        axis.text.y = element_text(size =8,
                                   family = "Arial",
                                   color = "black",
                                   face = colorado(reorder(factor(top10_up$term_name), top10_up$mlog10_pvalue, decreasing = F), 
                                                   c("circulatory sys dev",
                                                     "regulation of cell death",
                                                     "muscle structure dev"))),
        legend.position = "none",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.minor = element_blank(),  
        legend.key.width= unit(0.4, 'cm')) +
  coord_flip()
ggsave(filename = "Output/Ext Fig 3H - Top10_GO_BP-Terms_barplot.svg", width = 2, height = 2.3, units = "in")

########################## Fig. 3I - Progeny IZ ###############################

# Capitalize gene names
first_letters <- toupper(substr(degs_iz$gene, 1, 1))
name <- paste0(first_letters, substr(degs_iz$gene, 2, nchar(degs_iz$gene)))
degs_iz$names <- name

# Prepare datamatrix
progeny_matrix <- degs_iz %>% 
  dplyr::select(names, avg_log2FC) %>% 
  dplyr::filter(!is.na(avg_log2FC)) %>% 
  column_to_rownames(var = "names")%>%
  as.matrix()

# Run Progeny
PathwayActivity_counts <- progeny(progeny_matrix, 
                                  scale = T, 
                                  organism = "Mouse", 
                                  top = 100, perm = 10000, 
                                  z_scores = TRUE) %>% t
colnames(PathwayActivity_counts) <- "NES"

# Arrange Progeny output
PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_counts) %>% 
  rownames_to_column(var="Pathway") %>%
  dplyr::arrange(NES) %>% mutate(Pathway=factor(Pathway))

PathwayActivity_zscore_df <- rownames_to_column(PathwayActivity_zscore_df)
rownames(PathwayActivity_zscore_df) <- PathwayActivity_zscore_df$rowname

# Plot Progeny in Barplot
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, -NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") + 
  labs(title = "Progeny") +
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  scale_y_continuous(breaks = c(-2,0, 2))+
  theme_minimal() +
  theme(plot.title = element_text(face="plain", size = 9, hjust = 0.5), 
        axis.title = element_text(face = "plain", size = 8, family="Arial"),
        axis.text.x = element_text(hjust = 1, size =6, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial",color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.minor = element_blank(),
        legend.key.width= unit(0.4, 'cm'),
        legend.position = "none") +
  xlab("")+
  coord_flip()
ggsave(filename = "Output/Fig 3I - Visium ischemic zone Progeny.svg",  width = 1.6, height = 2.2)

####################### Ext data Fig. 3H - Dorothea IZ #########################

# Prepare Dataframe - only use regulons of confidence A-C
data(dorothea_mm, package="dorothea")
regulons <- dorothea_mm %>% dplyr::filter(confidence %in% c("A", "B", "C"))

current_df_matrix <- degs_iz %>% 
  dplyr::select(gene,avg_log2FC) %>% 
  dplyr::filter(!is.na(avg_log2FC)) %>%
  column_to_rownames(var="gene")%>%
  as.matrix()

# Run Dorothea
tf_activities_stat = dorothea::run_viper(current_df_matrix, regulons, 
                                         options = list(minsize=30, eset.filter=FALSE,
                                                        cores=1, verbose=FALSE, nes=TRUE))
# Reformat results into dataframe
tf_activities_stat_top25<- tf_activities_stat%>%
  as.data.frame() %>%
  rownames_to_column(var="gene") %>%
  dplyr::rename(NES="avg_log2FC") %>%
  dplyr::top_n(20, wt=abs(NES)) %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(GeneID=factor(gene))

# Plot Dorothea
ggplot(tf_activities_stat_top25, aes(y=reorder(GeneID, -NES), x=NES))+ 
  geom_bar(aes(fill=NES), stat="identity")+
  ggtitle("Dorothea")+
  geom_vline(xintercept = 0,
             color = "black", size = 0.3)+
  xlab("NES")+
  scale_fill_gradient2(low="darkblue", high="indianred", mid="whitesmoke", midpoint=0)+
  scale_x_continuous(breaks = c(-4,-2,0, 2))+
  theme_minimal()+
  theme(plot.title =  element_text(size = 9, hjust = 0.5), axis.title=element_text(size=12),
        axis.text.x=element_text(angle=0, hjust=1, size=6, color = "black"),
        axis.text.y = element_text(size=8, color ="black"),
        axis.title.x = element_text(size=8),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
ggsave(filename = "Output/Ext. Fig. 3I - Dorothea IZ.svg", width = 1.6, height = 2.8, bg = "white")



