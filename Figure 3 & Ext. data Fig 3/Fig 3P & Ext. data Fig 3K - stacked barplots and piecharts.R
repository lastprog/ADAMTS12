# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(ggplot2)
library(dplyr)
library(tidyr)
library(mgsub)
setwd("xxx") # Set your working directory here

# Reading dataframe of ADAMTS12 expression across zones and most probable cell types
whole_data <- read.csv("Input/compiled_all_MI_samples.csv")

# Define zone names
zones <- c("Control", "Ischemic", "Border", "Remote", "Fibrotic")

# Rename zones
whole_data$sample_category <- mgsub(whole_data$sample_category, c("CTRL", "IZ", "BZ", "RZ", "FZ"), 
                                    zones)

##################### Fig 3P - Spatial ADAMTS12 exp. in huMI ###################

# Aggregate samples (sum)
adamts12_exp_by_ct <- aggregate(whole_data$ADAMTS12 ~ whole_data$sample_category + whole_data$cell_type_original, 
                                data = whole_data, sum)

# Plot stacked barplot
level_order <- rev(zones)
ggplot(adamts12_exp_by_ct, aes(x = factor(adamts12_exp_by_ct$`whole_data$sample_category`, levels = level_order),
                               y = adamts12_exp_by_ct$`whole_data$ADAMTS12`/1000, 
                               fill = reorder(adamts12_exp_by_ct$`whole_data$cell_type_original`, adamts12_exp_by_ct$`whole_data$ADAMTS12`)))+
  geom_bar(stat = "identity")+
  ylab("total ADAMTS12\n[in 1000]")+
  scale_y_continuous(breaks = c(0, 2, 4))+
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
        axis.title.y = element_blank(),
        axis.text.x = element_text(size =6, angle = 0, hjust = 1, face= "plain", family= "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_blank(),
        legend.key.width= unit(0.2, 'cm'),
        legend.key.height= unit(0.1, 'cm'))+
  coord_flip()
ggsave("Output/huMI ds ADAMTS12 Expression_sum.svg", height = 1.8, width = 2.3)


################ Ext. Data Fig 3 - Cell type composition per zone ##############

# Count cell type occurence
ct_count <- list()
for(i in zones){ 
  subs <- whole_data[whole_data$sample_category == i,]
  ct_count[[i]] <- table(subs$cell_type_original)
}

# Extract unique celltype names
cell_types <- unique(names(c(ct_count$Control, ct_count$Ischemic, ct_count$Border, ct_count$Remote, ct_count$Remote, ct_count$Fibrotic)))

# Calculate percentage of each cell type by zone
df_cc <- data.frame(celltype = cell_types)
for(i in zones){
  tmp <- as.numeric(ct_count[[i]][match(df_cc$celltype, names(ct_count[[i]]))])
  tmp <- tmp/sum(tmp, na.rm = T) 
  df_cc <- cbind(df_cc, tmp)
}
colnames(df_cc) <- c("celltype", zones)

# Prepare dataframe for plotting
cell_count_long <- gather(df_cc, key = "sample", value = "percent")
cell_count_long <- cell_count_long[cell_count_long$sample != "celltype",]
cell_count_long$celltype <- rep(df_cc$celltype, nrow(df_cc)/length(df_cc$celltype))
cell_count_long <- na.omit(cell_count_long)
cell_count_long$percent <- as.numeric(cell_count_long$percent)*100

# Plot stacked barplot of cell type composition
ggplot(cell_count_long, aes(x = factor(cell_count_long$sample, levels = zones),
                            y =cell_count_long$percent, 
                            fill = reorder(cell_count_long$celltype, cell_count_long$percent)))+
  geom_bar(stat = "identity")+
  ylab("cell type composition [%]")+
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
ggsave("Output/Ext data Fig 3L - cell type composition_huMI.svg", height = 2.4, width = 3)

################### Fig 3P - Cell type composition pie charts ##################

for(i in zones){
  tmp <- cell_count_long[cell_count_long$sample == i,]
  
  ggplot(tmp, aes(x = factor(tmp$sample, levels = zones),
                  y =tmp$percent, 
                  fill = reorder(tmp$celltype, tmp$percent)))+
    geom_bar(stat = "identity")+
    coord_polar("y", start = 0)+
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
        "Lymphoid" ="#FF63B6"))+
    theme_void()+
    theme(plot.title = element_text(size = 9, face="plain", hjust = 0.5), 
          axis.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none")
  ggsave(paste0("Output/cell_type_comp pie_charts_", i, ".svg"), width = 0.2, height = 0.2)
}

