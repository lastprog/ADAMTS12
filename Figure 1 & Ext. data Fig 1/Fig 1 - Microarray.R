# Written in 2023 by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(affy)
library(limma)
library(mouse4302.db)
library(mouse4302cdf)
library(tidyr)
library(tibble)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggrepel)
library(dplyr)
library(stringr)
library(progeny)
source("Source/Pseudobulking_support_functions.R")
source("Source/my_vulcano_nice.R")
source("Source/colorado_bold text ggplot2.R")
setwd("xxx") # Set your working directory here

# Read in Data and list to annotate genenames
data <- ReadAffy(celfile.path = "Raw Data/CEL files")
df_4302 <- read_excel("D:/Doktorarbeit offline/R/Auswertungen ADAMTS12/Microarray Gli1 kidney/Affymetrix GeneChip Mouse Genome 430 2.0 Array Annotation file.xlsx")

# Normalise Data and get expression table
data.rma.norm <- rma(data)
rma <- exprs(data.rma.norm)

# Rename Columns
rma <- as.data.frame(rma)
samples <- c("Sham 1", "Sham 2", "Sham 3", "UUO 1", "UUO 2", "UUO 3")
colnames(rma) <- samples
rma <- rownames_to_column(rma, var="Probeset ID")

# Summarize groups and samples in a dataframe
targets <- data.frame(sample = samples,
                      condition = gsub("[0-9]","",samples))

# Annotate Data and remove rows that dont have rownames associated with them
df_4302 <- df_4302[,c("Probeset ID","Sequence ID","Gene Symbol")]
rma <- left_join(rma, df_4302, by="Probeset ID")
rma <- column_to_rownames(rma, var="Probeset ID")
rma <- drop_na(rma, "Gene Symbol")

######################## Ext data Fig 1A - VlnPlot ###########################

# Melt dataframe for plotting
melted_df <- reshape::melt(rma)
melted_df$condition <- str_split_fixed(melted_df$variable, " ",2)[,1]
melted_df$NO <- str_split_fixed(melted_df$variable, " ",2)[,2]

# Plot VlnPlots
ggplot(melted_df, aes(x = variable, y = value, group = variable, fill = condition))+
  geom_violin()+
  labs(fill = "surgery")+
  ylab("Normalized gene expr. \n(log2)")+
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15), expand = c(0,0))+
  theme_classic()+
  theme(plot.title = element_blank(),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size =8, 
                                   family = "Arial",
                                   colour = "black",
                                   angle = 45, 
                                   hjust = 1),
        axis.text.y = element_text(size =6, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none",
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =6, face= "plain", family = "Arial"),
        legend.key.width= unit(0.2, 'cm'))
ggsave("Output/Ext data Fig 1A - VlnPlot.svg", height = 1.8, width = 2,units = "in")

######################## Ext data Fig 1B - PCA ################################

# Calculate PCA
pca <- prcomp(t(rma[complete.cases(rma),samples]), center = T, scale. = F)

# Calculate Variance explained by each component
explained <- (pca$sdev)^2 / sum(pca$sdev^2)v

# Create dataframe for plotting 
data.to.plot <- as.data.frame(cbind(pca$x[,"PC1"], pca$x[,"PC2"]))
data.to.plot$sample <- targets$sample
data.to.plot$condition <- targets$condition
colnames(data.to.plot) <- c("pc1", "pc2", "sample", "condition")

# Point Plot PCA 
ggplot(data.to.plot, aes(x=pc1, y=pc2, color=condition)) +
  geom_point(size= 4, alpha = 0.3) +
  geom_text_repel(aes(label=samples), size = 2)+
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_classic() +
  xlab(paste("PC1 (",round(100*explained[1],digits=2),'%)',sep=''))+
  ylab(paste("PC2 (",round(100*explained[2],digits=2),'%)',sep='')) +
  ylim(c(-max(abs(pca$x[,2])),max(abs(pca$x[,2]))+10))+
  scale_x_continuous(position = "bottom")+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none")
ggsave("Output/Ext data Fig 1B - PCA.svg", height = 1.8, width = 2,units = "in")

####################### Differential analysis ################################

comparisons <- list("SHAMvsUUO"=c(2,-1)) # Subtract the first condition (-1) from the second one (2)
limmaRes <- runLimma(measurements = rma[,1:6], targets = targets, comparisons = comparisons) 
ttop_UUOvsSHAM <- ttopFormatter(topTable(limmaRes[[1]], coef = 1, number = length(rma[,1]), adjust.method = "fdr"))
names(ttop_UUOvsSHAM)[1]<-"Probeset ID"
ttop_UUOvsSHAM <- left_join(ttop_UUOvsSHAM, df_4302, by="Probeset ID")
ttop_UUOvsSHAM <- arrange(.data = ttop_UUOvsSHAM, -t)
ttop_UUOvsSHAM <- ttop_UUOvsSHAM[!duplicated(ttop_UUOvsSHAM$`Gene Symbol`),]
write.csv(ttop_UUOvsSHAM, "Output/Supp. table 1 - DEGs Microarray.csv")

############################# Fig 1B - Vulcano ################################

volcano_nice(df = ttop_UUOvsSHAM, hAss = 0.05, FCIndex = 2,pValIndex = 5, IDIndex = 9, 
             vAss = 0.5, label = F, straight = FALSE, nlabels = 20) +
  theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.1),
        panel.grid.minor = element_line(colour = "grey", size = 0.1))
 ggsave(filename ="Output/Fig 1B - Vulcano.tiff", width = 2.2, height =1.7, units = "in", dpi = 1000)

############################# Fig 1C - Barplot ################################

up_down <- rbind(head(ttop_UUOvsSHAM, 7), tail(ttop_UUOvsSHAM,3))
ggplot(up_down,
       aes(x = reorder(factor(up_down$`Gene Symbol`), up_down$t), y = up_down$t))+
  geom_bar(aes(fill = up_down$t), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred",
                       mid = "whitesmoke", midpoint = mean(up_down$logFC))+
  labs(title = "Top DEG") +
  xlab("")+
  ylab("t-value")+
  theme_minimal()+
  scale_y_continuous(breaks = c(-20, 0, 20), limits = c(-20,20))+
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        legend.position = "none",
        axis.text.y = element_text(size =8, face= colorado(reorder(factor(up_down$`Gene Symbol`), up_down$t), "Adamts12"),
                                   family = "Arial", colour = "black"),
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =6, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank(), 
        legend.key.width= unit(0.2, 'cm'))+
  guides(fill = guide_legend(title = "t"))+
  coord_flip()
ggsave("Output/Fig 1C - barplot DEGs.svg", width = 1.3, height = 1.8, units = "in")

########################## Ext data Fig 1C - Progeny ########################### 

# Prepare data matrix for Progeny
ttop_UUOvsSHAM_matrix <- ttop_UUOvsSHAM %>% 
  dplyr::select(`Gene Symbol`, t) %>% 
  dplyr::filter(!is.na(t)) %>% 
  column_to_rownames(var = "Gene Symbol")%>%
  as.matrix()

# Run Progeny
PathwayActivity_zscore <- progeny(ttop_UUOvsSHAM_matrix, scale=TRUE, organism="Mouse", top = 100, perm = 10000, z_scores = TRUE) %>%t()

# Arrange results
colnames(PathwayActivity_zscore) <- "NES"
PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var="Pathway") %>%
  dplyr::arrange(NES) %>% mutate(Pathway=factor(Pathway))  

# Plot Progeny
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") + labs(title = "Progeny") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  theme_minimal() +
  xlab("Pathways")+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5), 
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        axis.text.y = element_text(size =8, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank()) +
  coord_flip()
ggsave("Output/Ext data Fig 1C - Progeny.svg", width = 1.2, height = 2.1, units = "in")

########################## Ext data Fig 1D - Dorothea ##########################

# Load Dorothea Regulons
data(dorothea_mm, package = "dorothea")
regulons <- dorothea_mm %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

# Run Dorothea
tf_activities_stat <- dorothea::run_viper(ttop_UUOvsSHAM_matrix, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::top_n(10, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

# Plot Dorothea
ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  labs(title = "Dorothea") +
  scale_y_continuous(breaks = c(-6,-3, 0, 3))+
  theme_minimal() +
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  theme(axis.title = element_text(size =8, face= "plain", family = "Arial",color = "black"),
        plot.title = element_text(face="plain", size = 10, hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", color = "black"),
        axis.text.y = element_text(size =8, face= "plain", family = "Arial", color = "black"),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank(),  
        legend.key.height = unit(0.1, "inches"),
        legend.key.width = unit(0.05, "inches"),
        legend.position = "none") +
  coord_flip()
ggsave("Output/Ext data Fig 1D - Dorothea.svg", width = 1.2, height = 2.1, units = "in")




