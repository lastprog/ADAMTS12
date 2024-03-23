# Written by Hyojin Kim, adjusted by Lars Koch
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(vsn)
library(DESeq2)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
library(mgsub)
setwd("xxx") # Set your working directory here

# Read data
dat <- read.table("Input/filtered.count.txt")

# Select only Vehicle
dat <- dat[,grep("Veh", colnames(dat))]

# Define meta dataframe
design <- colnames(dat) %>% strsplit(., "_REP") %>% lapply(., function(x) x[1]) %>% unlist() %>% factor()
sample <- colnames(dat) %>% factor()
meta <- data.frame(sample, design)
colnames(meta) <- c("sample", "condition")
outdir <- "Output/"

###################### Supp. Fig 7D - Quality Control ######################

# Logarithmize
dat_log <- dat
dat_log[dat_log == 0] <- NA
dat_log <- log2(dat)
dat_log$ID <- rownames(dat_log)

# melt dataframe
melted_df <- reshape::melt(dat_log)

# Define Groups
melted_df$group <- str_split_fixed(melted_df$variable, pattern = "_",2)[,1]
melted_df$variable <- gsub("_Veh_REP" ," ",melted_df$variable)
melted_df$variable <- gsub("WT", "Active", melted_df$variable)
melted_df$variable <- factor(melted_df$variable, levels = c("Inact 4","Inact 3","Inact 2","Inact 1",
                                                            "Active 4","Active 3","Active 2","Active 1",
                                                            "KO 4","KO 3","KO 2","KO 1"))

# Plot Violin Plots
ggplot(melted_df, aes(x = variable, y = value, fill = group)) + 
  geom_violin()+
  ylab("log2 count")+
  scale_color_manual(
    values = c(
      "KO" = "#F8766D",
      "Active" = "#00BA38",
      "Inact" = "#619CFF"))+
  theme_classic()+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black",
                                   angle = 0, 
                                   hjust = 1),
        axis.text.y = element_text(size =8, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none",
        legend.text = element_text(size =6, face= "plain", family = "Arial"),
        legend.title = element_text(size =6, face= "plain", family = "Arial"),
        legend.key.width= unit(0.2, 'cm'))+
  coord_flip()
ggsave("Output/Supp Fig 4I - Violin Plots QC.svg", height = 2, width = 1.7, units = "in")

######################### Preprocess by DESEq2 #################################

dds <- DESeqDataSetFromMatrix(countData = round(dat,0), colData = meta, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

################### Supp. Fig. Fig 7E - PCA with Inact 4#########################

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Rename pcaData
pcaData$condition <- mgsub(as.character(pcaData$condition), c("Inact_Veh","KO_Veh","WT_Veh"), c("Inact", "KO", "Active"))
pcaData$name <- mgsub(as.character(pcaData$name),
                      c("Inact_Veh_REP1","Inact_Veh_REP2","Inact_Veh_REP3","Inact_Veh_REP4","KO_Veh_REP1","KO_Veh_REP2","KO_Veh_REP3","KO_Veh_REP4","WT_Veh_REP1","WT_Veh_REP2","WT_Veh_REP3","WT_Veh_REP4"),
                      c("Inact 1", "Inact 2", "Inact 3", "Inact 4", "KO 1", "KO 2", "KO 3", "KO 4", "Active 1", "Active 2", "Active 3", "Active 4"))

# Point Plot PCA 
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=3, alpha = 0.3) +
  geom_text_repel(aes(label=pcaData$name), size = 2)+
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_classic() +
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  scale_color_manual(
    values = c(
      "KO" = "#F8766D",
      "Active" = "#00BA38",
      "Inact" = "#619CFF"))+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none")
ggsave("Output/PCA_bulk_KO_WT_overexp.svg", height = 1.8, width = 2,units = "in")

################# Supp. Fig. 7F - PCA without Inact 4 #######################

# Remove "Inact_Veh_REP4"
dat <- select(dat, -Inact_Veh_REP4)

design <- colnames(dat) %>% strsplit(., "_REP") %>% lapply(., function(x) x[1]) %>% unlist() %>% factor()
sample <- colnames(dat) %>% factor()
meta <- data.frame(sample, design)
colnames(meta) <- c("sample", "condition")

######################### Preprocess by DESEq2 #################################
dds <- DESeqDataSetFromMatrix(countData = round(dat,0), colData = meta, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
################################################################################

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Rename pcaData
pcaData$condition <- mgsub(as.character(pcaData$condition), c("Inact_Veh","KO_Veh","WT_Veh"), c("Inact", "KO", "Active"))
pcaData$name <- mgsub(as.character(pcaData$name),
                      c("Inact_Veh_REP1","Inact_Veh_REP2","Inact_Veh_REP3","KO_Veh_REP1","KO_Veh_REP2","KO_Veh_REP3","KO_Veh_REP4","WT_Veh_REP1","WT_Veh_REP2","WT_Veh_REP3","WT_Veh_REP4"),
                      c("Inact 1", "Inact 2", "Inact 3", "KO 1", "KO 2", "KO 3", "KO 4", "Active 1", "Active 2", "Active 3", "Active 4"))

# Point Plot PCA 
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=3, alpha = 0.3) +
  geom_text_repel(aes(label=pcaData$name), size = 2)+
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_classic() +
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  scale_color_manual(
    values = c(
      "KO" = "#F8766D",
      "Active" = "#00BA38",
      "Inact" = "#619CFF"))+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none")
ggsave("Output/PCA_bulk_KO_WT_overexp.svg", height = 1.8, width = 2,units = "in")

