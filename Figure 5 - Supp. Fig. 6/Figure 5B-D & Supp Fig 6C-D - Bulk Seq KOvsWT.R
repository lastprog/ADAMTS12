# Adjusted by Lars Koch, on the basis of a script provided Hyojin Kim
# Contact: https://www.linkedin.com/in/lars-koch-93844023a/
# Publication: "ADAMTS12 promotes fibrosis by restructuring ECM to enable activation 
# and migration of injury-responsive fibroblasts"

library(dplyr)
library(tibble)
library(DESeq2)
library(ggplot2)
library(vsn)
library(stringr)
library(EnhancedVolcano)
library(progeny)
library(gprofiler2)
library(mgsub)
setwd("xxx") # Set your working directory here
source("Source/my_vulcano_nice.R")
source("Source/colorado_bold text ggplot2.R")

# Read Count Matrix
dat <- read.table("Input/salmon.merged.gene_counts.tsv", header = T)

# Arrange Count Matrix
# Exclude Genes which occure >1 time
dat <- dat[!(duplicated(dat$gene_name)),]
dat <- dat[, -grep("TGFb", colnames(dat))]
dat <- select(dat, -gene_id)
rownames(dat) <- NULL
dat <- column_to_rownames(dat, "gene_name")

# Create Meta df from Colnames
design <- colnames(dat) %>% strsplit(., "_REP") %>% lapply(., function(x) x[1]) %>% unlist() %>% factor()
sample <- colnames(dat) %>% factor()
meta <- data.frame(sample, design)
colnames(meta) <- c("sample", "condition")
outdir = "DESeq2/"

# Supp. Fig. 4C - Quality control

# Logarithmize
dat_log <- dat
dat_log[dat_log == 0] <- NA
dat_log <- log2(dat)
dat_log$ID <- rownames(dat_log)

# Melt dataframe
melted_df <- reshape::melt(dat_log)

# Define Groups
melted_df$group <- str_split_fixed(melted_df$variable, pattern = "_",2)[,1]
melted_df$variable <- gsub("_Vehicle_REP" ," ",melted_df$variable)
melted_df$variable <- factor(melted_df$variable, levels = c("KO 4", "KO 3",
                                                            "KO 2", "KO 1",
                                                            "WT 4", "WT 3",
                                                            "WT 2", "WT 1"))

# Plot Violin Plots
ggplot(melted_df, aes(x = variable, y = value, group = variable, fill = group)) + 
  geom_violin()+
  ylab("log2 count")+
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
 ggsave("Output/Ext. data Fig 4C - Quality control.svg", height = 2, width = 1.4, units = "in")

######################### Preprocess by DESEq2 #################################
 
dds <- DESeqDataSetFromMatrix(countData = round(dat,0), colData = meta, design = ~condition)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

# Supp. Fig. 6D - PCA

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c( "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$name <- gsub("_Vehicle_REP", " ", pcaData$name)

# Plot PCA 
ggplot(pcaData, aes(x=PC1, y=PC2, color= pcaData$group)) +
  geom_point(size=3, alpha = 0.3) +
  geom_text_repel(aes(label=pcaData$name), size = 2)+
  scale_alpha_discrete(range=c(0.3, 1.0)) +
  theme_classic() +
  xlab(paste0("PC1: ", percentVar[1], "% variance"))+
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  theme(plot.title = element_text(size = 9, face= "plain", hjust = 0.5),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        axis.text.x = element_text(size =6, 
                                   family = "Arial",
                                   colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", 
                                   family = "Arial", 
                                   colour = "black"),
        legend.position = "none")

ggsave(paste0(outdir,"/Ext data Fig 4D - PCA.svg"), height = 1.8, width = 2,units = "in")

############################# Calculate DEGs ##################################

get_DEgene <- function(dds,A,B,qval_cutoff, log2_cutoff, folder) {
  # ----------------------------------- #
  A_B <- c("condition", A, B)
  # ----------------------------------- #
  res_table_A_B <- results(dds, contrast=A_B, name=paste0(A, ".vs.", B), alpha = qval_cutoff)
  # ----------------------------------- #
  res_table_A_B["index"] = paste0(A, ".vs.", B)
  res_table_A_B_df = res_table_A_B %>% as.data.frame() %>% rownames_to_column(.) %>% as.data.frame() 
  # ----------------------------------- #
  res_table_A_B_df[ is.na(res_table_A_B_df) ] == ""	
  # ----------------------------------- #
  res = res_table_A_B_df
  x_max = res$log2FoldChange %>% max()
  x_min = res$log2FoldChange %>% min() %>% abs()
  x_lim_max <- max(x_max, x_min)
  # ----------------------------------- #
  p <- EnhancedVolcano(res, lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'padj',
                       xlim = c(x_lim_max*-1, x_lim_max),
                       pointSize = 3.0,
                       labSize = 1.0,
                       ylab = bquote(~-Log[10]~adjusted~italic(P)),
                       FCcutoff = log2_cutoff,
                       pCutoff = qval_cutoff
  )
  ggsave(paste0(outdir, folder, "/",  "DESeq2.volcano.log2FC", log2_cutoff, ".qval.", qval_cutoff, ".jpg"), plot = p,width=7, height=7)
  
  return (res_table_A_B_df)
}

## condition 
folders = "WT_Veh.vs.KO_Veh/"
B = "WT_Vehicle"
A = "KO_Vehicle"

## cutoff
qval_cutoff = 0.05
log2_cutoff = 0.05

# Calculate DEGs
A_ups <- get_DEgene(dds, A, B, qval_cutoff, log2_cutoff, folders)

# Save DEGs
deg <- na.omit(A_ups) # Remove genes with NA for p.values as these are very weakly expressed and therefore assigned NA values by DESeq comparison
colnames(deg) <- c("gene_name",  "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
deg <- deg[order(deg$t, decreasing = F),]

write.csv(deg, file=paste0(outdir, folders,"Supp. table 3 - DEGs Bulk Seq KOvsWT.csv"), row.names = F)

# Figure 5B - Vulcano Plot

# Volcano Plot with labeling the top 10 up- and downregulated genes ranked by T-Value
top10 <- deg %>% top_n(n = 10, wt = t) %>% pull(`gene_name`)
bottom10 <- deg  %>% top_n(n = -10, wt = t) %>% pull(`gene_name`)
labels <- c(top10,bottom10)

# For annotation
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 7, vAss = 0.3,
             IDIndex = 1, label = T,nlabels = 10, manual_labels = c(labels, "ADAMTS12"))+
  ggtitle("WT vs KO - DEG")+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0(outdir,"Vulcano_ann.jpg"), width = 20, height =20, units = "in", dpi = 1000)

# For Plotting
volcano_nice(df = deg, hAss = 0.05, FCIndex = 3,pValIndex = 7, 
             IDIndex = 1, vAss = 0.3,  nlabels = 10, label = F)+
  ggtitle("DEG - WT vs. ADAMTS12-KO")+
  labs(x="logFC", y="-log(P-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05))
ggsave(filename = paste0(outdir,"Fig  4B - Vulcano.svg"), width = 2.4, height =2, units = "in", dpi = 1000)

# Figure 4C - Progeny

# Prepare data matrix
df_matrix <- deg %>% 
  dplyr::select(gene_name,t) %>% 
  dplyr::filter(!is.na(t)) %>%
  column_to_rownames(var="gene_name")%>%
  as.matrix()

# Run Progeny
PathwayActivity_counts <- progeny(df_matrix, scale=FALSE, organism = "Human", top=500 , perm=10000,z_scores = TRUE) %>% t()
colnames(PathwayActivity_counts) = "NES"

# Arrange Progeny output
PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_counts) %>% 
  rownames_to_column(var="Pathway") %>%
  dplyr::arrange(NES) %>% mutate(Pathway=factor(Pathway))  

# Plot Progeny
ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, -NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") + 
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  labs(title = "Progeny")+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0)+
  scale_y_continuous(breaks = c(-16,-8,0))+
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
ggsave(filename = paste0(outdir, "/Fig 4C - Progeny.svg"),  width = 1.2, height = 2.1, units = "in")

# Fig 4D - GO:BP Terms

top100_down <- head(deg, 500)[,"gene_name"]

go_down <- gost(query = top100_down, 
                organism = "hsapiens", ordered_query = T, 
                multi_query = FALSE, significant = T, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
res_down <- go_down$result

# Filter GO:BP Terms with a term size <2000
res_down_filter <- res_down[res_down$term_size < 2000 & res_down$source =="GO:BP",]

# Prepare to 10 BP-GO Terms by -log10(p-value)
go_bp_down <- res_down_filter[res_down_filter$source == "GO:BP", ]
go_bp_down <- go_bp_down[order(go_bp_down$p_value),]
top10_go_bp <- head(go_bp_down, 10)
top10_go_bp$log10_pvalue <- log10(top10_go_bp$p_value)

# Abbreviate long names
top10_go_bp$term_name <- mgsub(top10_go_bp$term_name, 
                               c("tube development", "circulatory system development",
                                 "blood vessel development","vasculature development",
                                 "muscle organ development", "circulatory system process", "tube morphogenesis"),
                               c("tube dev", "circ system dev", "blood vessel dev",
                                 "vasc dev", "muscle dev", "circ system proc", "tube morph"))

# Plot bar Plot of Top 10 BP GO-terms downreg. in KO
ggplot(top10_go_bp,
       aes(x = reorder(factor(term_name), top10_go_bp$log10_pvalue, decreasing = T), y = top10_go_bp$log10_pvalue))+ 
  geom_bar(aes(fill = top10_go_bp$log10_pvalue), stat = "identity")+
  geom_hline(yintercept = 0,
             color = "black", size = 0.3)+
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = median(top10_go_bp$precision)) + 
  labs(title = "Top 10 BP-GO")+
  xlab("")+
  ylab("log10(p-value)")+
  theme_minimal()+
  theme(axis.title = element_text(size =8, face= "plain", family = "Arial"),
        plot.title = element_text(size = 9, face="plain",  hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial",
                                   color = "black"),
        axis.text.y = element_text(size =8, face= colorado(reorder(factor(top10_go_bp$term_name), 
                                                                   top10_go_bp$log10_pvalue, decreasing = T), 
                                                           c("cell adhesion",
                                                             "locomotion",
                                                             "cell migration")),
                                   family = "Arial",
                                   color = "black"),
        legend.position = "none",
        legend.text = element_text(size =8, face= "plain", family = "Arial"),
        legend.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_blank(),  
        legend.key.width= unit(0.4, 'cm'))+
  coord_flip()
ggsave(filename = paste0(outdir,"/Fig 4D - GO_BP Terms.svg"), width = 1.6, height = 2, units = "in")




