# PROGGENY AND GO Analysis of BULK RNA Script
library(ggplot2)
library(reshape2)
library(writexl)
library(readxl)
library(plyr)
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library(readr)
library(ggrepel)
library(gridExtra)
library(graphics)
library(fgsea)
library(stringr)
library(gprofiler2)

windowsFonts("Arial" = windowsFont("Arial"))
setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 Überexpression/2021 ADAMTS12 Überexpression/Experiments/Bulk-Seq/Analysis/ExceptInactRep4")
sample="ADAMTS12_OvExp_"
ActvsKO<-read.delim(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 Überexpression/2021 ADAMTS12 Überexpression/Experiments/Bulk-Seq/Analysis/ExceptInactRep4/ADAMTS12_OvExp_00_DESEQ2_WTvsKO.txt")
ActvsInact<-read.delim(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 Überexpression/2021 ADAMTS12 Überexpression/Experiments/Bulk-Seq/Analysis/ExceptInactRep4/ADAMTS12_OvExp_00_DESEQ2_WTvsInact.txt")
InactvsKO<-read.delim(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 Überexpression/2021 ADAMTS12 Überexpression/Experiments/Bulk-Seq/Analysis/ExceptInactRep4/ADAMTS12_OvExp_00_DESEQ2_InactvsKO.txt")

Comparisons <- list(WTvsKO, WTvsInact, InactvsKO)
Comparisons.names <- c("ActvsKO", "ActvsInact", "InactvsKO")

for (i in 1:3) {
df<-Comparisons[[i]]
Current_Name=Comparisons.names[i]

colnames(df)<-c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")

# Remove genes with NA for p.values as these are very weakly expressed and therefore assigned NA values by DESeq comparison
df<-na.omit(df)
df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$log2FC > 0.5 & df$pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FC < -0.5 & df$pval < 0.05] <- "DOWN"
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$ID[df$diffexpressed != "NO"]

# PROGENy Pathway Analysis based on DEG
paletteLength <-100
myColor <-colorRampPalette(c("darkblue", "whitesmoke", "indianred"))(paletteLength)

df_matrix <- df %>% 
  dplyr::select(ID,t) %>% 
  dplyr::filter(!is.na(t)) %>%
  column_to_rownames(var="ID")%>%  
  as.matrix()

PathwayActivity_counts = progeny(df_matrix, scale=FALSE, organism = "Human", top=500 , perm=10000,z_scores = TRUE)
rownames(PathwayActivity_counts) = Current_Name
Activity_counts=as.vector(PathwayActivity_counts)
progenyBreaks = c(seq(min(Activity_counts),0,length.out=ceiling(paletteLength/2)+1), seq(max(Activity_counts)/paletteLength, max(Activity_counts), length.out = floor(paletteLength/2)))
PathwayActivity_zscore <- t(PathwayActivity_counts)

if (i == 1) {
  PathwayActivity_df <- PathwayActivity_zscore
} else {
  PathwayActivity_df <- cbind(PathwayActivity_df,PathwayActivity_zscore)
}


# GO Term Analysis based on Top upregulated genes
df<-df %>% arrange(desc(t))
genes_up_df <-df[df$t>0]
#genes_up <- genes_up_df$ID
genes_up <- c(head(df, 500)[,"ID"])
go_up <- gost(query = genes_up, 
                organism = "hsapiens", ordered_query = T, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
p<-gostplot(go_up, capped = TRUE, interactive =F)
publish_gostplot(p, highlight_terms = go_up$result$term_id[1:10] ,width = 10, height = 12, filename = paste0(sample, "03_GO_Up_",Current_Name, ".jpeg"))
go_up.new.results=NULL
go_up.new.results = go_up$result
go_up.new.results$comparison<-Current_Name

#add results to overall dataframe
if (i == 1) {
  go.results <- go_up.new.results
  
} else {
  go.results <- rbind(go.results,go_up.new.results)
}
}
writexl::write_xlsx(go.results, path = paste0(sample, "04_TopGOTerms_Test.xlsx"))


#Fig4k: Plotting PROGENy pathway enrichment per cluster
library(pheatmap)
pheatmap(PathwayActivity_df, cluster_cols = F,cluster_rows = T, filename = paste0(sample, "06_Progeny_H.pdf"), width = 3, height=4, cutree_rows = 1, treeheight_col = 10)
dev.off()


# Fig 4l: Plotting Top 5 enrichecd GO:BP per comparison
go.results.filter <- go.results[go.results$term_size<2000,] #filter out large and often very general GOs
go.results.filter.bp <- go.results.filter[go.results.filter$source=="GO:BP",]# filter out all Terms from Term GF
go.results.top5 <- go.results.filter.bp %>% group_by(comparison) %>% top_n(n = -5,wt=p_value)
order_desc_up <- unique(go.results.top5$term_name)
go.results.top5 <- go.results.filter.bp[go.results.filter.bp$term_name %in% order_desc_up,]

ggplot(go.results.top5, aes(x=factor(term_name, level = order_desc_up), y=factor(comparison),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  geom_point() + 
  scale_size(range = c(1,4)) +
  coord_flip() +
  ylab("Cluster") + xlab("") + ggtitle("Top GO:BP based on upregulated DE genes")+
  theme(axis.text.y = ,axis.text.x = element_text(angle=45,hjust=1),
      legend.title = element_text(size = 8),plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(sample, "05_Top5_GO_BP.jpeg"), width=5.2, height = 3.7)
ggsave(filename = paste0(sample, "05_Top5_GO_BP.svg"),  width=5.2, height = 3.7)


# Fig 4m-n: Scoring of the different signatures in Fibroblast Reference Map from Peisker et al., Nat. com
library(Seurat)
sc <- readRDS(file = "C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/ADAMTS12 Überexpression/2021 ADAMTS12 ?berexpression/Experiments/Bulk-Seq/Analysis/Fibroblast_integrated_filtered_processed_annotated.rds")
DimPlot(sc, group.by = "Annotated_Subsets", label = TRUE)
ggsave(filename = paste0(sample,"06_PVM_Dimplot.jpeg"), width=10 , height = 10)

midpoint <- c(0.35 ,0.1, 0.25)
for (i in 1:3) {
  df<-Comparisons[[i]]
  Current_Name=Comparisons.names[i]
  m <- midpoint[i]
  colnames(df)<-c("ID", "base.Meanx", "log2FC", "lfc.SE", "t", "pval", "adj.pval", "index")
  df.filtered <- filter(df, adj.pval < 0.05 & t >0)
  features<- str_to_title(df.filtered [,1])
  features <- list(features)
  DefaultAssay(sc)="RNA"
  ctrl_genes = 50 #important

  sc = AddModuleScore(object = sc, features = features, name = Current_Name, ctrl = ctrl_genes)
  FeaturePlot(sc, features = paste0(Current_Name, '1'),pt.size = 1,order = T,label = FALSE) + scale_colour_gradient2(low = 'darkblue', mid = 'lightgrey', high = 'red', midpoint = m)
  ggsave(filename = paste0(sample,"06_PVM_", Current_Name, "_F.jpeg"), width=10 , height = 10)
  ggsave(filename = paste0(sample,"06_PVM_", Current_Name, "_F.svg"), width=10 , height = 10)
  FeaturePlot(sc, features = paste0(Current_Name, '1'),pt.size = 1,order = T,label = TRUE) + scale_colour_gradient2(low = 'darkblue', mid = 'lightgrey', high = 'red', midpoint = m)
  ggsave(filename = paste0(sample,"06_PVM_", Current_Name, "_FL.jpeg"), width=10 , height = 10)
  VlnPlot(sc, features = paste0(Current_Name, '1'),  group.by="Annotated_Subsets",pt.size = 0, sort = "decreasing") +NoLegend() + coord_flip() +
    geom_boxplot(width=0.1,position = position_dodge(1),outlier.shape = NA,coef=0)+theme(
      axis.title = element_blank(),
      text = element_text(size=8, family = "Arial"),
      axis.text.x = element_text(size=7, family = "Arial", angle=45, hjust=1),
      axis.text.y = element_text(size=7, family = "Arial"),
      plot.title = element_text(size=9, family="Arial", face = "plain"),
      axis.line = element_line(size = 0.2), 
      axis.ticks = element_line(size = 0.2))
  ggsave(filename = paste(sample,"06_PVM_", Current_Name, "_V.jpeg"), width=2, height = 2)
  ggsave(filename = paste(sample,"06_PVM_", Current_Name, "_V.svg"), width=2, height = 2)
}
