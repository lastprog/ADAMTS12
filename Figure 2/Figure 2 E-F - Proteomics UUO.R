library(ggplot2)
library(ggrepel)
library(writexl)
library(tidyverse)
library(tibble)
library(readxl)

sample <- "UUO_WTvsKO_"
setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/Revision/17) Proteomics of Kidneys/")

# Read dataset
df<-read.delim(file="RawData/231108 protein two sample ttest.txt")
df$Student.s.T.test.Difference.ko_wt<-as.numeric(gsub(",", ".", df$Student.s.T.test.Difference.ko_wt))
df$Student.s.T.test.q.value.ko_wt<-as.numeric(gsub(",", ".", df$Student.s.T.test.q.value.ko_wt))
df$X.Log.Student.s.T.test.p.value.ko_wt<-as.numeric(gsub(",", ".", df$X.Log.Student.s.T.test.p.value.ko_wt))

# Vulcano_nice function with settable y-axis limit, because -log(P-Val) = Infinity
volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA, ylimab = 20, xlimab = 5)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- xlimab
    ylimAbs <- ylimab
  }
  else {
    xlimAbs <- xlimab
    ylimAbs <- ylimab
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(P.Val), color = couleur)) +
      geom_point(alpha = 1, size = 0.1) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "indianred",
                                     "darkblue")) + theme_minimal() + theme(legend.position = "none")+
      theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
            axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(P.Val),
                          color = couleur)) + geom_point(alpha = 1, size = 0.1) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "indianred",
                                       "darkblue")) + theme_minimal() + theme(legend.position = "none")+
        theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
              axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(P.Val),
                          color = couleur)) + geom_point(alpha = 1, size = 0.1) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "indianred",
                                       "darkblue")) + theme_minimal() + theme(legend.position = "none")+
        theme(axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
              axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"))
    }
  }
  return(a)
}

#Plot DEP

# Proteomic
# Format Table
df2<-df[-c(1,2), ]
df2 <- df2[c("Gene.names", "Protein.IDs","Student.s.T.test.Difference.ko_wt","X.Log.Student.s.T.test.p.value.ko_wt", "Student.s.T.test.q.value.ko_wt")]
write_xlsx(df2, paste0(sample, "01_SuppTable.xlsx"))

# For plotting
volcano_nice(df = df2, hAss = 0.05, FCIndex = 2,pValIndex = 4, vAss = 0.3,
             IDIndex = 1, label = T,nlabels = 20, ylimab = 10, xlimab = 10)+
  ggtitle(sample)+
  labs(x="logFC", y="-log(Q-Val.)")+
  theme(plot.title = element_text(size = 9, family ="Arial", hjust = 0.5),
        axis.text.x = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.text.y = element_text(size =6, face= "plain", family = "Arial", colour = "black"),
        axis.title = element_text(size =8, face= "plain", family = "Arial"),
        panel.grid.major = element_line(colour = "grey", size = 0.05),
        panel.grid.minor = element_line(colour = "grey", size = 0.05),
        panel.background = element_rect(fill="white", colour = "white"))
ggsave(filename = paste0(sample, "02_VolcanoPlot_DEP.jpg"), width = 4, height = 4)
ggsave(filename = paste0(sample, "02_VolcanoPlot_DEP.svg"), width = 4, height = 4)


# Process murine matrisome genes from #http://matrisomeproject.mit.edu/
matrisome_mm_masterlist <- as.data.frame(read_excel("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/matrisome_mm_masterlist.xls"))
matrisome_mm_masterlist<-matrisome_mm_masterlist[c("Division","Category","Gene Symbol")]
matrisome_mm_masterlist$Division=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Division=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Division)
matrisome_mm_masterlist$Category=gsub(pattern = " ",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist$Category=gsub(pattern = "-",replacement = "_",x = matrisome_mm_masterlist$Category)
matrisome_mm_masterlist.2<-matrisome_mm_masterlist
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Matrisome_associated"] = "Matrisome_associated"
matrisome_mm_masterlist.2$Category[matrisome_mm_masterlist$Division=="Core_matrisome"] = "Core_matrisome"
matrisome_mm_masterlist<-rbind(matrisome_mm_masterlist,matrisome_mm_masterlist.2)
matrisome_mm_masterlist<-matrisome_mm_masterlist[matrisome_mm_masterlist$Division!="Retired",]
rm(matrisome_mm_masterlist.2)
matrisome_mm_genesetlist = list()
for (geneset in unique(matrisome_mm_masterlist$Category)) {
  matrisome_mm_genesetlist[[geneset]] = matrisome_mm_masterlist$`Gene Symbol`[matrisome_mm_masterlist$Category==geneset]
}

df.filtered<- df[df$Student.s.T.test.q.value.ko_wt<0.05,]
df.filtered <- df.filtered %>%  dplyr::select(Gene.names, Student.s.T.test.Difference.ko_wt) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene.names)
ranks <- deframe(df.filtered)

library(fgsea)
fgseaRes <- fgsea(pathways=matrisome_mm_genesetlist, stats=ranks)
fgseaResTidy <- fgseaRes %>%  as_tibble()
writexl::write_xlsx(fgseaResTidy, path = paste0(sample, "03_Matrisome_GSEA_Table.xlsx")) 

# Show in a nice table:
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES) %>%  arrange(padj) %>%  DT::datatable()
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Naba Matrisome Category", y="Normalized Enrichment Score",
       title=paste0("Matrisome NES based on DEP")) + 
  theme_minimal()
ggsave(paste0(sample, "04_Matrisome_GSEA_Table.jpg"), bg = "white", width = 4, height=4)
ggsave(paste0(sample, "04_Matrisome_GSEA_Table.svg"), bg = "white", width = 4, height=4)


df.filtered<-df.filtered %>% arrange(desc(Student.s.T.test.Difference.ko_wt))
genes_up_df <-df.filtered[df.filtered$Student.s.T.test.Difference.ko_wt>0,]
genes_down_df <-df.filtered[df.filtered$Student.s.T.test.Difference.ko_wt<0,]
genes_up <- genes_up_df$Gene.names
genes_down <- genes_down_df$Gene.names

library(gprofiler2)
go_down <- gost(query = genes_down, 
                organism = "mmusculus", ordered_query = T, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
go_down.results = go_down$result
write_xlsx(go_down.results, path = paste0(sample, "06_GO_Down.xlsx"))

# Plotting Top 10 down GO:BP + REACTOME
# Plotting Top 10 down GO:BP
go.results.filter <- go_down.results[go_down.results$term_size<2000,] #filter out large and often very general GOs
go.results.filter.bp <- go.results.filter[go.results.filter$source=="GO:BP",]# filter out all Terms from Term GF
go.results.top10 <- go.results.filter.bp %>% top_n(n = -10,wt=p_value)
order_desc_up <- unique(go.results.top10$term_name)
go.results.top10 <- go.results.filter.bp[go.results.filter.bp$term_name %in% order_desc_up,]
ggplot(go.results.top10, aes(x=factor(term_name, level = order_desc_up), y=factor(source),size=precision,fill=-log10(p_value),color=-log10(p_value))) + 
  geom_point() + 
  scale_size(range = c(1,6)) +
  coord_flip() +
  ylab("") + xlab("") + ggtitle("Top GO:BP based on downregulated DE Proteins")+
  theme(axis.text.y = ,axis.text.x = element_text(angle=45,hjust=1), legend.title = element_text(size = 8),plot.title.position = "plot", 
        plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(sample, "07_Top10_DownGO_BP.jpeg"), width=5, height = 3.7)
ggsave(filename = paste0(sample, "07_Top10_DownGO_BP.svg"),  width=5, height = 3.7)
