library(tidyverse)
library(fgsea)
library(readxl)
library(gprofiler2)
library(writexl)

setwd("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/Revision/FSGEA ECM In Vitro/")
sample <- "ADAMTS12_MassSpec_"
df<-read.delim(file="C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/3. Adamts12/in vitro/CRISPR Cloning/2020 huADAMTS12 CRISPR-KO/Experimente/Massenspektrometrie/KH21-05_ECMIsolation2/220629_DEG_ProteinLevel.txt")
df$Student.s.T.test.Difference.NT_ADAMTS12ko<-as.numeric(gsub(",", ".", df$Student.s.T.test.Difference.NT_ADAMTS12ko))
df$Student.s.T.test.q.value.NT_ADAMTS12ko<-as.numeric(gsub(",", ".", df$Student.s.T.test.q.value.NT_ADAMTS12ko))
write_xlsx(df, path = paste0(sample, "df.xlsx"))

# Process human matrisome genes from #http://matrisomeproject.mit.edu/
matrisome_hs_masterlist <- as.data.frame(read_excel("C:/Users/Konra/Dropbox/Kramann Arbeitsgruppe/1. Monocyte & Pericyte Interaction/R-Scripts/methods/matrisome_hs_masterlist.xls"))
matrisome_hs_masterlist<-matrisome_hs_masterlist[c("Division","Category","Gene Symbol")]
matrisome_hs_masterlist$Division=gsub(pattern = " ",replacement = "_",x = matrisome_hs_masterlist$Division)
matrisome_hs_masterlist$Division=gsub(pattern = "-",replacement = "_",x = matrisome_hs_masterlist$Division)
matrisome_hs_masterlist$Category=gsub(pattern = " ",replacement = "_",x = matrisome_hs_masterlist$Category)
matrisome_hs_masterlist$Category=gsub(pattern = "-",replacement = "_",x = matrisome_hs_masterlist$Category)
matrisome_hs_masterlist.2<-matrisome_hs_masterlist
matrisome_hs_masterlist.2$Category[matrisome_hs_masterlist$Division=="Matrisome_associated"] = "Matrisome_associated"
matrisome_hs_masterlist.2$Category[matrisome_hs_masterlist$Division=="Core_matrisome"] = "Core_matrisome"
matrisome_hs_masterlist<-rbind(matrisome_hs_masterlist,matrisome_hs_masterlist.2)
matrisome_hs_masterlist<-matrisome_hs_masterlist[matrisome_hs_masterlist$Division!="Retired",]
rm(matrisome_hs_masterlist.2)
matrisome_hs_genesetlist = list()
for (geneset in unique(matrisome_hs_masterlist$Category)) {
  matrisome_hs_genesetlist[[geneset]] = matrisome_hs_masterlist$`Gene Symbol`[matrisome_hs_masterlist$Category==geneset]
}

df.filtered <- df %>%  dplyr::select(PG.Genes, Student.s.T.test.Difference.NT_ADAMTS12ko) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(PG.Genes)

ranks <- deframe(df.filtered)

fgseaRes <- fgsea(pathways=matrisome_hs_genesetlist, stats=ranks)
fgseaResTidy <- fgseaRes %>%  as_tibble() %>%  arrange(desc(NES))
writexl::write_xlsx(fgseaResTidy, path = paste0(sample, "01_Matrisome_GSEA_Table.xlsx")) 

# Show in a nice table:
fgseaResTidy %>% dplyr::select(-leadingEdge, -ES) %>%  arrange(padj) %>%  DT::datatable()

# Supp. Fig 9F
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Naba Matrisome Category", y="Normalized Enrichment Score",
       title=paste0("Matrisome NES based on DEP")) + 
  theme_minimal()
ggsave(paste0(sample, "02_Matrisome_GSEA.svg"), bg = "white", width = 4, height=4)

######################## GO:BP Terms ######################## 
df<-df %>% arrange(desc(Student.s.T.test.Difference.NT_ADAMTS12ko))

genes_up_df <-df[df$Student.s.T.test.Difference.NT_ADAMTS12ko>0,]
genes_down_df <-df[df$Student.s.T.test.Difference.NT_ADAMTS12ko<0,]
genes_up <- genes_up_df$PG.Genes
genes_down <- genes_down_df$PG.Genes

go_up <- gost(query = genes_up, 
              organism = "hsapiens", ordered_query = T, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              domain_scope = "annotated", custom_bg = NULL, 
              numeric_ns = "", sources = NULL, as_short_link = FALSE)
go_up.results = go_up$result
write_xlsx(go_up.results, path = paste0(sample, "03_GO_Up.xlsx"))

go_down <- gost(query = genes_down, 
              organism = "hsapiens", ordered_query = T, 
              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, evcodes = FALSE, 
              user_threshold = 0.05, correction_method = "g_SCS", 
              domain_scope = "annotated", custom_bg = NULL, 
              numeric_ns = "", sources = NULL, as_short_link = FALSE)
go_down.results = go_down$result
write_xlsx(go_down.results, path = paste0(sample, "04_GO_Down.xlsx"))

# Supp. Fig 9E
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
ggsave(filename = paste0(sample, "05_Top10_DownGO_BP.svg"),  width=4.6, height = 3.8)
