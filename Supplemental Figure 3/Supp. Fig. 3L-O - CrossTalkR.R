#Load libraries
library(Seurat) 
library(tidyverse)
library(dplyr)
library(harmony)
library(zellkonverter)
library(SeuratDisk)
require(biomaRt)
require(SingleCellExperiment)
require(liana)
require(purrr)
library(ggplot2)
require(reticulate)
library(magrittr)
library(CrossTalkeR)
library(SeuratObject)
require(tibble)
library(patchwork)
library(rmarkdown)
library(tinytex)
#Read file
SeuratObj = readRDS('/path/to/all_samples.rds')
head(SeuratObj@meta.data)
Idents(SeuratObj) = SeuratObj$max_like_cell
#remove special characters
SeuratObj$condition = SeuratObj$orig.ident
SeuratObj$condition <- gsub("1", "", SeuratObj$condition)
SeuratObj$condition <- gsub("2", "", SeuratObj$condition)
#Converting seurat object to sce
data = as.SingleCellExperiment(SeuratObj)
head(SeuratObj@meta.data)
#separate data on basis of condition
gp <-assays(data)[["counts"]][,colData(data)$condition=='WT.']
pp <-assays(data)[["counts"]][,colData(data)$condition=='KO.']
#create metadata
metadata_gp <-tibble::tibble(rownames(colData(data))[colData(data)$condition=='WT.'],
                             as.character(colData(data)$max_like_cell[colData(data)$condition=='WT.']))
metadata_pp<- tibble::tibble(rownames(colData(data))[colData(data)$condition=='KO.'],
                             as.character(colData(data)$max_like_cell[colData(data)$condition=='KO.']))
#mapping data to human
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#exclude 
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(gp) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
mapping <- genesV2[genesV2$HGNC.symbol!='',]
mappinglist <- mapping$HGNC.symbol
names(mappinglist) <- mapping$MGI.symbol
rownames(gp) <- mappinglist[rownames(gp)]
gp <- gp[!is.na(rownames(gp)),]
rownames(pp) <- mappinglist[rownames(pp)]
pp <- pp[!is.na(rownames(pp)),]
#create seurat object for running liana for conditions separately
good = CreateSeuratObject(counts = gp)
Idents(good) = metadata_gp$`as.character(...)`
poor = CreateSeuratObject(counts = pp)
Idents(poor) = metadata_pp$`as.character(...)`
#Run liana for both conditions
good = NormalizeData(good, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test <-  liana_wrap(good)
poor = NormalizeData(poor, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test2 <- liana_wrap(poor)
saveRDS(liana_test, file = 'wildtype.rds')
saveRDS(liana_test2, file = 'knockoout.rds')

#reloading liana results
liana_test = readRDS('wildtype.rds')
liana_test2 = readRDS('knockoout.rds')
#Selecting cellphonedb results for downstream
gud <-liana_test$cellphonedb%>%
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 
pur <-liana_test2$cellphonedb %>% 
  filter(pvalue<=0.05) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 
#create a list for input in generating CCI
l = list()
l[['WT']] = gud
l[['KO']] = pur
#generate output
data_dis <- generate_report(lrpaths = l,
                            genes=NULL,
                            out_path='~/output/path/',
                            threshold=0,
                            out_file = 'spatial_.html',
                            output_fmt = "html_document",
                            report = T,
                            sel_columns = c('source','target','gene_A','gene_B','type_gene_A','type_gene_B','lr.mean'))


#read output LR results with scores to get individual plots  
data = readRDS('LR_data_final.Rds')

pdf(file = 'cci.pdf')
options(repr.plot.width=10, repr.plot.height=10)
plot_cci(graph = data@graphs$WT, ## Selecting the Cell Cell Network
         colors = data@colors[V(data@graphs$WT)], ## Setting the color Pallet
         plt_name = 'WT', ## Plot title
         coords = data@coords[V(data@graphs$WT),], ## Nodes Coordinates
         emax = NULL,
         leg = TRUE, ## Setting Legend
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12,
         pg=data@rankings$WT$Pagerank[V(data@graphs$WT)])

plot_cci(graph = data@graphs$KO, ## Selecting the Cell Cell Network
         colors = data@colors[V(data@graphs$KO)], ## Setting the color Pallet
         plt_name = 'KO', ## Plot title
         coords = data@coords[V(data@graphs$KO),], ## Nodes Coordinates
         emax = NULL,
         leg = TRUE, ## Setting Legend
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12,
         pg=data@rankings$KO$Pagerank[V(data@graphs$KO)])

plot_cci(graph = data@graphs$KO_x_WT, ## Selecting the Cell Cell Network
         colors = data@colors[V(data@graphs$KO_x_WT)], ## Setting the color Pallet
         plt_name = 'KO_x_WT', ## Plot title
         coords = data@coords[V(data@graphs$KO_x_WT),], ## Nodes Coordinates
         emax = NULL,
         leg = TRUE, ## Setting Legend
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12,
         pg=data@rankings$KO_x_WT$Pagerank[V(data@graphs$KO_x_WT)])


dev.off()
# comparative analysis
pdf(file = 'network.pdf')
plot_cci(graph = data@graphs$KO_x_WT, ## Selecting the Cell Cell Network
         colors = data@colors[V(data@graphs$KO_x_WT)], ## Setting the color Pallet
         plt_name = 'KO_x_WT', ## Plot title
         coords = data@coords[V(data@graphs$KO_x_WT),], ## Nodes Coordinates
         emax = NULL,
         leg = TRUE, ## Setting Legend
         low = 0,
         high = 0,
         ignore_alpha = FALSE,
         log = FALSE,
         efactor = 8,
         vfactor = 12,
         pg=data@rankings$KO_x_WT$Pagerank[V(data@graphs$KO_x_WT)])
dev.off()

x <- max(abs(data@pca$KO_x_WT_ggi$x[,1]))
y <- max(abs(data@pca$KO_x_WT_ggi$x[,2]))
z_x <- data@pca$KO_x_WT_ggi$x[,1]
z_y <- data@pca$KO_x_WT_ggi$x[,2]
ver_zx <- ifelse(abs(z_x)>=(2*data@pca$KO_x_WT_ggi$sdev[1]),1,0)
ver_zy <- ifelse(abs(z_y)>=(2*data@pca$KO_x_WT_ggi$sdev[2]),1,0)
## Setting max.overlaps=25 you can increase or decrease the labels
## All the labels needs to be at 2 stdev in both PCA axis
pdf('pca_category.pdf', width = 20, height = 20)
fviz_pca_biplot(data@pca$KO_x_WT_ggi,
                axes = c(1,2),
                pointshape = 21, pointsize = 0.5,labelsize = 10,
                repel = FALSE,max.overlaps=200,label='var')+
  geom_label_repel(aes(label=ifelse((ver_zx | ver_zy),rownames(data@pca$KO_x_WT_ggi$x),NA)),size = 5,max.overlaps=25)+                     
  xlim(-x, x)+
  ylim(-y, y)+
  ggtitle('PCA example')+
  theme(text = element_text(size = 7.5),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.5))
dev.off()


#sankey plots

data = readRDS('LR_data_final.Rds')
df = data@tables$KO_x_WT

pdf(file = 'sankeyplot_fibroblast.pdf', width = 12, height = 6)
plot_sankey(lrobj_tbl = df, ## Selecting the Table
            ligand_cluster = "Fibroblast",
            plt_name = "Fibroblast as ligand(source)",
            threshold = 10) ### 10 highest expressed (absolute value) pairs
plot_sankey(lrobj_tbl = df, ## Selecting the Table
            receptor_cluster = "Fibroblast",
            plt_name = "Fibroblast as receptor(target)",
            threshold = 10) ### 10 highest expressed (absolute value) pairs

dev.off()

pdf(file = 'sankeyplot_lig.pdf', width = 12, height = 6)
plot_sankey(lrobj_tbl = df, ## Selecting the Table
            target = 'COL1A1',
            plt_name = "COL1A1 as ligand",
            threshold = 10) ### 10 highest expressed (absolute value) pairs

plot_sankey(lrobj_tbl = df, ## Selecting the Table
            target = 'COL1A2',
            plt_name = "COL1A2 as ligand",
            threshold = 10) ### 10 highest expressed (absolute value) pairs

plot_sankey(lrobj_tbl = df, ## Selecting the Table
            target = 'FN1',
            plt_name = "FN1 as ligand",
            threshold = 10) ### 10 highest expressed (absolute value) pairs
dev.off()



