# ADAMTS12 in fibrosis
# GRN construction for identification of upstream regulators using scMEGA - part 2
# Felix Schreibing [2024]

library(Seurat)
library(Signac)
library(scMEGA)
library(ArchR)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyverse)
library(writexl)
library(ggplot2)
library(igraph)
library(GGally)
library(RColorBrewer)
set.seed(777)


# PRE-PROCESSING:
# KPMP SNARE-seq2 snRNA- and snATAC-seq data sets have been integrated (see part 1)
# and the integrated object has been subsetted for fibroblast populations

# define input and output directory
indir = '~/GRNs_KPMP/KPMP_SNARE_seq2_objects/'

outdir = '~/GRNs_KPMP/'



# load data
obj = readRDS(paste0(indir, 'combined_fibroblast_object.rds'))
DefaultAssay(obj) = 'RNA'



# PART 1: CREATE FIBROSIS BASED EMBEDDING
# we first create a "fibrosis-based" embedding of all cells
# which orders all cells according to their myofibroblast- and ECM-signature
# this allows a cluster-independent definition of myofibroblasts as the cells:
#   1.) with strongest ECM production
#   2.) with strongest myofibroblast signature
# load myofibroblast signature genes (TRAVAGLINI_LUNG_MYOFIBROBLAST_CELL from MSigDB)
myofib.sig = read.csv(paste0('~/TRAVAGLINI_LUNG_MYOFIBROBLAST_CELL/TRAVAGLINI_LUNG_MYOFIBROBLAST_CELL.txt'))
myofib.sig$X = NULL

# calculate myofibroblast gene signature
obj = AddModuleScore(obj,
                     features = list(myofib.sig$gene),
                     ctrl = length(myofib.sig$gene),
                     name = 'myofibroblast_signature',
                     search = TRUE)


# load core matrisome genes (NABA_CORE_MATRISOME)
matrisome = read.csv(paste0('~/NABA_CORE_MATRISOME/NABA_CORE_MATRISOME.v2023.1.Hs.txt'))
matrisome$X = NULL

# calculate myofibroblast gene signature
obj = AddModuleScore(obj,
                     features = list(matrisome$gene),
                     ctrl = length(matrisome$gene),
                     name = 'core_matrisome_score',
                     search = TRUE)

# get signature scores per cell
reduction = FetchData(obj, vars = c('myofibroblast_signature1', 'core_matrisome_score1'))
colnames(reduction) = c('FIB_1', 'FIB_2')

# add signature scores as embedding coordinates
reduction = CreateDimReducObject(embeddings = as.matrix(reduction),
                                 assay = DefaultAssay(obj))
obj@reductions$fib = reduction

# plot "fibrosis-based" cell embedding
plot.data = FetchData(obj, vars = c('FIB_1', 'FIB_2', 'myofibroblast_signature1', 'core_matrisome_score1'))
colnames(plot.data) = c('FIB_1', 'FIB_2', 'myofibroblast signature', 'NABA core matrisome score')
plot.data = gather(plot.data, key = 'signature', value = 'score', -c(FIB_1, FIB_2))

pdf(paste0(outdir, 'fibrosis_based_embedding_score_projection.pdf'), width = 7, height = 3.5)
ggplot(plot.data, aes(x = FIB_1, y = FIB_2, color = score)) +
  geom_point() +
  facet_wrap(~ signature) +
  theme_bw() + 
  scale_color_viridis(option = 'B')
dev.off()


# create new cell identities based on the "fibrosis embedding"
new.idents = FetchData(obj, vars = c('FIB_1', 'FIB_2'))
new.idents$type = ifelse(new.idents$FIB_1 > 0 & new.idents$FIB_2 > (max(new.idents$FIB_2)/2),
                         'myofibroblasts', 'fibroblasts')
new.idents = dplyr::select(new.idents, type)
obj = AddMetaData(obj, metadata = new.idents)
Idents(obj) = 'type'


# plot cell identity projection and relevant marker genes
plot.data = FetchData(obj, vars = c('FIB_1', 'FIB_2', 'type', 'ACTA2', 'POSTN', 'DCN', 'COL1A1', 'ADAMTS12'))

pdf(paste0(outdir, 'score_based_cell_idents.pdf'), width = 6.5, height = 5)
ggplot(plot.data, aes(x = FIB_1, y = FIB_2, fill = type)) +
  geom_point(shape = 21, color = 'black') +
  theme_bw() +
  scale_fill_brewer(palette = 'Reds') +
  labs(title = 'score-based cell identities')
dev.off()

plot.data = gather(plot.data, key = 'gene', value = 'expression', -c(FIB_1, FIB_2, type))

pdf(paste0(outdir, 'marker_expression_score_based_cell_idents.pdf'), width = 9, height = 3)
ggplot(plot.data, aes(x = type, y = expression, fill = type)) +
  geom_violin(color = 'black') +
  geom_jitter(size = 0.05) +
  facet_wrap(~ gene, nrow = 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_brewer(palette = 'Reds') +
  labs(title = 'score-based cell identities')
dev.off()



# PART 2: INFER TRAJECTORIES
# infer trajectories based on "fibrosis embedding"
obj = AddTrajectory(object = obj,
                    trajectory = c('fibroblasts', 'myofibroblasts'),
                    group.by = 'type',
                    reduction = 'fib',
                    dims = 1:2,
                    use.all = FALSE)

# we only keep the cells that are part of the trajectory
obj = obj[, !is.na(obj$Trajectory)]

# plot trajectory
plot.data = FetchData(obj, vars = c('FIB_1', 'FIB_2', 'Trajectory'))

pdf(paste0(outdir, 'fibroblast_trajectory.pdf'), width = 5, height = 4)
ggplot(plot.data, aes(x = FIB_1, y = FIB_2, fill = Trajectory)) +
  geom_point(shape = 21, color = 'black') +
  theme_bw() +
  scale_fill_distiller(palette = 'RdBu')
dev.off()



# PART 3: INFER TF ACTIVITY USING chromVAR
# get list of motif position frequency matrices from the JASPAR database
pfm = getMatrixSet(x = JASPAR2020,
                   opts = list(collection = "CORE",
                               tax_group = 'vertebrates',
                               all_versions = FALSE))

# add motif information
obj = AddMotifs(object = obj,
                genome = BSgenome.Hsapiens.UCSC.hg38,
                pfm = pfm,
                assay = 'ATAC')

# run chromVAR
obj = RunChromVAR(object = obj,
                  genome = BSgenome.Hsapiens.UCSC.hg38,
                  assay = 'ATAC')

saveRDS(obj, paste0(outdir, 'scMEGA_object.RDS'))



# PART 4: SELECT TFs AND GENES
# select transcription factors
sel.tfs = SelectTFs(object = obj,
                    return.heatmap = TRUE,
                    cor.cutoff = 0.4)

df.cor = sel.tfs$tfs

write_xlsx(df.cor, paste0(outdir, 'selected_transcription_factors.xlsx'))

# select genes
sel.genes = SelectGenes(object = obj,
                        var.cutoff.gene = 0.9,
                        labelTop1 = 0,
                        labelTop2 = 0,
                        fdr.cutoff = 0.05)

df.p2g = sel.genes$p2g

write_xlsx(df.p2g, paste0(outdir, 'peak_to_gene_correlations.xlsx'))



# PART 5: INFER GRNs
tf.gene.cor = GetTFGeneCorrelation(object = obj,
                                   tf.use = df.cor$tfs,
                                   gene.use = unique(df.p2g$gene),
                                   tf.assay = "chromvar",
                                   gene.assay = "RNA",
                                   trajectory.name = "Trajectory")

write_xlsx(tf.gene.cor, paste0(outdir, 'tf_gene_correlation.xlsx'))

# build GRNs
# if gene is regulated by a peak and this peak is bound by a TF
# this gene is considered to be a target of this TF
motif.matching = obj@assays$ATAC@motifs@data
colnames(motif.matching) = obj@assays$ATAC@motifs@motif.names
motif.matching = motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]

df.grn = GetGRN(motif.matching = motif.matching,
                df.cor = tf.gene.cor,
                df.p2g = df.p2g)

write_xlsx(df.grn, paste0(outdir, 'GRN_data_frame.xlsx'))



# PART 6: GRN VISUALIZATION
df.cor = df.cor[order(df.cor$time_point), ]
tfs.timepoint = df.cor$time_point
names(tfs.timepoint) = df.cor$tfs

# plot the whole network
df.grn2 = df.grn %>%
  subset(correlation > 0.5) %>%
  dplyr::select(c(tf, gene, correlation))

colnames(df.grn2)[colnames(df.grn2) == 'correlation'] = 'weights'

write_xlsx(df.grn2, paste0(outdir, 'GRN2_data_frame_for_plotting.xlsx'))

whole.network = GRNPlot(df.grn2,
                        tfs.timepoint = tfs.timepoint,
                        show.tf.labels = TRUE,
                        seed = 42,
                        plot.importance = FALSE,
                        min.importance = 2,
                        remove.isolated = FALSE)

pdf(paste0(outdir, 'GRN_whole_network.pdf'), width = 17, height = 15)
whole.network
dev.off()


# investigate activity of individual TFs
obj = AddTargetAssay(object = obj,
                     df.grn = df.grn2)

JUNB.plot = PseudotimePlot(object = obj, tf.use = "JUNB")
BACH1.plot = PseudotimePlot(object = obj, tf.use = "BACH1")

pdf(paste0(outdir, 'JUNB_activity_and_expression_plot.pdf'), width = 5.5, height = 4)
JUNB.plot +
  theme_bw() +
  labs(y = 'value')
dev.off()

pdf(paste0(outdir, 'BACH1_activity_and_expression_plot.pdf'), width = 5.5, height = 4)
BACH1.plot +
  theme_bw() +
  labs(y = 'value')
dev.off()


# visualize sub-networks for "JUNB" and "BACH1"
# subset the data frame for relevant transcription factors
df.grn3 = df.grn2 %>%
  dplyr::filter(weights > 0.8) %>%
  dplyr::filter(tf %in% c('JUNB', 'BACH1'))

# create a data frame which stores node meta data
vertex.data = data.frame(gene = unique(c(df.grn3$gene, df.grn3$tf)))

# define groups of genes, which will be highlighted in the network
vertex.data$type = ifelse(vertex.data$gene %in% matrisome$gene, 'core matrisome gene',
                          ifelse(vertex.data$gene %in% c('JUNB', 'BACH1'), 'central TF',
                                 ifelse(vertex.data$gene == 'ADAMTS12', 'ADAMTS12', 'other')))

# create an igraph-object
sub.network = graph_from_data_frame(df.grn3, vertices = vertex.data)

# visualize
sub.net = ggnet2(sub.network,
                 label = c('JUNB', 'BACH1', 'ADAMTS12'),
                 color = 'type',
                 palette = c('other' = 'gray',
                             'central TF' = '#67A9Cf',
                             'core matrisome gene' = '#F4A582',
                             'ADAMTS12' = '#B2182B'),
                 size = 'type',
                 size.palette = c('other' = 1,
                                  'central TF' = 10,
                                  'core matrisome gene' = 1,
                                  'ADAMTS12' = 5))

pdf(paste0(outdir, 'GRN_sub_network_JUNB_and_BACH1_v1.pdf'), width = 10, height = 7)
sub.net
dev.off()

pdf(paste0(outdir, 'GRN_sub_network_JUNB_and_BACH1_v2.pdf'), width = 15, height = 10)
sub.net
dev.off()


# plot expression of selected marker genes across pseudotime
data.plot = FetchData(obj, vars = c('Trajectory', 'ADAMTS12', 'COL1A1', 'POSTN', 'DCN', 'ACTA2'))

pdf(paste0(outdir, 'ADAMTS12_expression_across_pseudotime.pdf'), width = 5, height = 5)
ggplot(data.plot, aes(x = Trajectory, y = ADAMTS12)) +
  geom_smooth(color = '#B2182B') +
  theme_bw() +
  labs(title = 'ADAMTS12 expression',
       x = 'pseudotime',
       y = 'expression')
dev.off()

pdf(paste0(outdir, 'COL1A1_expression_across_pseudotime.pdf'), width = 5, height = 5)
ggplot(data.plot, aes(x = Trajectory, y = COL1A1)) +
  geom_smooth(color = '#B2182B') +
  theme_bw() +
  labs(title = 'COL1A1 expression',
       x = 'pseudotime',
       y = 'expression')
dev.off()

pdf(paste0(outdir, 'POSTN_expression_across_pseudotime.pdf'), width = 5, height = 5)
ggplot(data.plot, aes(x = Trajectory, y = POSTN)) +
  geom_smooth(color = '#B2182B') +
  theme_bw() +
  labs(title = 'POSTN expression',
       x = 'pseudotime',
       y = 'expression')
dev.off()

pdf(paste0(outdir, 'ACTA2_expression_across_pseudotime.pdf'), width = 5, height = 5)
ggplot(data.plot, aes(x = Trajectory, y = ACTA2)) +
  geom_smooth(color = '#B2182B') +
  theme_bw() +
  labs(title = 'ACTA2 expression',
       x = 'pseudotime',
       y = 'expression')
dev.off()

pdf(paste0(outdir, 'DCN_expression_across_pseudotime.pdf'), width = 5, height = 5)
ggplot(data.plot, aes(x = Trajectory, y = DCN)) +
  geom_smooth(color = '#B2182B') +
  theme_bw() +
  labs(title = 'DCN expression',
       x = 'pseudotime',
       y = 'expression')
dev.off()
