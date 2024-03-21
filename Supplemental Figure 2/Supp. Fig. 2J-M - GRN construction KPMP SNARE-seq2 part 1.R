# ADAMTS12 in fibrosis
# GRN construction for identification of upstream regulators using scMEGA - part 1
# Felix Schreibing [2024]

library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyverse)
library(MOJITOO)


# define input and output directory
indir = '~/GSE183279_Superseries/GSE183273_KPMP_SNARE_seq2/'

outdir = '~/GRNs_KPMP/KPMP_SNARE_seq2_objects/'


# PART 1: DATA PROCESSING
# we use the SNARE-seq 2 data set from:
# Lake, B.B., Menon, R., Winfree, S. et al. An atlas of healthy and injured cell states and niches in the human kidney. 
# Nature 619, 585–594 (2023). https://doi.org/10.1038/s41586-023-05769-3
# load data
ATAC = readRDS(paste0(indir, 'GSE183273_Kidney_Healthy-Injury_Cell_Atlas_SNARE2-AC_Peak-Counts_03282022.RDS'))
RNA = readRDS(paste0(indir, 'GSE183273_Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA_Counts_03282022.RDS'))
meta = read.delim(paste0(indir, 'GSE183273_Kidney_Healthy-Injury_Cell_Atlas_SNARE2-RNA-AC_Metadata_03282022.txt'))


# preprocess the snRNA-seq data
obj = CreateSeuratObject(counts = RNA,
                         assay = 'RNA',
                         meta.data = meta)

# subset snRNA-seq data for relevant cell types (fibroblasts and myofibroblasts)
obj = subset(obj, subset = subclass.l3 %in% c('FIB', 'M-FIB', 'MYOF', 'aFIB', 'dFIB'))
Idents(obj) = 'subclass.l3'

# add the snATAC-seq data
# subset peak-barcode matrix for relevant cells
ATAC = ATAC[, colnames(ATAC) %in% colnames(obj)]

# process the snATAC-seq data
obj[['ATAC']] = CreateChromatinAssay(counts = ATAC,
                                     sep = c(":", "-"),
                                     min.cells = 1,
                                     genome = 'hg38',
                                     fragments = paste0(indir, 'GSE183273_BUKMAP.fragments.sort.tsv.gz'))


# normalization and PCA of snRNA-seq data
DefaultAssay(obj) = 'RNA'

obj = obj %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = F) %>%
  ScaleData() %>%
  RunPCA(npcs = 50, reduction.name = 'RNA_PCA')


# gene annotation of snATAC-seq data
DefaultAssay(obj) = "ATAC"

# extract gene annotations from EnsDb
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86,
                                  verbose = FALSE)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) = 'UCSC'

# add the gene information to the object
Annotation(obj) = annotations

# normalization and PCA of snATAC-seq data
obj = obj %>%
  RunTFIDF() %>%
  FindTopFeatures(min.cutoff = 'q0', verbose = FALSE) %>%
  RunSVD()

# save the object as rds file
saveRDS(obj, paste0(outdir, 'combined_fibroblast_object.rds'))



# PART 2: RUN MOJITOO (not used)
# run MOJITOO
obj = mojitoo(object = obj,
              reduction.list = list("RNA_PCA", "lsi"),
              dims.list = list(1:50, 2:50), ## exclude 1st dimension of LSI
              reduction.name = 'MOJITOO',
              assay = "RNA")

DefaultAssay(obj) = "RNA"

embedd = Embeddings(obj[["MOJITOO"]])

obj = RunUMAP(obj,
              reduction = "MOJITOO",
              reduction.name = "MOJITOO_UMAP",
              dims = 1:ncol(embedd),
              verbose = FALSE)

DimPlot(obj,
        group.by = "subclass.l2",
        label = TRUE,
        reduction = "MOJITOO_UMAP")

saveRDS(obj, paste0(outdir, 'combined_fibroblast_object_MOJITOO.rds'))