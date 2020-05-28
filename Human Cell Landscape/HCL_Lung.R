######################################################################
#Analysis of Lung scRNAseq data from the Human Cell Landscape
######################################################################

#Load libraries and processed data

reticulate::use_python('/usr/bin/python3')
library(reticulate)
library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(plotly)

setwd("~/Human Cell Landscape") #Change working directory to a dedicated folder

download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134355&format=file', 
              destfile = '~/Human Cell Landscape/GSE134355_RAW.tar', method = 'wget')
untar(tarfile = 'GSE122960.tar', exdir = "data")

counts1 <- read.table('data/GSM4008628_Adult-Lung1_dge.txt', sep = '\t', row.names = 1, header = T)
counts2 <- read.table('data/GSM4008629_Adult-Lung2_dge.txt', sep = '\t', row.names = 1, header = T)
counts31 <- read.table('data/GSM4008630_Adult-Lung3-1_dge.txt', sep = '\t', row.names = 1, header = T)
counts32 <- read.table('data/GSM4008631_Adult-Lung3-2_dge.txt', sep = '\t', row.names = 1, header = T)
counts33 <- read.table('data/GSM4008632_Adult-Lung3-3_dge.txt', sep = '\t', row.names = 1, header = T)
counts34 <- read.table('data/GSM4008633_Adult-Lung3-4_dge.txt', sep = '\t', row.names = 1, header = T)

dat1 <- CreateSeuratObject(counts1)
dat2 <- CreateSeuratObject(counts2)
dat31 <- CreateSeuratObject(counts31)
dat32 <- CreateSeuratObject(counts32)
dat33 <- CreateSeuratObject(counts33)
dat34 <- CreateSeuratObject(counts34)

dat <- merge(dat1, list(dat2, dat31, dat32, dat33, dat34))
#Add some metadata
s1 <- rep('Sample 1', times = ncol(dat1))
s2 <- rep('Sample 2', times = ncol(dat2))
s31 <- rep('Sample 3-1', times = ncol(dat31))
s32 <- rep('Sample 3-2', times = ncol(dat32))
s33 <- rep('Sample 3-3', times = ncol(dat33))
s34 <- rep('Sample 3-4', times = ncol(dat34))
sample_label <- c(s1, s2, s31, s32, s33, s34)
names(sample_label) <- colnames(dat)
dat <- AddMetaData(dat, metadata = sample_label, col.name = 'Sample')

######################################################################
#QC and filtering 
######################################################################

sets <- SplitObject(dat, split.by = 'Sample')
for (i in 1:length(sets)){
  ob <- sets[[i]]
  
  ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^MT-", assay = 'RNA')
  percent.mito <- PercentageFeatureSet(ob, pattern = "^MT-", assay = 'RNA')
  percent.mito <- percent.mito$nCount_RNA
  
  counts_per_cell <- Matrix::colSums(ob)
  counts_per_gene <- Matrix::rowSums(ob)
  genes_per_cell <- Matrix::colSums(ob@assays$RNA@counts >0)
  
  plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
  plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)', )
  plot(sort(percent.mito), xlab='cell', log='y', main='percent.mito (ordered)')
  
  sets[[i]] <- ob
}

sets[[1]] <- subset(sets[[1]], subset = nFeature_RNA > 250 & nFeature_RNA < 1000 &
                                nCount_RNA > 500 & nCount_RNA < 1800 &
                                percent.mt < 20)

sets[[2]] <- subset(sets[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 1000 &
                      nCount_RNA > 500 & nCount_RNA < 1800 &
                      percent.mt < 20)

sets[[3]] <- subset(sets[[3]], subset = nFeature_RNA > 100 & nFeature_RNA < 1000 &
                      nCount_RNA > 400 & nCount_RNA < 2000 &
                      percent.mt < 20)

sets[[4]] <- subset(sets[[4]], subset = nFeature_RNA > 100 & nFeature_RNA < 1000 &
                      nCount_RNA > 200 & nCount_RNA < 2000 &
                      percent.mt < 20)

sets[[5]] <- subset(sets[[5]], subset = nFeature_RNA > 100 & nFeature_RNA < 1000 &
                      nCount_RNA > 200 & nCount_RNA < 1800 &
                      percent.mt < 20)

sets[[6]] <- subset(sets[[6]], subset = nFeature_RNA > 100 & nFeature_RNA < 1000 &
                      nCount_RNA > 200 & nCount_RNA < 1800 &
                      percent.mt < 20)

######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 5000)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

dat <- merge(sets[[1]], list(sets[[2]], sets[[3]], sets[[4]], sets[[5]], sets[[6]]))
genes <- rownames(dat)

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat, features = rownames(dat))
dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 4000)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2)

dat <- RunUMAP(object = dat, dims = 1:50, min.dist = 0.5)
UMAPPlot(dat, group.by = 'Sample') #There's no strong batch effect!

######################################################################
#Embedd with dbMAP
######################################################################
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
umap <- reticulate::import('umap')
genes <- dat@assays$integrated@var.features
data <- t(dat@assays$integrated@data[genes,])
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(300), knn = as.integer(30))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 169 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

# Add to Seurat. For convenience, the UMAP layout of the diffusion components can be done with Seurat 'RunUMAP()'.
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

# Run dbMAP by using Seurat UMAP layout on the multiscaled and normalized diffusion components
dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.3, spread = 1.5, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

######################################################################
# Cluster, plot and save
######################################################################

dat <- FindNeighbors(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)))
dat <- FindClusters(dat, resolution = 0.8)

DimPlot(dat, reduction = 'dbmap', group.by = c('Sample', 'seurat_clusters'), pt.size = 0.5)
  
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE2', pt.size = 0.5, min.cutoff = 0.1) 
VlnPlot(dat, features = 'ACE2', group.by = 'seurat_clusters')

saveRDS(dat, 'HCL_Lung_Integrated.Rds')
