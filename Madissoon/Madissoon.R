######################################################################
#Analysis of Lung scRNAseq data from Madissoon et al
######################################################################

#Load libraries and processed data

reticulate::use_python('/usr/bin/python3')
library(reticulate)
library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(plotly)
library(cerebroApp)

# setwd("~/Madissoon")  #Change working directory to a dedicated folder

download.file('https://cellgeni.cog.sanger.ac.uk/tissue-stability/tissue-stability/lung_ts.rds', 
              destfile = '~/Madissoon/lung_ts.rds', method = 'wget')

dat1 <- readRDS('lung_ts.rds')
counts <- dat1@assays$RNA@counts
meta <- dat1@meta.data

dat <- CreateSeuratObject(counts, meta.data = meta)

######################################################################
#QC and filtering 
######################################################################

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-", assay = 'RNA')
counts_per_cell <- Matrix::colSums(dat)
counts_per_gene <- Matrix::rowSums(dat)
genes_per_cell <- Matrix::colSums(dat@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)', )
plot(sort(percent.mito), xlab='cell', log='y', main='percent.mito (ordered)')

VlnPlot(object = dat, features = "nCount_RNA", group.by = 'Donor')
VlnPlot(object = dat, features = "nCount_RNA", group.by = 'Time')
VlnPlot(object = dat, features = "nFeature_RNA", group.by = 'Donor')
VlnPlot(object = dat, features = "nFeature_RNA", group.by = 'Time')
VlnPlot(object = dat, features = "percent.mt", group.by = 'Donor')
VlnPlot(object = dat, features = "percent.mt", group.by = 'Time')

dat <- subset(dat, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
                nCount_RNA > 1800 & nCount_RNA < 35000 &
                percent.mt < 10)

######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
genes <- rownames(dat)
sets <- SplitObject(dat, split.by = 'Donor')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]])
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat, features = rownames(dat))
dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat)
dat <- RunUMAP(object = dat, dims = 1:50, min.dist = 0.5)
UMAPPlot(dat, group.by = 'Time') #There's no strong batch effect!
UMAPPlot(dat, group.by = 'Donor') #There's no strong batch effect!

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

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(200), knn = as.integer(15))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 147 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

# Add to Seurat. For convenience, the UMAP layout of the diffusion components can be done with Seurat 'RunUMAP()'.
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

# Run dbMAP by using Seurat UMAP layout on the multiscaled and normalized diffusion components
dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.6, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

######################################################################
# Cluster, plot and save
######################################################################

#Clustering
dat <- FindNeighbors(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings))) #Cluster on the diffusion structure information graph
dat <- FindClusters(dat)

#Some plots
DimPlot(dat, reduction = 'dbmap', group.by = 'Donor', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'Time', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'Celltypes', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5) 

#Check ACE2 expression
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE2', pt.size = 0.5, min.) 
VlnPlot(dat, features = 'ACE2', group.by = 'Celltypes')

# Save the Seurat file as RDS so that we can use it later
saveRDS(dat, 'Lung_Madissoon.Rds')
