######################################################################
#Analysis of Lung scRNAseq data (Travaglini et al)
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

# Set the working directory to the Travaglini et al dedicated folder
setwd("~/Travaglini") 

# Obtain the raw data. For this, you'll need to install synapser, and to create a synapse account (https://www.synapse.org/)
install.packages("synapser", repos=c("http://ran.synapse.org", "http://cran.fhcrc.org"))
library(synapser) 
synLogin('username','password') 

# Obtain synapse pointers and download the data 
syn21560510 <- synGet(entity= 'syn21560510',downloadLocation = 'data')  
syn21560409 <- synGet(entity= 'syn21560409', downloadLocation = 'data') 
 
counts10x <- as.matrix(read.csv('krasnow_hlca_10x_UMIs.csv', header = T, row.names = 1, stringsAsFactors = F))
meta10x <- read.csv('krasnow_hlca_10x_metadata.csv', header = T, row.names = 1, stringsAsFactors = F)

dat <- CreateSeuratObject(counts = counts10x, meta.data = meta10x)

######################################################################
#QC and filtering 
######################################################################

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-", assay = 'RNA')
percent.mito <- PercentageFeatureSet(dat, pattern = "^MT-", assay = 'RNA')
percent.mito <- percent.mito$nCount_RNA
dat <- AddMetaData(dat, metadata = percent.mito, col.name = 'percent.mt')
counts_per_cell <- Matrix::colSums(dat)
counts_per_gene <- Matrix::rowSums(dat)
genes_per_cell <- Matrix::colSums(dat@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')
plot(sort(percent.mito+1), xlab='cell', log='y', main='percent.mito (ordered)') #Empty!

VlnPlot(object = dat, features = "nCount_RNA", group.by = 'patient')
VlnPlot(object = dat, features = "nCount_RNA", group.by = 'channel')
VlnPlot(object = dat, features = "nCount_RNA", group.by = 'sample')
VlnPlot(object = dat, features = "nFeature_RNA", group.by = 'patient')
VlnPlot(object = dat, features = "nFeature_RNA", group.by = 'sample')

dat <- subset(dat, subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
                nCount_RNA < 40000)

######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
genes <- rownames(dat)
sets <- SplitObject(dat, split.by = 'sample')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]])
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3500)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)

rm(dat, sets)
gc()

dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)

dat <- ScaleData(dat, features = rownames(dat))
dat <- FindVariableFeatures(dat, selection.method = 'disp')

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2)

dat <- RunUMAP(object = dat, dims = 1:30, min.dist = 0.5)
UMAPPlot(dat, group.by = 'patient') #There's no strong batch effect
UMAPPlot(dat, group.by = 'sample') #There's no strong batch effect
UMAPPlot(dat, group.by = 'channel') #There's no strong batch effect
UMAPPlot(dat, group.by = 'free_annotation') 

FeaturePlot(dat, features = 'ACE2', order = T, min.cutoff = 0, max.cutoff = 0.1, pt.size = 0.5)
FeaturePlot(dat, features = 'KNG1', order = T, min.cutoff = 0, max.cutoff = 0.2, pt.size = 0.5)
FeaturePlot(dat, features = c('ACE2', 'KNG1'), order = T, blend = T, min.cutoff = 0, max.cutoff = 0.2)

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
plot(evals) #Select meaningful diffusion components. Used 221 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

#Add to Seurat
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

dat <- FindNeighbors(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 0.8, graph.name = 'dbgraph')

######################################################################
#dbMAP and Adjustment
######################################################################
dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.3, spread = 2, learning.rate = 2, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'db', n.components = 3, dims = 1:ncol(dat@reductions$db@cell.embeddings), min.dist = 0.3, spread = 2, learning.rate = 2, reduction.key = 'dbMAP3D_', reduction.name = 'dbmap3d', init = 'spectral')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5)

DimPlot(dat, reduction = 'dbmap', group.by = 'Sample', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'free_annotation', pt.size = 0.5, label = T, repel = F)
FeaturePlot(dat, features = 'ACE2', reduction = 'dbmap', order = T, min.cutoff = 0, max.cutoff = 0.1)
FeaturePlot(dat, features = 'KNG1', reduction = 'dbmap', order = T, min.cutoff = 0, max.cutoff = 0.001)

plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "sample", 'free_annotation', 'ACE2', 'KNG1'
))
plot.data$label <- dat$free_annotation
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~KNG1, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=2),
        text=~free_annotation,
        hoverinfo="text") 

saveRDS(dat, 'Travaglini_integrated.Rds')
