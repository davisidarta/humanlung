######################################################################
#Analysis of Lung scRNAseq data from Reyfman et al
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

# Set the working directory to the Reyfman et al dedicated folder
setwd("~/Reyfman") 

# Obtain the raw data
library(utils)
download.file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122960&format=file', 
              destfile = 'GSE122960.tar', method = 'wget')
untar(tarfile = 'GSE122960.tar', exdir = "data")

######################################################################
#Load data
######################################################################

counts1 <- Read10X_h5('data/GSM3489182_Donor_01_raw_gene_bc_matrices_h5.h5')
counts2 <- Read10X_h5('data/GSM3489185_Donor_02_raw_gene_bc_matrices_h5.h5')
counts3 <- Read10X_h5('data/GSM3489187_Donor_03_raw_gene_bc_matrices_h5.h5')
counts4 <- Read10X_h5('data/GSM3489189_Donor_04_raw_gene_bc_matrices_h5.h5')
counts5 <- Read10X_h5('data/GSM3489191_Donor_05_raw_gene_bc_matrices_h5.h5')
counts6 <- Read10X_h5('data/GSM3489193_Donor_06_raw_gene_bc_matrices_h5.h5')
counts7 <- Read10X_h5('data/GSM3489195_Donor_07_raw_gene_bc_matrices_h5.h5')
counts8 <- Read10X_h5('data/GSM3489197_Donor_08_raw_gene_bc_matrices_h5.h5')

dat1 <- CreateSeuratObject(counts1)
dat2 <- CreateSeuratObject(counts2)
dat3 <- CreateSeuratObject(counts3)
dat4 <- CreateSeuratObject(counts4)
dat5 <- CreateSeuratObject(counts5)
dat6 <- CreateSeuratObject(counts6)
dat7 <- CreateSeuratObject(counts7)
dat8 <- CreateSeuratObject(counts8)

dat <- merge(dat1, list(dat2, dat3, dat4, dat5, dat6, dat7, dat8))
#Add some metadata
s1 <- rep('Sample 1', times = ncol(dat1))
s2 <- rep('Sample 2', times = ncol(dat2))
s3 <- rep('Sample 3', times = ncol(dat3))
s4 <- rep('Sample 4', times = ncol(dat4))
s5 <- rep('Sample 5', times = ncol(dat5))
s6 <- rep('Sample 6', times = ncol(dat6))
s7 <- rep('Sample 7', times = ncol(dat7))
s8 <- rep('Sample 8', times = ncol(dat8))
sample_label <- c(s1, s2, s3, s4, s5, s6, s7, s8)
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
} #Filter very low quality cells
sets[[1]] <- subset(sets[[1]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[2]] <- subset(sets[[2]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[3]] <- subset(sets[[3]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[4]] <- subset(sets[[4]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[5]] <- subset(sets[[5]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[6]] <- subset(sets[[6]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[7]] <- subset(sets[[7]], subset = nFeature_RNA > 300 & nCount_RNA > 300)
sets[[8]] <- subset(sets[[8]], subset = nFeature_RNA > 300 & nCount_RNA > 300)

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

sets[[1]] <- subset(sets[[1]], subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &
                                        nCount_RNA > 500 & nCount_RNA < 15000 &
                                        percent.mito < 20)
sets[[2]] <- subset(sets[[2]], subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &
                                        nCount_RNA > 1000 & nCount_RNA < 15000 &
                                        percent.mito < 10)
sets[[3]] <- subset(sets[[3]], subset = nFeature_RNA > 300 & nFeature_RNA < 2500 &
                                        nCount_RNA > 500 & nCount_RNA < 10000 &
                                        percent.mito < 10)
sets[[4]] <- subset(sets[[4]], subset = nFeature_RNA > 300 & nFeature_RNA < 3500 &
                                        nCount_RNA > 500 & nCount_RNA < 10000 &
                                        percent.mito < 10)
sets[[5]] <- subset(sets[[5]], subset = nFeature_RNA > 300 & nFeature_RNA < 4500 &
                                        nCount_RNA > 500 & nCount_RNA < 20000 &
                                        percent.mito < 10)
sets[[6]] <- subset(sets[[6]], subset = nFeature_RNA > 300 & nFeature_RNA < 5000 &
                                        nCount_RNA > 500 & nCount_RNA < 15000 &
                                        percent.mito < 10)
sets[[7]] <- subset(sets[[7]], subset = nFeature_RNA > 300 & nFeature_RNA < 4500 &
                                        nCount_RNA > 500 & nCount_RNA < 15000 &
                                        percent.mito < 10)
sets[[8]] <- subset(sets[[8]], subset = nFeature_RNA > 500 & nFeature_RNA < 1000 &
                                        nCount_RNA > 1000 & nCount_RNA < 2000 &
                                        percent.mito < 10)


######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
dat <- merge(sets[[1]], list(sets[[2]], sets[[3]], sets[[4]], sets[[5]], sets[[6]], sets[[7]], sets[[8]]))
sets <- SplitObject(dat, split.by = 'Sample')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 5000)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

dat <- merge(sets[[1]], list(sets[[2]], sets[[3]], sets[[4]], sets[[5]], sets[[6]], sets[[7]], sets[[8]]))
genes <- rownames(dat)

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat, features = rownames(dat))


######################################################################
#Default workflow
######################################################################

dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 4000)
VariableFeaturePlot(dat)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2)

dat <- RunUMAP(object = dat, dims = 1:30, min.dist = 0.5)
UMAPPlot(dat, group.by = 'Sample') #There's no strong batch effect towards Sample 8!

plot <- UMAPPlot(dat, group.by = 'Sample') #There's no strong batch effect towards Sample 8!
HoverLocator(plot, information = FetchData(dat, vars = 'Sample'))

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
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE2', pt.size = 0.5, min.cutoff = 0.2) 
VlnPlot(dat, features = 'ACE2', group.by = 'seurat_clusters')

#Diet it down and save
dat <- DietSeurat(dat, assays = DefaultAssay(dat), counts = F,  data = T, scale.data = F)
saveRDS(dat, 'Reyfman_integrated.Rds')
