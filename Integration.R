######################################################################
#Integrative analysis of Lung scRNAseq data for the manuscript analyzing possible mechanisms of COVID patogenesis
######################################################################
#Load libraries and processed data
reticulate::use_python('/usr/bin/python3')
library(reticulate)
library(Matrix)
library(Seurat)
library(plotly)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Lung")

######################################################################
#Import data
######################################################################

M <- readRDS('Madissoon/Lung_analyzed.Rds')
Y <- readRDS('Reyfman/Reyfman_integrated.Rds')
H <- readRDS('HCL/HCL_Lung_Integrated.Rds')

DefaultAssay(M) <- 'integrated'
DefaultAssay(Y) <- 'integrated'
DefaultAssay(H) <- 'integrated'

M <- DietSeurat(M, assays = DefaultAssay(M), counts = F,  data = T, scale.data = F)
Y <- DietSeurat(Y, assays = DefaultAssay(Y), counts = F,  data = T, scale.data = F)
H <- DietSeurat(H, assays = DefaultAssay(H), counts = F,  data = T, scale.data = F)

#Add some metadata
sm <- rep('Madisson & Willbrey-Clark et al', times = ncol(M))
sy <- rep('Reyfman et al', times = ncol(Y))
sh <- rep('Human Cell Landscape', times = ncol(H))
sample_label <- c( sm, sy, sh)

dat <- merge(M, list(Y, H))
dat <- AddMetaData(dat, metadata = sample_label, col.name = 'Study')
dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 20000)
genes <- VariableFeatures(dat) #Because integration using all genes is expensive, let's use only the top 20,000 variable genes
'ACE2' %in% genes #Check that ACE2 is contained in this gene list
'KNG1' %in% genes #Check that KNG1 is contained in this gene list
genes <- append(genes, 'KNG1') #Append KNG1 expression, as the KKS is of interest

######################################################################
#Merge with CCA
######################################################################
sets <-SplitObject(dat, split.by = 'Study')
rm(dat, H, M, Y) # let's preserve our memory, we'll need it for integration.
gc()
feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000, assay = c('integrated', 'integrated', 'integrated'))
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)

dat <- ScaleData(dat)
dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2)

dat <- RunUMAP(object = dat, dims = 1:50, min.dist = 0.5)

UMAPPlot(dat, group.by = 'Study') # There's no strong batch effect among studies
        
######################################################################
# Use Travaglini et al as a training set for FindTransferAnchors to help with cell subtype annotation.
######################################################################

Tr <- readRDS('Travaglini/Travaglini_integrated_diet.Rds')
Tr <- AddMetaData(Tr, metadata = rep('Travaglini et al'), times = ncol(Tr), col.name = 'Study')
Tr <- FindVariableFeatures(Tr, selection.method = 'disp', nfeatures = 5000)

######################################################################
#Embedd Travaglini et al data with dbMAP
######################################################################
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
umap <- reticulate::import('umap')
genes <- Tr@assays$integrated@var.features
data <- t(Tr@assays$integrated@data[genes,])
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(300), knn = as.integer(50))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 191 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

#Add to Seurat
Tr@reductions$db <- Tr@reductions$pca
rownames(db) <- colnames(Tr)
Tr@reductions$db@cell.embeddings <- db

######################################################################
#Embedd Atlas data with dbMAP
######################################################################
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
umap <- reticulate::import('umap')

dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 5000)

genes <- dat@assays$integrated@var.features
data <- t(dat@assays$integrated@data[genes,])
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(300), knn = as.integer(50))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 169 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

#Add to Seurat
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

######################################################################
#Find anchors for annotation
######################################################################

anchors <- FindTransferAnchors(reference = Tr, query = dat, 
                               dims = 1:30, reduction = 'cca')
predictions <- TransferData(anchorset = anchors, refdata = Tr$free_annotation, 
                            dims = 1:30, weight.reduction = 'cca')


dat <- AddMetaData(dat, metadata = predictions)


######################################################################
#Embedd with dbMAP diffusion structure
######################################################################
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
umap <- reticulate::import('umap')
genes <- dat@assays$integrated@var.features
data <- t(dat@assays$integrated@data[genes,])
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(300), knn = as.integer(50))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 169 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

# Add to Seurat. For convenience, the UMAP layout of the diffusion components can be done with Seurat 'RunUMAP()'.
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

dat <- FindNeighbors(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 0.8, graph.name = 'dbgraph')
dat$db.cluster <- dat$seurat_clusters

######################################################################
#dbMAP and Adjustment
######################################################################

# Run dbMAP by using Seurat UMAP layout on the multiscaled and normalized diffusion components

dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.3, spread = 2, learning.rate = 2, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'db', n.components = 3, dims = 1:ncol(dat@reductions$db@cell.embeddings), min.dist = 0.3, spread = 2, learning.rate = 2, reduction.key = 'dbMAP3D_', reduction.name = 'dbmap3d', init = 'spectral')

######################################################################
#Plot
######################################################################

DimPlot(dat, reduction = 'dbmap', group.by = 'Study', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'predicted.id', pt.size = 0.5) #Labels predicted from Travaglini et al data
DimPlot(dat, reduction = 'dbmap', group.by = 'Celltypes', pt.size = 0.5) #Original metadata from Madissoon et al
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5) + DarkTheme() + NoLegend()
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.2, order = T) 

plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "Donor", 'Time', 'Celltypes', 'seurat_clusters', 
                                              'ACE2', 'KNG1', 'BDKRB2', 'BDKRB1', 'IL6'
))
plot.data$label <- dat$seurat_clusters
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~seurat_clusters, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~seurat_clusters,
        hoverinfo="text", plot_bgcolor = 'black') 

######################################################################
#Now thatwe've annotated the dataset, lets change cell type labels
######################################################################

current.cluster.ids <- levels(Idents(dat))

new.ids <-  c("Alveolar type 2", 
              "Macrophage 1", 
              "Alveolar type 2 Basal",
              "T CD4+ Memory/Effector",
              "Macrophage 2",
              'T CD8+ Memory/Effector',
              "Macrophage 3",
              'NK 1',
              "Alveolar type 2.3",
              'T CD4+ Naive',
              'Fibroblast/Perycyte',
              'Monocytes 1',
              'Monocytes 2',
              'Mast',
              'Endothelial Capillary',
              'Alveolar type 2 Club',
              'NK 2',
              'Endothelial Vessel',
              'Macrophage Proliferating',
              'Signaling Alveolar type 2.1',
              'DC 2',
              'Alveolar type 1',
              'Doublets 1',
              'Signaling Alveolar type 2.2',
              'Alveolar type 2/1 Intermediate',
              'B Cells',
              'Alveolar Fibroblast',
              'Monocytes 3',
              'Muscle',
              'Alveolar type 1 Basal',
              'Ciliated',
              'Plasma',
              'Macrophage 5',
              'T CD8 2',
              'Macrophage 6',
              'Lymph vessel',
              'NK dividing',
              'DC 1',
              'Macrophage 7',
              'Macrophage 8',
              'T CD4/CD8',
              'Doublets 2',
              'DC/Mono dividing',
              'DC activated',
              'NK 3',
              'T CD8 3',
              'T CD4 3',
              'Alveolar/Ciliated',
              'T CD4 NA',
              'Unassigned1',
              'Unassigned2',
              'Unassigned3'
              )


names(new.ids) <- levels(dat)
dat <- RenameIdents(dat, new.ids)
DimPlot(dat, label = T, repel = T, label.size = 3.8, reduction = 'dbmap')

dat@meta.data$seurat_clusters <- plyr::mapvalues(x = dat@meta.data$seurat_clusters, from = current.cluster.ids, to = new.ids)

######################################################################
#Remove doublets
######################################################################
Idents(dat) <- 'seurat_clusters'
dat <- subset(dat, cells = WhichCells(dat, idents = c('Doublets 1', 'Doublets 2', 'Doublets 3', 'Doublets 4', 'Doublets 5'), invert = T))

# Final dbMAP plot - Figure 2B
DimPlot(dat, label = T, repel = T, label.size = 6, reduction = 'dbmap')

# Learned cell types from Travaglini et al - Supplementary Figure 2C
DimPlot(dat, label = T, repel = T, label.size = 6, reduction = 'dbmap', group.by = 'predicted.id') + NoLegend() 

# Learned cell types from Travaglini et al - Supplementary Figure 2D
DimPlot(dat, label = T, repel = T, label.size = 6, reduction = 'umap')

######################################################################
# Plots for figures 3, 4, 5 and 6
######################################################################
# ACE2
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.2, order = T) 

#DimPlots for Figure 3
FeaturePlot(dat, reduction = 'dbmap', features = 'CTSL', pt.size = 0.5, min.cutoff = 0, max.cutoff = 3.5, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'TMPRSS2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'TPCN2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.5, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'PIKFYVE', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.5, order = T) 
#DotPlot for Figure 3
DotPlot(dat, features = c('ACE2','TMPRSS2','TPCN2','PIKFYVE','CTSL'), scale.by = 'radius', col.min = 0) + RotatedAxis()

#DimPlots for Figure 4
FeaturePlot(dat, reduction = 'dbmap', features = 'KNG1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.02, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'KLKB1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'BDKRB2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.5, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'BDKRB1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.5, order = T) 
#DotPlot for Figure 4
DotPlot(dat, features = c('KNG1','KLKB1','BDKRB1','ACE','BDKRB2', 'ACE2'), col.min = 0, col.max = 6) + RotatedAxis()

#DimPlots for Figure 5
FeaturePlot(dat, reduction = 'dbmap', features = 'AGT', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'REN', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'ACE', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'AGTR1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.5, order = T) 
#DotPlot for Figure 5
DotPlot(dat, features = c('ACE2','REN','AGT','ACE', 'AGTR1'), col.min = 0, col.max = 3) + RotatedAxis()

#DimPlots for Figure 6
FeaturePlot(dat, reduction = 'dbmap', features = 'KLKB1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'SERPINE1', pt.size = 0.5, min.cutoff = 0, max.cutoff = 2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'PLAT', pt.size = 0.5, min.cutoff = 0, max.cutoff = 2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'FGG', pt.size = 0.5, min.cutoff = 0, max.cutoff = 2, order = T) 
#DotPlot for Figure 6
DotPlot(dat, features = c('ACE2','PLAT','FGG','SERPINE1', 'KLKB1'), col.min = 0, col.max = 5) + RotatedAxis()

######################################################################
# Marker genes - Supplementary Figures 3, 4 and 5
######################################################################

markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5, assay = 'integrated')
write.csv(markers, 'Lung_markers.csv', row.names = T, col.names = T)
top2 <- markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Heatmap - Supplementary Figure 5
DoHeatmap(subset(dat, downsample = 100), features = top2$gene, assay = 'integrated', size = 3) + NoLegend()

# Dotplot - Supplementary Figure 3
genes <- top$gene
DotPlot(dat, features = make.unique(genes), group.by = 'seurat_clusters', scale.by = 'radius') + RotatedAxis() + FontSize(main = 1)

# Gene expression panel - Supplementary Figure 4
FeaturePlot(dat, reduction = 'dbmap', features = c('SFTPC', 'C1QA', 'FCN1', 'GZMH', 'CLDN5', 
                                                   'MS4A1', 'IGLC2','LUM', 'TFF3', 'TPSB2'),
            pt.size = 0.01, min.cutoff = c(rep(0, times = 10)),
            max.cutoff = c(14, 14, 14, 14 , 14, 14, 14, 14, 14), order = T, ncol = 2, 
            combine = T) 

######################################################################
#Export to Cerebro
######################################################################
# Let's allow Cerebro to deal with corrected, normalized data
dat@assays$integrated@counts <- dat@assays$integrated@data

#Add percentage of mitochondrial and ribossomal genes
dat <- addPercentMtRibo(dat, organism = 'hg', gene_nomenclature = 'name', assay = 'integrated')

#Add marker genes for each cluster and study
dat <- getMarkerGenes(dat, organism = 'hg', column_sample = 'Study', column_cluster = 'seurat_clusters', assay = 'integrated', min_pct = 0.7)

# Perform functional enrichment analysis for each cluster and sample. We perform a broad look over a series of databases
dat <- getEnrichedPathways(dat, column_sample = 'db.cluster', column_cluster = 'seurat_clusters',
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))

# Export the cerebro (.crb) file
cerebro <- exportFromSeurat(dat, file = 'Lung_integrated.crb', experiment_name = 'Healthy Human Lung', organism = 'hg',
                            column_sample = 'Study', column_cluster = 'seurat_clusters',
                            column_nUMI = 'nCount_RNA', column_nGene = 'nFeature_RNA', assay = 'integrated')
# Save the final Seurat object with functional enrichment analysis
saveRDS(dat, 'Lung_analyzed.Rds')
launchCerebro(maxFileSize = 500000)

