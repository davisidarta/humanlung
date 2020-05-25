# Human Lung Integrated Cell Atlas
  Code for the human lung integrated cell atlas generation in Sidarta-Oliveira et al. In this work, we performed a meta-analysis of single-cell RNA sequencing of the human lung to perform insights into ACE2 cell function. For this, we analyzed each dataset individually with Seurat v3(), and merged them into an comprehensive integrated dataset that serves as a reference for further single-cell lung studies. We provide the complete Seurat object in (), as well as a loom file (), which is less computationally demanding. 

# Interactive exploration of data with Cerebro
  We provide the atlas as a fully interactive Cerebro object at (link_atlas). This atlas can be explored by non-biologists to easily generate publication-level figures from lung gene expression, and also allow for searching cluster marker gene expression, functional annotations from a variey of databases (GO, KEGG, etc). If you have no experience with coding, we suggest we explore the data in the atlas or into a local session of Cerebro, which can be installed here (link_cerebro). 

# Docker 
  To make the reviewing and reproduction process easier, we provide a docker image (link_imagem docker) which can be used to reproduce our findings in a high-performance computing cluster in a containerized fashion. This is particularly advised if you do not have rights to install R and python packages in your machine. If you choose to use the docker image, you can skip steps 1 and 2 of the Reproduction Steps.

# Reproduction Steps
  Please follow these steps on a high-performance computing cluster, cloud server or server-grade workstation:
  
## 1- Install dbMAP
  Prior to installing dbMAP, make sure you have scikit-build and cmake available in your system. These are required for installation.
  ```
  sudo apt-get install cmake
  pip3 install scikit-build
  ```
dbMAP has been implemented in Python3, and can be installed using:

```
git clone git://github.com/davisidarta/dbMAP.git
cd dbMAP
sudo -H pip3 install .
```
## 2- Install Seurat
  In a R terminal, notebook or a RStudio-based session, install Seurat. We'll need it for data integration and as a general toolkit for data analysis.

```
install.packages('Seurat')
```
## 3- Obtain data
  As detailed in our Methods section (link_preprint), we obtained data from the following sources:
  ```
  # Reyfman et al. data
  
  # Human Cell Atlas
  
  # Madissoon et al. data
  
  # Travaglini et al. data
  
  ```
  
  Run each study under the 'Preprocess' directory to analyze all data with a batch-correction of intra-study sample heterogeneity. This usually takes a long time and should be done one study at a time. Be sure to have disk space to store the raw data and the processed Seurat objects.
  ```
  R (ver como fazer isso em R)
  
  ```
  
## 4- Reproduce computations
  To reproduce data integration and label transfering, run the scripts under the 'process' directory.
  ```
  R (ver como fazer isso em R
  ```

## 5- Reproduce plots

  To generate high-resolution images that reproduce all the plots shown in the manuscript, please run the `plots.R` script.

# Citation
Sidarta-Oliveira, D. et al. Citation.
Reyfman et al. citation
Madissoon et al. citation.
Human Cell Landscape citation

# Troubleshooting
  For technical issues, please use the 'issues' section of this repository. For more detailed questions, please contact davisidarta@fcm.unicamp.br .
