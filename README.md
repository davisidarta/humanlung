[![License: GPL-3.0](https://img.shields.io/badge/License-GNU--GLP%20v3.0-green.svg)](https://opensource.org/licenses/GPL-3.0)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/DaviSidarta.svg?label=Follow%20%40DaviSidarta&style=social)](https://twitter.com/DaviSidarta)


# Human Lung Integrated Cell Atlas
  Code for the human lung integrated cell atlas generation in Sidarta-Oliveira et al. In this work, we performed a meta-analysis of single-cell RNA sequencing of the human lung to achieve insights into ACE2 cellular role in the lung. For this, we analyzed each dataset individually with Seurat v3(), and merged them into an comprehensive integrated dataset that serves as a reference for further single-cell lung studies. We provide the complete Seurat object in (), as well as a loom file (), which is less computationally demanding. 
    ![Human Lung Integrated Cell Atlas](https://github.com/davisidarta/humanlung/blob/master/Lung.png)


# Interactive data exploration with Cerebro
  We provide the atlas as a fully interactive Cerebro object at (link_atlas). This atlas can be explored by non-biologists to easily generate publication-level figures from lung gene expression, and also allow for searching cluster marker gene expression and functional annotations from a variey of databases (GO, KEGG, etc) which can be exported as tables (.csv, opens in Excel and other table or text readers). If you have no experience with coding, we suggest you explore the data in the atlas or into a local session of Cerebro, which can be installed here (link_cerebro). One you've installed Cerebro, you can explore the lung atlas data with this (cerebro_file_link) object. 

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
## 3- Setup the directories and run the analysis
  
  Please clone and extract this repository. From the repository home folder, run:
  
```  
  R < Madissoon/Madissoon.R
  R < Reyfman/Reyfman.R
  R < HCL/HCL_Lung.R
  R < Integration.R
 ```  
  This reproduces all the results. Alternatively, users can execute these R scripts in the provided order on an R IDE such as RStudio.
 
# Citation
Sidarta-Oliveira, D. et al. 

# Troubleshooting
  For technical issues, please use the 'issues' section of this repository or contact davisidarta@fcm.unicamp.br. 
