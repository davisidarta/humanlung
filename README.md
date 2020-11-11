[![License: GPL-3.0](https://img.shields.io/badge/License-GNU--GLP%20v3.0-green.svg)](https://opensource.org/licenses/GPL-3.0)
[![Twitter](https://img.shields.io/twitter/url/https/twitter.com/DaviSidarta.svg?label=Follow%20%40DaviSidarta&style=social)](https://twitter.com/DaviSidarta)


# Human Lung Integrated Cell Atlas
  Code for the human lung integrated cell atlas generation in [Sidarta-Oliveira et al](https://www.nature.com/articles/s41598-020-76488-2). In this work, we performed a meta-analysis of single-cell RNA sequencing of the human lung to achieve insights into ACE2 cellular role in the lung. For this, we analyzed each dataset individually with [Seurat v3](https://github.com/satijalab/seurat), and merged them into an comprehensive integrated dataset that serves as a reference for further single-cell lung studies. 
    ![Human Lung Integrated Cell Atlas](https://github.com/davisidarta/humanlung/blob/master/Lung.png)


# Interactive data exploration with Cerebro
  We provide a subset of 10,000 randomly sampled cells from the atlas as an online interactive Cerebro object at https://humanlung.iqm.unicamp.br . This atlas can be explored by non-biologists to easily generate publication-level figures from lung gene expression, and also allow for searching cluster marker gene expression and functional annotations from a variey of databases (GO, KEGG, etc) which can be exported as tables (.csv, opens in Excel and other table or text readers). If you have no experience with coding, we suggest you explore the data in the atlas or into a local session of Cerebro, which can be installed  at https://github.com/romanhaa/Cerebro. One you've installed Cerebro, you can explore the lung atlas data with this [object](https://figshare.com/s/169a53cccfe5b341d1fb) 

# Docker 
  To make the reviewing and reproduction process easier, we provide a docker image which can be used to reproduce our findings in a high-performance computing cluster in a containerized fashion. This is particularly advised if you do not have rights to install R and python packages in your machine. If you choose to use the docker image, you can skip steps 1 and 2 of the Reproduction Steps. If you do not have docker installed, check this [page](https://docs.docker.com/get-docker/).
  Once docker is installed, run on a terminal:
  
  ```
  docker pull davisidarta/humanlung
  ```
  
  # Get docker image

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
Sidarta-Oliveira, D., Jara, C.P., Ferruzzi, A.J. et al. SARS-CoV-2 receptor is co-expressed with elements of the kinin–kallikrein, renin–angiotensin and coagulation systems in alveolar cells. Sci Rep 10, 19522 (2020). https://doi.org/10.1038/s41598-020-76488-2

# Troubleshooting
  For technical issues, please use the 'issues' section of this repository or contact davisidarta@fcm.unicamp.br. 
