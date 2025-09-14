########## Pairwise_spatial_enrichment_results_transfer
library(data.table)
library(patchwork)
library(Matrix)
library(Seurat)
library(dendextend)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(tidyverse)
library(viridis)   #色盲色
library(stringr)
library(tibble)
library(tidyr)
library(reshape2)
library(limma)
library(MOVICS)
########################################################################
'%!in%' <- Negate('%in%')   #### Not in function
############################ Here begin
##################################################### Tumor
setwd('H:\\Reng')
pbmc <- read.csv('./Gset_analysis_score.csv')














