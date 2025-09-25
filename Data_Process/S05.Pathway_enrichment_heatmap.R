library("magrittr")
library("sscVis")
library("R.utils")
library("ggpubr")
library("ggplot2")
library("plyr")
library("grid")
library("cowplot")
library("ggrepel")
library("data.table")
library("tidyverse")
library(viridis)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
######################
setwd('H:\\Reng')
Score_matrix <- read.csv('./Gset_analysis_score.csv')
Pathway_meta <- read.csv('./Pathway.csv')
Cluster_num <- read.csv('./Cluster_num.csv')
####################### Adjust
colnames(Score_matrix)[1] <- 'Metabolic'
colnames(Pathway_meta)[1] <- 'Metabolic'
Score_matrix$Type <- Pathway_meta$Classification[match(Score_matrix$Metabolic, Pathway_meta$Metabolic)]
Score_matrix <- Score_matrix[order(Score_matrix$Type),]
Plot_matrix <- Score_matrix
expre <- Plot_matrix[,2:61]
expre <- t(expre) %>% as.data.frame()
expre$clust <- Cluster_num$clust
expre <- expre[order(expre$clust),]
expre <- expre[,1:86]
expre <- t(expre) %>% as.data.frame()
Plot_matrix <- cbind(Plot_matrix$Type, expre)
colnames(Plot_matrix)[1] <- 'Type'
############################################## Plot
Row_split <-  Plot_matrix$Type
Order_cluster <- Cluster_num[order(Cluster_num$clust),'clust']
Order_cluster <- as.character(Order_cluster)
Col_split <- Order_cluster

expre <- scale(expre,center = TRUE)

Color1 = colorRamp2(seq(-3, 3, length = 10), colorRampPalette(c("blue","white","red"))(10))

p3 <- Heatmap(expre, name = "Expression", col = Color1, row_split = Row_split, column_names_rot = 45,
              column_split = Col_split, cluster_rows = TRUE,cluster_columns = TRUE, show_column_names = TRUE,show_row_names = TRUE, #rect_gp = gpar(col = "black", lwd = 2),
              cluster_column_slices = FALSE,row_title_rot = 0)










