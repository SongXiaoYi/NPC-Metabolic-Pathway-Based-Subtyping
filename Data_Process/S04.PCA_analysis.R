#加载包
#####################################
library(ggplot2)
library(scales)
library(jjAnno)
library(dplyr)
library(ggalluvial)
library(patchwork)
library(Seurat)
library(data.table)
library(Matrix)
library(dendextend)
library(tidyr)
library(reshape2)
library(magrittr)
library(tidyverse)
library(viridis)   #色盲色
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(dendsort)
################################
library(ggpubr)
library(limma)
library(magrittr)
library(destiny)
library(factoextra)
library(FactoMineR)
###########################################
###########################
'%!in%' <- Negate('%in%')   #### Not in function
#############################
setwd('H:\\Reng')
pbmc <- read.csv('./Gset_analysis_score.csv',header = T)
rownames(pbmc) <- pbmc$X
pbmc <- pbmc[,-1]
################# Cluster
Cluster_meta <- read.csv('Cluster_num.csv')
pbmc <- t(pbmc)
data_table <- as.data.frame(pbmc)
data_table$Patient_Match <- rownames(data_table)
data_table$Cluster <- as.character(Cluster_meta$clust)
#############################################
########################################
########################## Palette
Type_colors <- c('red', 'blue')   ######### 30个月
ROI_colors <- c('#FF4500', '#0000FF', '#008000')

res.pca <- PCA(data_table[,1:86],scale.unit = T,ncp=5,graph = T)


########################## Treatment
PCA_plot <- fviz_pca_ind(res.pca,
                  geom.ind = "point", pointsize = 2, 
                  col.ind = data_table$Cluster, # color by groups
                  palette = ROI_colors,
                  addEllipses = TRUE, # Concentration ellipses
                  legend.title = "Cluster",mean.point=F) 

p1 <- PCA_plot +  theme(axis.text.y= element_text(size=15,color="black", vjust=0.5, hjust=0.5)) + theme(axis.text.x= element_text(size=15,color="black", vjust=0.5, hjust=0.5)) + theme(axis.title= element_text(size=15, color="black", vjust=0.5, hjust=0.5)) + theme(legend.key.width = unit(1, "cm"),legend.key.height = unit(1, "cm")) + theme(legend.text= element_text(size=15,color="black",vjust=0.5, hjust=0.5)) + theme(legend.title=element_text(size=14)) + labs(title='Sample individuals - PCA') + scale_shape_manual(values = c(15, 19, 17, 18)) +
ggtitle('Cluster') + theme(plot.title = element_text(hjust = 0.5))

setwd('H:\\Reng')
pdf(file = 'PCA_analysis.pdf', width = 4.5, height = 3.65)
p1
dev.off()
