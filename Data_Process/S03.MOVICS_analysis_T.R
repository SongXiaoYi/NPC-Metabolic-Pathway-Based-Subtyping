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
pbmc <- read.csv('./Gset_analysis_score.csv',header = T)
rownames(pbmc) <- pbmc$X
pbmc <- pbmc[,-1]
pbmc1 <- pbmc
pbmc2 <- pbmc
#################################
########################################################
mo.data <- list(pbmc1 = pbmc1, pbmc2 = pbmc2)
############################################################################## MOVICS analysis
optk <- getClustNum(data = mo.data,
                    is.binary   = c(F,F), # note: the 4th data is somatic mutation which is a binary matrix
                    try.N.clust = 2:8, # try cluster number from 2 to 8
                    fig.name    = "CLUSTER NUMBER OF Allregion")
###############################################################################
#moic.res.list <- getMOIC(data = mo.data,
#                         methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "ConsensusClustering", "IntNMF", "iClusterBayes"),
#                         N.clust = 3,
#                         type = c("gaussian", "gaussian", "gaussian", "binomial"))

moic.SNF <- getMOIC(data = mo.data,
                         methodslist = "SNF",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))


moic.CIMLR <- getMOIC(data = mo.data,
                         methodslist = "CIMLR",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.PINSPlus <- getMOIC(data = mo.data,
                         methodslist = "PINSPlus",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.NEMO <- getMOIC(data = mo.data,
                         methodslist = "NEMO",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.COCA <- getMOIC(data = mo.data,
                         methodslist = "COCA",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.MoCluster <- getMOIC(data = mo.data,
                         methodslist = "MoCluster",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.LRAcluster <- getMOIC(data = mo.data,
                         methodslist = "LRAcluster",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.ConsensusClustering <- getMOIC(data = mo.data,
                         methodslist = "ConsensusClustering",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.IntNMF <- getMOIC(data = mo.data,
                         methodslist = "IntNMF",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))

moic.iClusterBayes <- getMOIC(data = mo.data,
                         methodslist = "iClusterBayes",
                         N.clust = 3,
                         type = c("gaussian", "gaussian"))


moic.res.list <- list(moic.SNF, moic.CIMLR, moic.PINSPlus, moic.NEMO, moic.COCA, moic.MoCluster, moic.LRAcluster, moic.ConsensusClustering, moic.IntNMF, moic.iClusterBayes)
names(moic.res.list) <- c("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "ConsensusClustering", "IntNMF", "iClusterBayes")

cmoic <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "ward.D")    # average

############################################# Silhouette
getSilhouette(sil      = cmoic$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 4,
              width    = 3.55)
############################################# Plot
indata <- mo.data
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2), # no truncation for mutation
                     centerFlag = c(T,T), # no center for mutation
                     scaleFlag  = c(T,T)) # no scale for mutation

feat1 <- rownames(pbmc1)
feat2 <- rownames(pbmc2)
annRow <- list(feat1, feat2)

Freq.col <- c("#6699CC", "white", "#FF3C38")
Density.col <- c("#6699CC", "white", "#FF3C38")
#Districts.col <- c("#6699CC", "white", "#FF3C38")
col.list   <- list(Freq.col, Density.col)

getMoHeatmap(data          = plotdata,
             row.title     = c("1","2"),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("Enrichment","Enrichment"),
             clust.res     = cmoic$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

setwd('H:\\Reng')
write.csv(cmoic$clust.res, file = 'Cluster_num.csv')
######################################################
Combine_clust <-  data.frame( SNF = as.factor(moic.SNF$clust.res[,2]), 
                              CIMLR = as.factor(moic.CIMLR$clust.res[,2]),
                              PINSPlus = as.factor(moic.PINSPlus$clust.res[,2]),
                              NEMO = as.factor(moic.NEMO$clust.res[,2]),
                              COCA = as.factor(moic.COCA$clust.res[,2]),
                              MoCluster = as.factor(moic.MoCluster$clust.res[,2]),
                              LRAcluster = as.factor(moic.LRAcluster$clust.res[,2]),
                              ConsensusClustering = as.factor(moic.ConsensusClustering$clust.res[,2]),
                              IntNMF = as.factor(moic.IntNMF$clust.res[,2]),
                              iClusterBayes = as.factor(moic.iClusterBayes$clust.res[,2]))

                              
                             
########################################## Annotated meta new
annCol <- Combine_clust

################################################
annColors   <-  list("SNF" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "CIMLR" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "PINSPlus" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "NEMO" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "COCA" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "MoCluster" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "LRAcluster" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "ConsensusClustering" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "IntNMF" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4"),
                     "iClusterBayes" = c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4")
                    )

getMoHeatmap(data          = plotdata,
             row.title     = c("",""),
             is.binary     = c(F,F), # the 4th data is mutation which is binary
             legend.name   = c("Enrichment","Enrichment"),
             clust.res     = cmoic$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = annCol, # no annotation for samples
             annColors     = annColors, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES single clust")

####################################################################
########################################## Annotated meta old
setwd('H:\\ZF-IMC-Formal-analysis\\Experiment\\S05.MOVICS_cluster_analysis')
Meta_data <- read.csv('./Meta_data_old.csv')
colnames(Meta_data)[1] <- 'Patient_Match'
filter <- read.csv('./SampleID_filter.csv')
filter <- filter$Patient
Meta_data <- Meta_data[Meta_data$Patient_Match %in% filter,]
rownames(Meta_data) <- Meta_data$Patient_Match
Meta_data <- Meta_data[rownames(cmoic$clust.res),]

Meta_data <- Meta_data[,2:4]
annCol <- Meta_data
################################################
annColors   <-  list("Treatment" = c("C+R" = "#F8766D", "C" = "#7CAE00", "C+I" = "#00BFC4", "I" = "#C77CFF"),
                 "Response" = c("R" = "#2DB600FF", "NR" = "#E8C32EFF"),
                 "Outcome" = c( "CR" = "#80C1C0", "PR" = "#D695FE"))

getMoHeatmap(data          = plotdata,
             row.title     = c("Freq","Density","Districts"),
             is.binary     = c(F,F,F), # the 4th data is mutation which is binary
             legend.name   = c("Enrichment","Enrichment","Enrichment"),
             clust.res     = cmoic$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = annCol, # no annotation for samples
             annColors     = annColors, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")










