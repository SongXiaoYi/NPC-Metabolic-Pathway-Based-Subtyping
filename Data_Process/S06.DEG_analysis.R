library(biomaRt)
library(SummarizedExperiment)
library(TCGAbiolinks)
####### Dietary Restriction GSVA Score
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(data.table)
##################################################

check_package <- function(package){
    if (!requireNamespace(package, quietly = TRUE)) {
        stop(package, " package is needed for this function to work. Please install it.",
             call. = FALSE)
    }
}

get.GRCh.bioMart <- function(
    genome = c("hg19", "hg38"),
    as.granges = FALSE
) {

    genome <- match.arg(genome)
    # Since the amount of users complaining about the access we
    # also added the data into a data package
    check_package("TCGAbiolinksGUI.data")
    if (genome == "hg19") {
        gene.location <- get(
            data("gene.location.hg19",
                 package = "TCGAbiolinksGUI.data",
                 envir = environment())
        )

        if (as.granges) {
            gene.location$strand[gene.location$strand == 1] <- "+"
            gene.location$strand[gene.location$strand == -1] <- "-"
            gene.location$chromosome_name <- paste0("chr",gene.location$chromosome_name)
            gene.location <- makeGRangesFromDataFrame(
                gene.location, seqnames.field = "chromosome_name",
                start.field = "start_position",
                end.field = "end_position",
                keep.extra.columns = TRUE
            )
        }
    } else {
        gene.location <- get(
            data(
                "gencode.v36.annotation.genes",
                package = "TCGAbiolinksGUI.data",
                envir = environment()
            )
        )
        if(!as.granges) gene.location <- as.data.frame(gene.location)
    }

    return(gene.location)
}

##########################################
##########################################
#Perform

#############################################
setwd('H:\\Reng')
pbmc <- fread('./Gene_counts.csv') %>% as.data.frame()
rownames(pbmc) <- pbmc$id
genome = "hg38"
gene.location <- get.GRCh.bioMart(genome)
gene.location <- gene.location[, c("gene_name","gene_type")]
gene.location <- gene.location[gene.location$gene_type == 'protein_coding',]
pbmc <- pbmc[rownames(pbmc) %in% gene.location$gene_name,]
######################################

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
#############################]
################# Cluster
Cluster_meta <- read.csv('Cluster_num.csv')
Cluster_meta$clust <- paste0('Cluster', Cluster_meta$clust)
##############################
expr <- pbmc[,-1]
##################################################################
# 亚型名称
n.sub.label <- unique(Cluster_meta$clust) 
# 亚型个数
n.sub <- length(table(Cluster_meta$clust)) 
#创建配对比较的列表信息 
group <- Cluster_meta$clust
names(group) <- Cluster_meta$samID
# 创建需要配对比较的列表函数,创建了三个分组。
createList <- function(group=NULL) {
  tumorsam <- names(group)
  sampleList = list()
  treatsamList =list()
  treatnameList <- c()
  ctrlnameList <- c()
  #A-1: 类1 vs 其他
  sampleList[[1]] = tumorsam
  treatsamList[[1]] = intersect(tumorsam, names(group[group=="Cluster1"])) # 亚型名称需要根据情况修改
  treatnameList[1] <- "Cluster1" # 该亚型的命名
  ctrlnameList[1] <- "Others" # 其他亚型的命名
  #A-2: 类2 vs 其他
  sampleList[[2]] = tumorsam
  treatsamList[[2]] = intersect(tumorsam, names(group[group=="Cluster2"]))
  treatnameList[2] <- "Cluster2"
  ctrlnameList[2] <- "Others"
  #A-3: 类3 vs 其他
  sampleList[[3]] = tumorsam
  treatsamList[[3]] = intersect(tumorsam, names(group[group=="Cluster3"]))
  treatnameList[3] <- "Cluster3"
  ctrlnameList[3] <- "Others"
  #如果有更多类，按以上规律继续写
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}
complist <- createList(group=group)
# 配对DESeq2函数
twoclassDESeq2 <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  options(warn=1)
  for (k in 1:length(sampleList)) { # 循环读取每一次比较的内容
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    compname <- paste(treatname, "_vs_", ctrlname, sep="") # 生成最终文件名
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { # 因为差异表达分析较慢，因此如果文件存在，在不覆盖的情况下（overwt=F）不再次计算差异表达
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    # 差异表达过程，具体参数细节及输出结果解释，请参阅相关document
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula("~ Type")) # 设计矩阵仅包含亚型信息
    dds$Type <- relevel(dds$Type,ref = "control")
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    #将分析结果转化为数据框
resData <- as.data.frame(res[order(res$padj),])
    #将行名作为id列
resData$id <- rownames(resData)
    #提取想要的列数据
resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    #修改列名
colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")
    #输出到文件
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}
# 差异表达分析过程比较慢请耐心等待
twoclassDESeq2(res.path = ".", #所有配对差异表达结果都会输出在res.path路径下
               countsTable = expr[,intersect(colnames(expr),Cluster_meta$samID)],
               prefix = "NPC", #文件名以SKCM开头
               complist = complist,
               overwt = F)

