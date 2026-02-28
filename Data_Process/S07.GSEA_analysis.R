#####################
library(Seurat)
library(dplyr)
library(stringr)
library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GseaVis)
library(dplyr)
library(stringr)
library(patchwork)
library(knitr)
library(msigdbr) 
library(fgsea)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(limma)
library(clusterProfiler)
library(DESeq2)
library(ggplotify)
library(cowplot)
library(edgeR)
library(tinyarray)
library(gridExtra)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(enrichplot)
library(data.table)
###########################
setwd('H:\\Reng\\DEG_analysis')
###########################
h.gmt <- msigdbr(species = 'human', category = "H") %>%
            dplyr::select(gs_name, gene_symbol)
#################################### Cluster 1
Cluster3_DEG <- fread('./NPC_deseq2_test_result.Cluster3_vs_Others.txt')

DEG <- Cluster3_DEG
colnames(DEG)[1] <- 'Gene'
colnames(DEG)[3] <- 'logFC'
colnames(DEG)[7] <- 'adj.P.Val'
##################################### GSEA
rownames(DEG) <- DEG$Gene
IR <- DEG[!duplicated(DEG$Gene),c('Gene','logFC','adj.P.Val')]
#########################
library(org.Hs.eg.db)
gene <- bitr(IR$Gene,     #转换的列是nrDEG的列名
             fromType = "SYMBOL",     #需要转换ID类型
             toType =  "ENTREZID",    #转换成的ID类型
             OrgDb = org.Hs.eg.db)    #对应的物种，小鼠的是org.Mm.eg.db
#让基因名、ENTREZID、logFC对应起来
gene$logFC <- IR$logFC[match(gene$SYMBOL,IR$Gene)]
head(gene)
geneList=gene$logFC
names(geneList)=gene$SYMBOL 
geneList <- geneList[!duplicated(names(geneList))]
#按照logFC的值来排序geneList
geneList=sort(geneList,decreasing = T)
head(geneList)

Plot_GeneSet <- h.gmt
colnames(Plot_GeneSet) <- c('term','gene')

set.seed(1234)
egmt <- GSEA(geneList, TERM2GENE = Plot_GeneSet, 
             minGSSize = 10,maxGSSize = 10000,
             pvalueCutoff = 1,pAdjustMethod="fdr", 
                    seed=TRUE, by="fgsea")

gesa_res <- as.data.frame(egmt@result)
###############################################
#gseaplot2(egmt, "HALLMARK_INFLAMMATORY_RESPONSE", color = "firebrick",
#            rel_heights=c(1, .2, .6), subplots=1:2)

############################################### Inflammatory Res
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INFLAMMATORY_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('H:\\Reng\\DEG_analysis')
pdf(file="./Cluster3_InflammatoryRes.pdf",width= 10,height= 4)
p1
dev.off()

############################################### IL6_JAK_STAT3_SIGNALING
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_IL6_JAK_STAT3_SIGNALING",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('H:\\Reng\\DEG_analysis')
pdf(file="./Cluster3_IL6_JAK.pdf",width= 10,height= 4)
p1
dev.off()

############################################### INTERFERON_ALPHA_RESPONSE
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('H:\\Reng\\DEG_analysis')
pdf(file="./Cluster3_INTERFERON_ALPHA.pdf",width= 10,height= 4)
dev.off()

############################################### HALLMARK_INTERFERON_GAMMA_RESPONSE
p1 <- gseaplot2(egmt,#数据
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",#画哪一列的信号通路
            #title = "Prion disease",#标题
            base_size = 15,#字体大小
            color = "green",#线条的颜色
            pvalue_table = TRUE,#加不加p值
            ES_geom="line",subplots=1:2)#是用线，还是用d点

setwd('H:\\Reng\\DEG_analysis')
pdf(file="./Cluster3_INTERFERON_GAMMA.pdf",width= 10,height= 4)
p1
dev.off()



