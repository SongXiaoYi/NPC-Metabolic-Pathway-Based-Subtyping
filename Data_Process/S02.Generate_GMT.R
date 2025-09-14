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
######################################
setwd('H:\\Reng')
######################################Building GMT
gset <- fread('./GeneSet.csv') %>% as.data.frame()
gset <- split(gset, gset$`Metabolic pathways involved`)
for (n in 1:length(gset)) {
  gset[[n]] <- as.list(unique(gset[[n]]$`Gene symbol`))
}

file=("./gset_total.gmt")

######################################Output GMT function
output_gmt<-function(geneset,file){
sink(file) #将输出结果重定向到file
lapply(names(geneset),function(i){
#按照列名，将要连接的变量转化为向量型，用制表符连接成一个字符串
cat(paste(c(i,'tmp',geneset[[i]]),collapse='\t')) 
cat('\n') #输出后新起一行
})
sink() #结束重定向
}

######################################## Conduct
output_gmt(gset,file)
#######################################################
gmt <- getGmt('gset_total.gmt')
pbmc <- pbmc[,-1]
TCGA <- as.matrix(pbmc)
TCGA <- avereps(TCGA)
gs.exp <- gsva(TCGA,gmt,method = "ssgsea",kcdf = "Poisson", min.sz = 1)

#DR_score <- gs.exp[190,] - gs.exp[191,]
setwd('H:\\Reng')
gs.exp <- as.data.frame(gs.exp)
write.csv(gs.exp, file = 'Gset_analysis_score.csv', row.names = FALSE)


