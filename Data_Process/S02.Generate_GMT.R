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
  gset[[n]] <- as.list(gset[[n]]$`Gene symbol`)
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




