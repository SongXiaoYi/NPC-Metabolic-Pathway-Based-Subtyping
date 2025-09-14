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
#接下来运行S02


