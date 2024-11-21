##############################################################################
# analogous to how David assigns indices to his mutants
# egl-9, 2, b -- daf-18
# hif-1, 4, c -- daf-16
# egl-9 hif-1, 7, f -- double
##############################################################################

# load libraries
library(tidyverse)
library(magrittr)
library(DESeq2)
library(readxl)

# load data
ws273_GeneInfo <- read_csv("ws273_gene_info.csv", col_names = TRUE) %>% as.data.frame()
colnames(ws273_GeneInfo)[1] <- "gene_id"

bg <- read_tsv("background_gene_ids.txt", col_names = FALSE)
colnames(bg) <- "gene_id"

tmp1 <- read_csv("20200902_RNAseq_raw_reads_merged.csv", col_names = TRUE)
tmp1 %<>% filter(gene_id %in% bg$gene_id)

# creat countData tables
dbl_vs_n2 <- 
  as.data.frame(tmp1[, 1]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "N2_starv")]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "double")])
colnames(dbl_vs_n2)
dbl_vs_n2 <- dbl_vs_n2[, c(1, 3,5,4, 9,6,7)]

d18_vs_n2 <- 
  as.data.frame(tmp1[, 1]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "N2_starv")]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "daf-18")])
colnames(d18_vs_n2)
d18_vs_n2 <- d18_vs_n2[, c(1, 3,5,4, 8,9,6)]

d16_vs_n2 <- 
  as.data.frame(tmp1[, 1]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "N2_starv")]) %>%
  cbind(tmp1[, str_detect(colnames(tmp1), pattern = "daf-16")])
colnames(d16_vs_n2)
d16_vs_n2 <- d16_vs_n2[, c(1, 3,5,4, 8,9,7)]

# create metadata tables
dbl_meta <- data.frame(id = colnames(dbl_vs_n2)[2:7],
                       dex = c(rep("control", 3),
                               rep("treated", 3)),
                       condition = c(rep("N2_starv", 3),
                                     rep("double", 3)))

d18_meta <- data.frame(id = colnames(d18_vs_n2)[2:7],
                       dex = c(rep("control", 3),
                               rep("treated", 3)),
                       condition = c(rep("N2_starv", 3),
                                     rep("daf-18", 3)))

d16_meta <- data.frame(id = colnames(d16_vs_n2)[2:7],
                       dex = c(rep("control", 3),
                               rep("treated", 3)),
                       condition = c(rep("N2_starv", 3),
                                     rep("daf-16", 3)))

# DE analysis for double vs N2_starv
countData <- dbl_vs_n2
rownames(countData) <- countData$gene_id
metaData <- dbl_meta
metaData$dex %<>% as.factor()
levels(metaData$dex)

# tidy = TRUE: the first column of countData is row names of the count matrix
# count matrix: the input of the function "DESeqDataSetFromMatrix"
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~ dex,
                              tidy = TRUE)
dds <- DESeq(dds)

# log2FoldChange -- comparison: the LAST level of the last variable in the design formula OVER the FIRST level
# levels(metaData$dex)
# tidy = TRUE outputs the results table with row names as the first column
res <- results(dds, tidy = TRUE)
colnames(res)[1] <- "gene_id"
# head(res)
# summary(res)
res <- res[order(res$padj), ]
head(res)

# blind = TRUE is used to perform sample quality assurance, etc.
# blind = FALSE is used to transform data for downstream analysis and makes use of the design matrix
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "dex")

# convert to the format that David's analysis takes in
DE_genes <- merge(ws273_GeneInfo, res, by = "gene_id", all = FALSE, sort = FALSE)
DE_genes <- DE_genes[, c(1:3, 15,16,19)]
colnames(DE_genes) <- c("ens_gene", "ext_gene", "target_id",
                        "b", "se_b", "qval")
DE_genes$genotype <- "daf-16;daf-18"
DE_genes$sorter <- "7"
DE_genes$code <- "f"
DE_genes_dbl <- DE_genes

# DE analysis for daf-18 vs N2_starv
countData <- d18_vs_n2
rownames(countData) <- countData$gene_id
metaData <- d18_meta
metaData$dex %<>% as.factor()
levels(metaData$dex)

# tidy = TRUE: the first column of countData is row names of the count matrix
# count matrix: the input of the function "DESeqDataSetFromMatrix"
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~ dex,
                              tidy = TRUE)
dds <- DESeq(dds)

# log2FoldChange -- comparison: the LAST level of the last variable in the design formula OVER the FIRST level
# levels(metaData$dex)
# tidy = TRUE outputs the results table with row names as the first column
res <- results(dds, tidy = TRUE)
colnames(res)[1] <- "gene_id"
# head(res)
# summary(res)
res <- res[order(res$padj), ]
head(res)

# blind = TRUE is used to perform sample quality assurance, etc.
# blind = FALSE is used to transform data for downstream analysis and makes use of the design matrix
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "dex")

# convert to the format that David's analysis takes in
DE_genes <- merge(ws273_GeneInfo, res, by = "gene_id", all = FALSE, sort = FALSE)
DE_genes <- DE_genes[, c(1:3, 15,16,19)]
colnames(DE_genes) <- c("ens_gene", "ext_gene", "target_id",
                        "b", "se_b", "qval")
DE_genes$genotype <- "daf-18"
DE_genes$sorter <- "2"
DE_genes$code <- "b"
DE_genes_d18 <- DE_genes

# DE analysis for daf-16 vs N2_starv
countData <- d16_vs_n2
rownames(countData) <- countData$gene_id
metaData <- d16_meta
metaData$dex %<>% as.factor()
levels(metaData$dex)

# tidy = TRUE: the first column of countData is row names of the count matrix
# count matrix: the input of the function "DESeqDataSetFromMatrix"
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~ dex,
                              tidy = TRUE)
dds <- DESeq(dds)

# log2FoldChange -- comparison: the LAST level of the last variable in the design formula OVER the FIRST level
# levels(metaData$dex)
# tidy = TRUE outputs the results table with row names as the first column
res <- results(dds, tidy = TRUE)
colnames(res)[1] <- "gene_id"
# head(res)
# summary(res)
res <- res[order(res$padj), ]
head(res)

# blind = TRUE is used to perform sample quality assurance, etc.
# blind = FALSE is used to transform data for downstream analysis and makes use of the design matrix
vsdata <- vst(dds, blind = FALSE)
plotPCA(vsdata, intgroup = "dex")

# convert to the format that David's analysis takes in
DE_genes <- merge(ws273_GeneInfo, res, by = "gene_id", all = FALSE, sort = FALSE)
DE_genes <- DE_genes[, c(1:3, 15,16,19)]
colnames(DE_genes) <- c("ens_gene", "ext_gene", "target_id",
                        "b", "se_b", "qval")
DE_genes$genotype <- "daf-16"
DE_genes$sorter <- "4"
DE_genes$code <- "c"
DE_genes_d16 <- DE_genes

# save into a CSV that will go into David's analysis
DE_genes <- rbind(DE_genes_d16, DE_genes_dbl, DE_genes_d18)
write_csv(DE_genes, "DE_genes_jingxian.csv", col_names = TRUE)
