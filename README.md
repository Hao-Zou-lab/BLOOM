<p align="center">
  <img src="man/figures/logo.png" alt="BLOOM logo" width="200"/>
</p>

# BLOOM

1.BLOOM,R version 4.1.2(This version or any version higher than it is acceptable).

2.BLOOM can transfers plant phenotype-associated signals from bulk RNA-seq to the single-cell level.

3.If you need all.count.txt, mesophyll_cells_all.rds, please download the files in inst/exdata in the BLOOM_0.1.0.tar.gz package to obtain the data you need.



## Step1: installation

1.install 'Seurat','pls', 'MASS','ggplot2','readxl','tcltk','wrMisc' R package from Cran.

2.installing STutility from github: dowaload the file 'STutility-1.1.1tar.gz'/ 'STutility-1.1.1.zip'and install the package from local path.

3.Installing scSTAR from github:

&nbsp;	i. We can install scSCTAR2 by downloading the installation file (BLOOM_0.1.0.tar.gz) from the Assets page of Releases and installing it:

        install.packages("scSTAR2_0.1.1.tar.gz", repos = NULL, type="source")

&nbsp;	ii. Another method is to download the file "BLOOM_0.1.0.tar.gz" and install the package from a local path.



## Step2: run BLOOM

```r

rm(list = ls())

gc()

graphics.off()

library(Seurat)

library(pls)

library(MASS)

library(ggplot2)

library(readxl)

library(tcltk)

library(STutility)

library(wrMisc)

library(scSTAR2)

library(circlize)

library(monocle) 
```

### load scRNA-seq data

```r

setwd("./your path")

mesophyll_cells <- readRDS("./mesophyll_cells_all.rds")

```

### load bulk RNA-seq data with phenotype information

```r

data_bulk = read.table("./all.count.txt", sep = "\\t", header = T)

geneList_bulk = rownames(data_bulk)

data_matrix_bulk <- as.matrix(data_bulk)

data_matrix_bulk <- matrix(as.numeric(data_matrix_bulk), nrow = nrow(data_matrix_bulk))

rownames(data_matrix_bulk) <- geneList_bulk

colnames(data_matrix_bulk) <- colnames(data_bulk)

data_matrix_bulk <- data_matrix_bulk[!is.na(rownames(data_matrix_bulk)), ]

data_matrix_bulk <- log2(data_matrix_bulk + 1)

```

### keep the shared genes between datasets

```r

location <- c(rep('TL0', 3), rep('TL1', 3))

data_matrix_sc <- as.matrix(mesophyll_cells@assays$RNA@counts)

geneList_sc = rownames(data_matrix_sc)

rownames(data_matrix_sc) <- geneList_sc

colnames(data_matrix_sc) <- colnames(mesophyll_cells)

geneList <- intersect(geneList_bulk, geneList_sc)

ia <- match(geneList, geneList_bulk)

ib <- match(geneList, geneList_sc)

data_matrix_bulk <- data_matrix_bulk[ia, ]

data_matrix_sc <- data_matrix_sc[ib, ]

data_matrix_sc <- log2(data_matrix_sc + 1)

geneList_OGFSC_share <- geneList

```

### PLS1 projection

```r

NCV <- 5

minNC <- 2

PLScomp <- 2

MODEL <- PLSconstruct(t(data_matrix_sc), t(data_matrix_bulk), 'mc', NCV, PLScomp, minNC)

```

### reconstruct scRNA-seq data

```r

temp <- t(data_matrix_sc)

temp = temp - matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

temp = temp - temp %*% MODEL$XL %*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

S_sc <- t(temp)

```

### reconstruct bulk

```r
temp = t(data_matrix_bulk)

temp = temp - matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

temp = temp - temp %*% MODEL$XL %*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

S_bulk = t(temp)
```

### meta prediction model training

```r

PLScomp2 = 3 # by default, 3. Might be slightly adjusted to 4 or 5

FCV = 3

data_1 = S_bulk[, location == 'TL0']

data_2 = S_bulk[, location == 'TL1']
```

### train discriminatory model on bulk data

```r

MODEL <- PLSconstruct(t(data_1), t(data_2), 'mc', NCV, PLScomp2, minNC)
```

### apply the model on SC data

```r

temp <- t(S_sc)

temp = temp - matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

temp = temp %*% MODEL$XL %*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %*% MODEL$mu

SS_sc <- t(temp)
```

### save results

```r
geneList_OGFSC_share <- as.list(geneList_OGFSC_share)

location <- as.list(location)

BLOOM_mesophyll <- list(

data_matrix_sc = data_matrix_sc, 

data_matrix_bulk = data_matrix_bulk, 

S_bulk = S_bulk, 

S_sc = S_sc, 

SS_sc = SS_sc, 

geneList_OGFSC_share = geneList_OGFSC_share, 

location = location)

names(BLOOM_mesophyll) <- c('data_matrix_sc', 'data_matrix_bulk', 'S_bulk', 'S_sc', 'SS_sc', 'geneList_OGFSC_share', 'location')

saveRDS(BLOOM_mesophyll, file = './BLOOM_sc.rds')

```

### load data

```r

data <- readRDS("./BLOOM_sc.rds")

SS_sc <- data$SS_sc

SS_sc <- matrix(as.numeric(SS_sc), nrow = nrow(SS_sc), dimnames = dimnames(SS_sc))

genelist <- unlist(data$geneList_OGFSC_share)

rownames(SS_sc) <- genelist

colnames(SS_sc) <- colnames(data_matrix_sc)

data_matrix_sc <- data$data_matrix_sc

data_matrix_sc <- matrix(as.numeric(data_matrix_sc), nrow = nrow(data_matrix_sc), dimnames = dimnames(data_matrix_sc))

rownames(data_matrix_sc) <- genelist

colnames(data_matrix_sc) <- paste0("Cell", 1:ncol(data_matrix_sc))

S_sc <- data$S_sc

S_sc <- matrix(as.numeric(S_sc), nrow = nrow(S_sc), dimnames = dimnames(S_sc))

rownames(S_sc) <- genelist
```

### clustering and identify marker genes of each cluster

```r

Integrated <- CreateSeuratObject(counts = SS_sc, project = "Integrated", min.cells = 0, min.features = 0)

Integrated <- NormalizeData(Integrated, verbose = FALSE)

Integrated <- FindVariableFeatures(Integrated, selection.method = "vst", nfeatures = 2000)

Integrated <- ScaleData(Integrated, verbose = FALSE)

Integrated <- RunPCA(Integrated, npcs = 20, verbose = FALSE)

Integrated <- RunTSNE(Integrated, dims = 1:20)

Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:20)

Integrated <- FindClusters(Integrated, resolution = 0.3)

n <- length(unique(Integrated$seurat_clusters)) - 1

scS_cluster <- paste0('scS_', Integrated$seurat_clusters)

Integrated.markers <- FindAllMarkers(Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)

markers_filtered <- Integrated.markers %>% filter(p_val_adj < 0.05)

write.csv(markers_filtered, "./MarkerGenesByCluster_BLOOM.csv", row.names = FALSE)

library(tidyverse)

library(cowplot)

library(ggplot2)

data1 <- Integrated@meta.data

table(data1$seurat_clusters)

data2 <- Integrated@reductions[["tsne"]]@cell.embeddings%>%as.data.frame()

mydata <- merge(data2,data1[,c(1,5,ncol(data1))],by=0,all.x = T)%>%  

column_to_rownames("Row.names")

mydata$seurat_clusters <- paste0("BM", mydata$seurat_clusters)

sort(unique(mydata$seurat_clusters))

cellcolors <- c("#fb9e9a","#d688a1","#a7789f","#e0845d","#0075c3","#47c4f1","#008eb9","#c06967","#798945")

ggplot(data = mydata, aes(tSNE_1, tSNE_2, fill = seurat_clusters, colour = seurat_clusters)) +

geom_point(shape = 21, size = 1.5, alpha = 0.15) +  scale_fill_manual(values = cellcolors) +

scale_colour_manual(values = cellcolors) +

theme_bw(base_rect_size = 1) +   

theme(

axis.title = element_text(size = 15, family = "sans", face = "bold"),     

axis.text = element_text(size = 12, family = "sans", face = "bold"),      

panel.grid = element_blank(),

legend.title = element_blank(),

legend.text = element_text(size = 12, family = "sans", face = "bold"),    

legend.position = "right",

legend.key.height = unit(1, 'cm'),

legend.key.width = unit(0.5, 'cm'),

plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "sans")    

) +

ggtitle("Mesophyll BM Cluster") +

guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.8)))

ggsave("./Mesophyll_BM_cluster.tiff", height = 5,width = 6)
```

### load bulk data

```r
data_matrix_bulk <- data$data_matrix_bulk

data_matrix_bulk <- matrix(as.numeric(data_matrix_bulk), nrow = nrow(data_matrix_bulk), dimnames = dimnames(data_matrix_bulk))

genelist <- unlist(data$geneList_OGFSC_share)

rownames(data_matrix_bulk) <- genelist

colnames(data_matrix_bulk) <- paste0("Patient", c(1:dim(data_matrix_bulk)[2]))

```

### phenotype data

```r

location <- unlist(data$location)
```

### load gene markers

```r
markers <- read.csv("./MarkerGenesByCluster\_BLOOM.csv")

```

### Processing of gene marker data may involve extracting related genes by cluster

This part of the code involves selecting specific genes from the clustering information for subsequent analysis.

```r

clusters <- markers$cluster

clusters_uni <- sort(unique(clusters))

FC <- markers$avg_log2FC

selectedGenes <- list()

for (i in clusters_uni) {

idx <- which(clusters == i)

x <- sort(markers$avg_log2FC[idx], decreasing = TRUE, index.return = TRUE)

selectedGenes[[paste0("cluster", i)]] <- as.vector(markers$gene[idx[x$ix]])

}
```

This part evaluates genetic similarity between each cluster and bulk data with different meta information.

```r
buffer <- data$S_bulk

buffer <- matrix(as.numeric(buffer), nrow = nrow(buffer), dimnames = dimnames(buffer))

genelist <- unlist(data$geneList_OGFSC_share)

rownames(buffer) <- genelist

ngenes <- 200

```

For each cluster, extract the specified number of genes from the bulk data and the single-cell data of the corresponding cluster, then calculate the Pearson correlation between these genes, and save the results in the list corrMat

```r

corrMat <- list()

for (i in 0:n) {

cluster_name <- paste0("scS_", i)

X <- SS_sc[, which(scS_cluster == cluster_name)]

selected_genes <- selectedGenes[[paste0("cluster", i)]][1:ngenes]

common_genes <- intersect(rownames(buffer), selected_genes)

if (length(common_genes) < 10) {

cat("Skip cluster", i, "due to too few genes.\\n")

next

}

Cor <- cor(buffer[common_genes, ], X[common_genes, ], method = "pearson")

corrMat[[paste0("Cor", i)]] <- Cor

}

corrMat[["location"]] <- as.list(location)

saveRDS(corrMat, file = './corrMat.rds')
```

### This code mainly implements the OPLS-DA (Partial Least Squares Discriminant Analysis) process

```r

data <- readRDS("./corrMat.rds")

pattern <- do.call(rbind, data$location)

idx_AC <- which(pattern == "TL0", arr.ind = TRUE)

idx_SA <- which(pattern == "TL1", arr.ind = TRUE)

cor_names <- names(data)[grepl("^Cor\\\\d+$", names(data))]

cor_indices <- as.integer(gsub("Cor", "", cor_names))

for (i in cor_indices) {

gc()

CorMat <- data[[paste0("Cor", i)]]

data1 <- CorMat[idx_AC[, 1], ]

data2 <- CorMat[idx_SA[, 1], ]

cellindex <- matrix(seq(1, ncol(CorMat)), nrow = 1)

Sys.sleep(1) 

set.seed(123)

model <- OPLSDA(data1, data2,cellindex,nrcv=3)

ggsave(file = paste0('C', i, '.tiff'), dpi = 300, compression = 'lzw', width = 6, height = 5, units = "in")

dev.off()

}
```




































&nbsp;	

