# \# BLOOM <img src="man/figures/logo.png" align="right" height="138" />



\##BLOOM

1.BLOOM,R version 4.1.2(This version or any version higher than it is acceptable).

2.BLOOM can transfers plant phenotype-associated signals from bulk RNA-seq to the single-cell level.

3.If you need all.count.txt, mesophyll\_cells\_all.rds, please download the files in inst/exdata in the BLOOM\_0.1.0.tar.gz package to obtain the data you need.



\##Step1: installation

1.install 'Seurat','pls', 'MASS','ggplot2','readxl','tcltk','wrMisc' R package from Cran.

2.installing STutility from github: dowaload the file 'STutility-1.1.1tar.gz'/ 'STutility-1.1.1.zip'and install the package from local path.

3.Installing scSTAR from github:

&nbsp;	i. We can install scSCTAR2 by downloading the installation file (BLOOM\_0.1.0.tar.gz) from the Assets page of Releases and installing it:

&nbsp;	```r

&nbsp;	install.packages("scSTAR2\_0.1.1.tar.gz", repos = NULL, type="source")

&nbsp;	```

&nbsp;	ii. Another method is to download the file "BLOOM\_0.1.0.tar.gz" and install the package from a local path.



\##Step2: run BLOOM

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

\###load scRNA-seq data

```r

setwd("./your path")

mesophyll\_cells <- readRDS("./mesophyll\_cells\_all.rds")

```

\### load bulk RNA-seq data with phenotype information

```r

data\_bulk = read.table("./all.count.txt", sep = "\\t", header = T)

geneList\_bulk = rownames(data\_bulk)

data\_matrix\_bulk <- as.matrix(data\_bulk)

data\_matrix\_bulk <- matrix(as.numeric(data\_matrix\_bulk), nrow = nrow(data\_matrix\_bulk))

rownames(data\_matrix\_bulk) <- geneList\_bulk

colnames(data\_matrix\_bulk) <- colnames(data\_bulk)

data\_matrix\_bulk <- data\_matrix\_bulk\[!is.na(rownames(data\_matrix\_bulk)), ]

data\_matrix\_bulk <- log2(data\_matrix\_bulk + 1)

```

\###keep the shared genes between datasets

```r

location <- c(rep('TL0', 3), rep('TL1', 3))

data\_matrix\_sc <- as.matrix(mesophyll\_cells@assays$RNA@counts)

geneList\_sc = rownames(data\_matrix\_sc)

rownames(data\_matrix\_sc) <- geneList\_sc

colnames(data\_matrix\_sc) <- colnames(mesophyll\_cells)

geneList <- intersect(geneList\_bulk, geneList\_sc)

ia <- match(geneList, geneList\_bulk)

ib <- match(geneList, geneList\_sc)

data\_matrix\_bulk <- data\_matrix\_bulk\[ia, ]

data\_matrix\_sc <- data\_matrix\_sc\[ib, ]

data\_matrix\_sc <- log2(data\_matrix\_sc + 1)

geneList\_OGFSC\_share <- geneList

```

\###PLS1 projection

```r

NCV <- 5

minNC <- 2

PLScomp <- 2

MODEL <- PLSconstruct(t(data\_matrix\_sc), t(data\_matrix\_bulk), 'mc', NCV, PLScomp, minNC)

```

\###reconstruct scRNA-seq data

```r

temp <- t(data\_matrix\_sc)

temp = temp - matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

temp = temp - temp %\*% MODEL$XL %\*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

S\_sc <- t(temp)

```

\###reconstruct bulk

```r
temp = t(data\_matrix\_bulk)

temp = temp - matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

temp = temp - temp %\*% MODEL$XL %\*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

S\_bulk = t(temp)
```

\###meta prediction model training

```r

PLScomp2 = 3 # by default, 3. Might be slightly adjusted to 4 or 5

FCV = 3

data\_1 = S\_bulk\[, location == 'TL0']

data\_2 = S\_bulk\[, location == 'TL1']
```

\###train discriminatory model on bulk data

```r

MODEL <- PLSconstruct(t(data\_1), t(data\_2), 'mc', NCV, PLScomp2, minNC)
```

\###apply the model on SC data

```r

temp <- t(S\_sc)

temp = temp - matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

temp = temp %\*% MODEL$XL %\*% ginv(MODEL$XL)

temp = temp + matrix(1, dim(temp)\[1], 1) %\*% MODEL$mu

SS\_sc <- t(temp)
```

\###save results

```r
geneList\_OGFSC\_share <- as.list(geneList\_OGFSC\_share)

location <- as.list(location)

BLOOM\_mesophyll <- list(

&nbsp; data\_matrix\_sc = data\_matrix\_sc, 

&nbsp; data\_matrix\_bulk = data\_matrix\_bulk, 

&nbsp; S\_bulk = S\_bulk, 

&nbsp; S\_sc = S\_sc, 

&nbsp; SS\_sc = SS\_sc, 

&nbsp; geneList\_OGFSC\_share = geneList\_OGFSC\_share, 

&nbsp; location = location)

names(BLOOM\_mesophyll) <- c('data\_matrix\_sc', 'data\_matrix\_bulk', 'S\_bulk', 'S\_sc', 'SS\_sc', 'geneList\_OGFSC\_share', 'location')

saveRDS(BLOOM\_mesophyll, file = './BLOOM\_sc.rds')

```

\###load data

```r

data <- readRDS("./BLOOM\_sc.rds")

SS\_sc <- data$SS\_sc

SS\_sc <- matrix(as.numeric(SS\_sc), nrow = nrow(SS\_sc), dimnames = dimnames(SS\_sc))

genelist <- unlist(data$geneList\_OGFSC\_share)

rownames(SS\_sc) <- genelist

colnames(SS\_sc) <- colnames(data\_matrix\_sc)

data\_matrix\_sc <- data$data\_matrix\_sc

data\_matrix\_sc <- matrix(as.numeric(data\_matrix\_sc), nrow = nrow(data\_matrix\_sc), dimnames = dimnames(data\_matrix\_sc))

rownames(data\_matrix\_sc) <- genelist

colnames(data\_matrix\_sc) <- paste0("Cell", 1:ncol(data\_matrix\_sc))

S\_sc <- data$S\_sc

S\_sc <- matrix(as.numeric(S\_sc), nrow = nrow(S\_sc), dimnames = dimnames(S\_sc))

rownames(S\_sc) <- genelist
```

\###clustering and identify marker genes of each cluster

```r

Integrated <- CreateSeuratObject(counts = SS\_sc, project = "Integrated", min.cells = 0, min.features = 0)

Integrated <- NormalizeData(Integrated, verbose = FALSE)

Integrated <- FindVariableFeatures(Integrated, selection.method = "vst", nfeatures = 2000)

Integrated <- ScaleData(Integrated, verbose = FALSE)

Integrated <- RunPCA(Integrated, npcs = 20, verbose = FALSE)

Integrated <- RunTSNE(Integrated, dims = 1:20)

Integrated <- FindNeighbors(Integrated, reduction = "pca", dims = 1:20)

Integrated <- FindClusters(Integrated, resolution = 0.3)

n <- length(unique(Integrated$seurat\_clusters)) - 1

scS\_cluster <- paste0('scS\_', Integrated$seurat\_clusters)

Integrated.markers <- FindAllMarkers(Integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)

markers\_filtered <- Integrated.markers %>% filter(p\_val\_adj < 0.05)

write.csv(markers\_filtered, "./MarkerGenesByCluster\_BLOOM.csv", row.names = FALSE)

library(tidyverse)

library(cowplot)

library(ggplot2)

data1 <- Integrated@meta.data

table(data1$seurat\_clusters)

data2 <- Integrated@reductions\[\["tsne"]]@cell.embeddings%>%as.data.frame()

mydata <- merge(data2,data1\[,c(1,5,ncol(data1))],by=0,all.x = T)%>%  

&nbsp; column\_to\_rownames("Row.names")

mydata$seurat\_clusters <- paste0("BM", mydata$seurat\_clusters)

sort(unique(mydata$seurat\_clusters))

cellcolors <- c("#fb9e9a","#d688a1","#a7789f","#e0845d","#0075c3","#47c4f1","#008eb9","#c06967","#798945")

ggplot(data = mydata, aes(tSNE\_1, tSNE\_2, fill = seurat\_clusters, colour = seurat\_clusters)) +

&nbsp; geom\_point(shape = 21, size = 1.5, alpha = 0.15) +   scale\_fill\_manual(values = cellcolors) +

&nbsp; scale\_colour\_manual(values = cellcolors) +

&nbsp; theme\_bw(base\_rect\_size = 1) +   

&nbsp; theme(

&nbsp;   axis.title = element\_text(size = 15, family = "sans", face = "bold"),     

&nbsp;   axis.text = element\_text(size = 12, family = "sans", face = "bold"),      

&nbsp;   panel.grid = element\_blank(),

&nbsp;   legend.title = element\_blank(),

&nbsp;   legend.text = element\_text(size = 12, family = "sans", face = "bold"),    

&nbsp;   legend.position = "right",

&nbsp;   legend.key.height = unit(1, 'cm'),

&nbsp;   legend.key.width = unit(0.5, 'cm'),

&nbsp;   plot.title = element\_text(size = 16, face = "bold", hjust = 0.5, family = "sans")    

&nbsp; ) +

&nbsp; ggtitle("Mesophyll BM Cluster") +

&nbsp; guides(fill = guide\_legend(override.aes = list(size = 3, alpha = 0.8)))

ggsave("./Mesophyll\_BM\_cluster.tiff", height = 5,width = 6)
```

\### load bulk data

```r
data\_matrix\_bulk <- data$data\_matrix\_bulk

data\_matrix\_bulk <- matrix(as.numeric(data\_matrix\_bulk), nrow = nrow(data\_matrix\_bulk), dimnames = dimnames(data\_matrix\_bulk))

genelist <- unlist(data$geneList\_OGFSC\_share)

rownames(data\_matrix\_bulk) <- genelist

colnames(data\_matrix\_bulk) <- paste0("Patient", c(1:dim(data\_matrix\_bulk)\[2]))

```

\###phenotype data

```r

location <- unlist(data$location)
```

\###load gene markers

```r
markers <- read.csv("./MarkerGenesByCluster\_BLOOM.csv")

```

\###Processing of gene marker data may involve extracting related genes by cluster

This part of the code involves selecting specific genes from the clustering information for subsequent analysis.

```r

clusters <- markers$cluster

clusters\_uni <- sort(unique(clusters))

FC <- markers$avg\_log2FC

selectedGenes <- list()

for (i in clusters\_uni) {

&nbsp; idx <- which(clusters == i)

&nbsp; x <- sort(markers$avg\_log2FC\[idx], decreasing = TRUE, index.return = TRUE)

&nbsp; selectedGenes\[\[paste0("cluster", i)]] <- as.vector(markers$gene\[idx\[x$ix]])

}
```

This part evaluates genetic similarity between each cluster and bulk data with different meta information.

```r
buffer <- data$S\_bulk

buffer <- matrix(as.numeric(buffer), nrow = nrow(buffer), dimnames = dimnames(buffer))

genelist <- unlist(data$geneList\_OGFSC\_share)

rownames(buffer) <- genelist

ngenes <- 200

```

For each cluster, extract the specified number of genes from the bulk data and the single-cell data of the corresponding cluster, then calculate the Pearson correlation between these genes, and save the results in the list corrMat

```r

corrMat <- list()

for (i in 0:n) {

&nbsp; cluster\_name <- paste0("scS\_", i)

&nbsp; X <- SS\_sc\[, which(scS\_cluster == cluster\_name)]

&nbsp; selected\_genes <- selectedGenes\[\[paste0("cluster", i)]]\[1:ngenes]

&nbsp; common\_genes <- intersect(rownames(buffer), selected\_genes)

&nbsp; if (length(common\_genes) < 10) {

&nbsp;   cat("Skip cluster", i, "due to too few genes.\\n")

&nbsp;   next

&nbsp; }

&nbsp; Cor <- cor(buffer\[common\_genes, ], X\[common\_genes, ], method = "pearson")

&nbsp; corrMat\_Epidermal\[\[paste0("Cor", i)]] <- Cor

}

corrMat\[\["location"]] <- as.list(location)

saveRDS(corrMat, file = './corrMat.rds')
```

\###This code mainly implements the OPLS-DA (Partial Least Squares Discriminant Analysis) process

```r

data <- readRDS("./corrMat.rds")

pattern <- do.call(rbind, data$location)

idx\_AC <- which(pattern == "TL0", arr.ind = TRUE)

idx\_SA <- which(pattern == "TL1", arr.ind = TRUE)

cor\_names <- names(data)\[grepl("^Cor\\\\d+$", names(data))]

cor\_indices <- as.integer(gsub("Cor", "", cor\_names))

for (i in cor\_indices) {

&nbsp; gc()

&nbsp; CorMat <- data\[\[paste0("Cor", i)]]

&nbsp; data1 <- CorMat\[idx\_AC\[, 1], ]

&nbsp; data2 <- CorMat\[idx\_SA\[, 1], ]

&nbsp; cellindex <- matrix(seq(1, ncol(CorMat)), nrow = 1)

&nbsp; Sys.sleep(1) 

&nbsp; set.seed(123)

&nbsp; model <- OPLSDA(data1, data2,cellindex,nrcv=3)

&nbsp; ggsave(file = paste0('C', i, '.tiff'), dpi = 300, compression = 'lzw', width = 6, height = 5, units = "in")

&nbsp; dev.off()

}
```




































&nbsp;	

