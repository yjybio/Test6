SRSTCK: An Interpretable Algorithm for the Single-Cell Recognition of Senescent Tumor Cells Based on Kullback-Leibler Divergence
===
Senescence has been listed as one of the cancer hallmarks, and numerous anti-tumor treatments could cause tumor cell senescence, including radiotherapy, chemotherapy, and targeted therapies. Meanwhile, senescent tumor cells (STCs) have been found in diverse cancer stages. STCs could induce immunosuppression of cancer, and then promote the development of tumor. It is therefore critical to understand how STCs contribute to cancer progression by identifying STCs in cancer. To solve this problem, we proposed SRSTCK, a supervised STCs identifier that could accurately identify STCs in cancer single-cell RNA-seq data based on Kullback-Leibler Divergence. In this vignette, we will predict senescent tumor cells by calculating 10X single-cell RNA data using the {SRSTCK} R package.<br>

![image text](https://github.com/yjybio/Test5/blob/main/workflows/workflows.png)

Step 1: installation
--
Installing SRSTCM from GitHub <br>
```R
library(devtools) 
install_github("yjybio/SRSTCK")
```
We start by loading the packages needed for the analyses. Please install them if you haven't.<br>
```R
library("copykat")
library("Seurat")
library("SeuratObject")
```
Step 2: identify malignant cells in single-cell RNA-seq transcriptome data
--
```R
scRNA_tumor_cell <- as.matrix(data@assays$RNA@counts)
scRNA_tumor_cell <- apply(scRNA_tumor_cell, 2, as.numeric)
rownames(scRNA_tumor_cell) <- rownames(data@assays$RNA@counts)
colnames(scRNA_tumor_cell) <- colnames(data@assays$RNA@counts)
TME_cells <-  colnames(data)[data$Cell_Type %in% cell_type]
copykat.test <- copykat(rawmat = scRNA_tumor_cell, sam.name = cancer, norm.cell.names = TME_cells)
sample.name <- paste(cancer, "_copykat_", sep = "")
res <- paste(sample.name, "prediction.txt", sep = "")
copykat.prediction <- read.table(res)
copykat.prediction <- copykat.prediction[which(copykat.prediction[, 2] == "aneuploid"), ]
copykat.prediction <- data@assays$RNA@counts[ ,which(colnames(data@assays$RNA@counts) %in% copykat.prediction[, 1])]
```
Step 3: running SRSTCK
--
```R
copykat.prediction <- copykat.prediction[which(rownames(copykat.prediction) %in% sen.gene), ]
dir <- which(sen.gene %in% rownames(copykat.prediction))
if (length(dir) == 0 || length(dir) == 1)
    stop("too few senescence feature genes;cannot be calculated")
if (length(dir) > 1) {
    copykat.prediction <- copykat.prediction[match(sen.gene[dir], rownames(copykat.prediction)), ]
    p.sen.model <- sen.model[dir]
	p.sen.model <- as.numeric(p.sen.model)
	p.non.sen.model <- non.sen.model[dir]
	p.non.sen.model <- as.numeric(p.non.sen.model)
    q_Cancer <- apply(copykat.prediction, 2, function(x) (x + 0.1) / sum(x + 0.1))
    argmin <- function(p,q){
        result <- 0
        for(i in 1:length(p)){
            result <- result + (q[i] * log(p[i]) - q[i] * log(q[i]))
        }
        return(-result)
       }
    sen_Cancer <- matrix(NA, ncol(copykat.prediction), 2)
	for(i in 1:ncol(copykat.prediction)){
        if(argmin(p.sen.model, q_Cancer[, i]) < argmin(p.non.sen.model, q_Cancer[, i])){
            sen_Cancer[i, 2] <- "senescent"
          } else{
            sen_Cancer[i, 2] <- "non_senescent"
          }
		    sen_Cancer[i, 1] <- colnames(copykat.prediction)[i]
	}
	colnames(sen_Cancer) <- c("Cell_ID", "State")
	sen_Cancer <- as.data.frame(sen_Cancer)
	write.table(sen_Cancer, paste(cancer, "_STC_prediction.txt", sep = ""), sep = "\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
    }
```
Step 4: if the data type is bulk, the senescence score is calculated
--
```R
name <- rownames(data)
data <- apply(data, 2, as.numeric)
rownames(data) <- name
data <- data[which(rownames(data) %in% sen.gene), ]
dir <- which(sen.gene %in% rownames(data))
if (length(dir) == 0 || length(dir) == 1)
    stop("too few senescence feature genes;cannot be calculated")
if (length(dir) > 1) {
    data <- data[match(sen.gene[dir], rownames(data)), ]
	p.non.sen.model <- non.sen.model[dir]
	p.non.sen.model <- as.numeric(p.non.sen.model)
	Scores <- c()
    for(i in 1:ncol(data)){
        Scores <- c(Scores, cor(data[,i], p.non.sen.model))
        }
	names(Scores) <- colnames(data)
    Scores <- cbind(Scores, names(Scores))
    Scores <- as.data.frame(Scores)
    colnames(Scores)[1:2] <- c("Score", "sample")
    Scores$Score <- as.numeric(Scores$Score)
	write.table(Scores, paste(cancer, "_senescence_score.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
```
