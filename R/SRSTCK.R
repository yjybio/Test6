#' SRSTCK main_func.
#'
#' @param data raw data matrix; if the data is scRNA-seq, genes in rows; cell names in columns; if the data is bulk, genes in rows; sample names in columns.
#' @param cancer short for cancer: such as "OVCA";"CRC";SKCM" etc.
#' @param type if the data is scRNA-seq, put "sc"; if the data is bulk RNA-seq, put "bulk".
#' @param cell_type if the data is scRNA-seq, "cell_type" represents a collection of TME cells, including immune cells, fibroblasts, and endothelial cells. Please convert these cell types in your data into a vector and pass it to "cell_type".
#' @param multi.sample if scRNA-seq is a single sample, put "FALSE"; if scRNA-seq is multi-sample integration, put "TRUE".
#' @param minGene the minimum number of genes in a cell.
#' @param maxGene the maximum number of genes in a cell.
#' @param maxUMI the maximum UMI count in a cell.
#' @param pctMT mitochondrial gene ratio.
#' @param id.type gene id type: Symbol or Ensemble.
#' @param cell.line if the data are from pure cell line,put "yes"; if cell line data are a mixture of tumor and normal cells, still put "no".
#' @param ngene.chr minimal number of genes per chromosome for cell filtering.
#' @param win.size minimal window sizes for segmentation.
#' @param KS.cut segmentation parameters, input 0 to 1; larger looser criteria.
#' @param distance distance methods include euclidean, and correlation converted distance include pearson and spearman.
#' @param n.cores number of cores for parallel computing.
#'
#' @return if data is scRNA-seq: 1) aneuploid/diploid prediction results; 2)senescent tumor cells prediction results; 3)annotated cluster diagram of senescent tumor cells; if data is bulk: the senescence score results.
#' @export
#'
#' @examples
#' TME_cell_types <- c("Fibroblast", "Endothelial_cell", "T_NK", "B_cell", "Myeloid_cell")
#' test.SRSTCK <- SRSTCK(data = data, cancer = "CRC", type = "sc", cell_type = TME_cell_types, multi.sample = "FALSE", minGene = 500, maxGene = 5000, maxUMI = 40000, pctMT = 20, id.type = "S", cell.line = "no", ngene.chr = 5,  win.size = 25, KS.cut = 0.15, distance = "euclidean", n.cores = 1)
SRSTCK <- function (data, cancer, type, cell_type, multi.sample = "FALSE", minGene = 500,
                    maxGene = 5000, maxUMI = 40000, pctMT = 20, id.type = "S",
                    cell.line = "no", ngene.chr = 5, win.size = 25, KS.cut = 0.15,
                    distance = "euclidean", n.cores = 1)
{ 
  file_path <- system.file("Result_Data", "Model.xlsx", package = "SRSTCK")
  Model <- readxl::read_excel(file_path)
  Model <- as.data.frame(Model)
  rownames(Model) <- Model[, 1]
  Model <- Model[, -1]
  sen.gene <- rownames(Model)
  sen.model <- Model$sen_model
  non.sen.model <- Model$non_sen_model

  if (type == "sc") {
    print("Step1: Predictive tumor cells ...")
	scRNA_tumor_cell <- as.matrix(data@assays$RNA@counts)
    scRNA_tumor_cell <- apply(scRNA_tumor_cell, 2, as.numeric)
    rownames(scRNA_tumor_cell) <- rownames(data@assays$RNA@counts)
    colnames(scRNA_tumor_cell) <- colnames(data@assays$RNA@counts)
	
	if (!"Cell_Type" %in% colnames(data@meta.data)){
    stop("Error: The input data does not contain a column named 'Cell_Type'.\n",
         "Please rename your cell annotation column to 'Cell_Type' and try again.\n",
         "Current column names in your data are: ", paste(colnames(data@meta.data), collapse = ", "))
    }
	TME_cells <-  colnames(data)[data$Cell_Type %in% cell_type]
	
    if(length(TME_cells) >= 50) {
      copykat.test <- copykat(rawmat = scRNA_tumor_cell, sam.name = cancer, norm.cell.names = TME_cells)

    } else {
    warning("The number of TME cells is insufficient. It is recommended to check the cell type annotation")
      copykat.test <- copykat(rawmat = scRNA_tumor_cell, sam.name = cancer)
    }
    sample.name <- paste(cancer, "_copykat_", sep = "")
    res <- paste(sample.name, "prediction.txt", sep = "")
    copykat.prediction <- read.table(res)
    copykat.prediction <- copykat.prediction[which(copykat.prediction[, 2] == "aneuploid"), ]
    copykat.prediction <- data@assays$RNA@counts[ ,which(colnames(data@assays$RNA@counts) %in% copykat.prediction[, 1])]
    
	print("Step2: Predicting STCs ...")
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
    data@meta.data$sen_type <- data@meta.data$Cell_Type
    data@meta.data$sen_type <- as.character(data@meta.data$sen_type)
    data@meta.data$sen_type[which(rownames(data@meta.data) %in% sen_Cancer$Cell_ID)] <- sen_Cancer$State
    data@meta.data$sen_type <- as.factor(data@meta.data$sen_type)
    saveRDS(data, file = paste(cancer, "_SRSTCK_sc_prediction.rds", sep = ""))
    pdf(paste(cancer, "_note_STCs.pdf", sep = ""), height = 6, width = 8)
    p1 <- DimPlot(data, reduction = "umap", label = TRUE, group.by = "sen_type", label.size = 3)
    print(p1)
    dev.off()
	return(sen_Cancer)
  }
  
  if (type == "bulk") {
    print("Read and analyze bulk RNA-seq Data ...")
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
		return(Scores)
      }
    }
  }

