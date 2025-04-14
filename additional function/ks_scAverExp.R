ks_scAverExp <- function(obj,#one of seurat obj or Archr obj
                         features,#genes which your want to plot
                         CellTypes){ #cell annotation
  
  
  dataType = class(obj)
  
  if(dataType == "Seurat"){
    require(Seurat)
    # Dot plot of marker genes 
    count_mat <- GetAssayData(object=obj, slot="counts")
    
    min_pct = 5 
  }
  
  
  if(dataType == "ArchRProject"){
    require(ArchR)
    # Dot plot of GeneScoreMatrix cluster markers
    GSM_se <- getMatrixFromProject(obj, useMatrix="GeneScoreMatrix")
    count_mat <- assays(GSM_se)$GeneScoreMatrix
    rownames(count_mat) <- rowData(GSM_se)$name
    
    min_pct = 10
  }
  
  
  #average
  
  sparse_Matrix <- t(t(count_mat)/Matrix::colSums(count_mat)) * 10^4
  
  if(is(sparse_Matrix, "sparseMatrix")){
    matsum <- NMF::summary(sparse_Matrix) 
    logx <- log2(matsum$x + 1)
    logmat <- Matrix::sparseMatrix(i = matsum$i, j = matsum$j, 
                                   x = logx, dims = dim(sparse_Matrix),
                                   dimnames = dimnames(sparse_Matrix))
  }else{
    
    logmat <- log2(sparse_Matrix + 1) 
    
  }
  
  
  final_Matrix <- lapply(CellTypes, function(x) {
    
    rowMeans(logmat[, which(obj$celltype == x), drop = F], na.rm = TRUE)
    
  }) %>% Reduce("cbind", .)
  
  
  colnames(final_Matrix) <- CellTypes
  final_Matrix <- log2(final_Matrix + 1)
  
  avgExpr <- final_Matrix[rowSums(final_Matrix) > 0,]
  avgExpr <- avgExpr / apply(avgExpr,1,max)
  
  
  #pct
  groupFun <- function (mat, fun, groups = NULL, ...){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
      fun(mat[, which(groups == x), drop = F], ...)
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
  }
  
  
  pctExpr <- function(mat){
    nexpr <- apply(mat, 1, function(x) sum(x > 0))
    nexpr / ncol(mat)
  }
  
  pctExprMat <- groupFun(count_mat, pctExpr, groups=obj$celltype) * 100
  pctExprMat <- pctExprMat[apply(pctExprMat, 1, function(x) max(x) > 0),]
  
  
  expMelt <- reshape2::melt(avgExpr)
  colnames(expMelt) <- c("feature", "group", "avgExpr")
  pctMelt <- reshape2::melt(pctExprMat)
  colnames(pctMelt) <- c("feature", "group", "pctExpr")
  
  avgPctMat <- base::merge(expMelt, pctMelt, by=c("feature", "group"))
  
  subGenes <- featureSets %>% do.call("c",.)
  avgPctMat <- avgPctMat[avgPctMat$feature %in% subGenes,]
  
  avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- min_pct
  
  return(avgPctMat)
  
  
}

