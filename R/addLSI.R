#' add an LSI result to an SCE (NOTE: this assumes prefiltering was done!) 
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF") 
#' @param name          name of reducedDim to store the LSI embedding ("LSI") 
#' @param scaleDims     Z-scale the dimensions? (TRUE)
#' @param corCutOff     drop dimensions where cor(LSIdim, nFrag) > corCutOff
#' @param excludeChr    chroms to exclude by default (chrM, chrX, chrY) 
#' @param nDimensions   number of dimensions for the SVD (30) 
#' @param depth         colData column holding depth for correlations ('nFrags')
#' @param subsetLSI     subset to mcols(x)$usedForLSI? (TRUE) 
#'
#' @return              SingleCellExperiment with reducedDim(x, name)
#'
#' @details             This function assumes that ArchR selected the features!
#'                      Here we do a simple, noniterative one-pass SVD on TFIDF
#'
#' @seealso             addTfIdf
#' @seealso             ArchR::addIterativeLSI
#'
#' @import              SingleCellExperiment
#' @import              GenomeInfoDb
#' @import              irlba
#'
#' @export
#'
addLSI <- function(x, useMatrix=c("TFIDF","counts"), name="LSI", scaleDims=TRUE, corCutOff=0.75, excludeChr=c("chrM","chrX","chrY"), nDimensions=30, depth="nFrags", subsetLSI=TRUE, ...) { 

  stopifnot(depth %in% names(colData(x)))
  useMatrix <- match.arg(useMatrix) 
  keep <- setdiff(seqlevels(x), excludeChr)
  if (subsetLSI) {
    stopifnot("usedForLSI" %in% names(mcols(x)))
    idx <- which(mcols(x)$usedForLSI)
  } else {
    idx <- seq_len(nrow(x))
  }
  message("Subsetting TF-IDF matrix...") 
  mat <- assay(keepSeqlevels(x[idx,], keep, pruning.mode="coarse"), useMatrix)
  message("Running SVD...")
  svd <- irlba::irlba(mat, nDimensions + 5, nDimensions + 5)
  svdDiag <- matrix(0, nrow=nDimensions + 5, ncol=nDimensions + 5)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  message("Checking for depth-correlated columns...")
  toKeep <- which(cor(matSVD, colData(x)[[depth]])[, 1] < corCutOff)
  matSVD <- matSVD[, toKeep][, seq_len(nDimensions)]
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  reducedDim(x, name) <- matSVD
  message("Done.")
  return(x) 

}
