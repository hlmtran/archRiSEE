#' given a matrix of data and its LSI representation, check cor(projected, LSI)
#'
#' @param mat   binarized tilematrix from an ArchR project, a project, or an SCE
#' @param LSI   if 'mat' is not an ArchR project, a complete ArchR LSI object
#' @param samp  how many cells to sample for the comparison? (1000)
#' @param MSE   return mean squared error instead of correlations? (FALSE) 
#' 
#' @return      column-wise correlations or diffs of projected data vs. matSVD
#'
#' @details     merge this into archRiSEE 
#'
#' @export
#'
checkProjection <- function(mat, LSI=NULL, samp=1000, MSE=FALSE) {

  if (!is(mat, "Matrix")) { # {{{
    stopifnot(require(ArchR))
    if (is(mat, "ArchRProject")) {
      LSI <- getReducedDims(mat, returnMatrix=FALSE)
      SE <- getMatrixFromProject(mat, useMatrix="TileMatrix", binarize=TRUE)
      mat <- assay(SE)
    } else if (is(mat, "SummarizedExperiment")) {
      LSI <- metadata(mat)$LSI
      stopifnot("svd" %in% names(LSI))
      mat <- assay(mat, "counts")
    } else {
      stop("Please provide a Matrix AND LSI object _OR_ ArchRProject _OR_ SCE!")
    }
  } # }}}

  nonOutliers <- setdiff(colnames(mat), LSI$outliers)
  samp <- min(length(nonOutliers), samp)
  message("Sampling ", samp, " random cells to re-project...")
  randomCells <- sample(nonOutliers, samp)
  original <- LSI$matSVD[randomCells, ]
  projected <- projectLSI(mat[, randomCells], LSI)
  stopifnot(ncol(projected) == ncol(original))
  
  LSIdims <- seq_len(ncol(projected))
  names(LSIdims) <- colnames(original)
  
  if (MSE) { 
    sapply(LSIdims, function(i) mean((projected[, i] - original[, i])**2))
  } else { 
    sapply(LSIdims, function(i) cor(projected[, i], original[, i]))
  }

}

