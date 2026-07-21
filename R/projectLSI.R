#' stolen directly from ArchR, with simplifications for the default settings
#'
#' @param   mat     a sparse Matrix, M rows x N columns
#' @param   LSI     an ArchR LSI object with nDimensions = K
#'
#' @return          an N x K projection matrix 
#'
#' @import  Matrix
#'
#' @export
#'
projectLSI <- function(mat, LSI) { 

  mat <- mat[LSI$idx,]                                          # check idx?
  mat@x[mat@x > 0] <- 1                                         # binarize
  colSm <- Matrix::colSums(mat)                                 # check for 0
  if(any(colSm == 0)){                                          # exclude zeros 
    exclude <- which(colSm==0)
    mat <- mat[,-exclude]
    colSm <- colSm[-exclude]
  }
  stopifnot(all(colSm > 0))                                     # fail on zeros
  mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))          # TF normalize
  idf <- as(LSI$nCol / LSI$rowSm, "sparseVector")               # pull the IDF
  mat <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat # reweight mat
  mat@x <- log(mat@x * LSI$scaleTo + 1)                         # log-TF-IDF
  idxNA <- Matrix::which(is.na(mat), arr.ind=TRUE)              # add a warning!
  if(length(idxNA) > 0) mat[idxNA] <- 0                         # patch any NAs
  V <- Matrix::t(mat) %*% LSI$svd$u %*% diag(1/LSI$svd$d)       # eigenvectors  
  svdDiag <- matrix(0,nrow=LSI$nDimensions,ncol=LSI$nDimensions)# D in UDt(V)
  diag(svdDiag) <- LSI$svd$d                                    # as above
  matSVD <- as.matrix(Matrix::t(svdDiag %*% Matrix::t(V)))      # project mat
  rownames(matSVD) <- colnames(mat)                             # cell names
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))       # LSIdim names
  return(matSVD)                                                # if all is well

}
