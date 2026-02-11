#' ArchR-style log1p(TF-IDF) transformation but without mandatory binarization
#' 
#' @param mat       a (sparse) matrix of counts (cells as columns) or an SE
#' @param prune     prune (cells >= prune)? (10, by default; 0 to disable)
#' @param idf       a pre-trained idf or logTfIdf result (NULL; compute idf)
#' @param binarize  binarize the values for TF-IDF? (FALSE)
#' 
#' @return          a log1p(TF-IDF) version of the (counts) matrix (see Details)
#'
#' @details         TF == (term / sum(terms_in_document));
#'                  IDF == (1 / sum(documents_with_term_not_zero));
#'                  currently only assays(mat)[1] gets transformed for SE/SCEs;
#'                  provided or generated model idf placed in attr(mat, 'idf')
#'
#' @seealso         text2vec::TfIdf
#'
#' @import          Matrix
#'
#' @export
#'
logTfIdf <- function(mat, prune=10, idf=NULL, binarize=TRUE) { 

  if (is(mat, "SummarizedExperiment")) {
    return(logTfIdf(assay(mat), prune=prune, idf=idf, binarize=binarize))
  }
  if (is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  if (binarize) mat@x <- ifelse(mat@x > 0, 1, 0)
  origrows <- rownames(mat)

  # directly borrowed from ArchR
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  message("Term frequency normalization... ", appendLF=FALSE)
  mat@x <- (mat@x / rep.int(colSm, Matrix::diff(mat@p)))
  message("done.")

  if (!is.null(idf)) {
    message("Using precomputed IDF (inverse document frequency) table...")
    idfnames <- attr(idf, 'names')
    if (!all(idfnames == rownames(mat))) {
      stop("rownames(mat) != attr(idf, 'names'). Subset your matrix first.")
    }
  } else { 
    rowSm2 <- rowSums(mat > 0)
    toPrune <- rowSm2 < prune
    if (sum(toPrune) > 0) {
      message("Minimum observed feature frequency: ", min(rowSm2[rowSm2 > 0]))
      message("Minimum allowable feature frequency: ", prune)
      message("Zeroing out ",sum(toPrune)," of ",length(toPrune)," features.")
      keep <- as(Matrix::Diagonal(x=!toPrune), "sparseMatrix")
      rownames(keep) <- rownames(mat) 
      mat <- keep %*% mat
      rowSm2 <- rowSums(mat > 0)
      message("Minimum feature frequency is ", min(rowSm2[rowSm2 > 0]),
              " (", sum(keep), "/", length(toPrune), " features retained).")
    }
    message("Computing IDF (inverse document frequency) table... ", appendLF=0)
    idf <- as(ncol(mat)/rowSm2, "sparseVector")
    attr(idf, 'names') <- rownames(mat)
    message("done.")
  }
  mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
  message("Computing log1p(TF-IDF)... ", appendLF=FALSE)
  mat@x <- log(mat@x * 10001)
  rownames(mat) <- origrows
  attr(mat, 'idf') <- idf
  message("done.")
  return(mat)

}
