#' ArchR-style log1p(TF-IDF) transformation but without mandatory binarization
#' 
#' @param mat       a (sparse) matrix of counts (cells as columns) or an SE
#' @param prune     prune (cells >= prune)? (10, by default; 0 to disable)
#' @param scaleTo   scaling factor for library sizes (10000, i.e., cp10k)
#' @param idf       a pre-trained idf or logTfIdf result (NULL; compute idf)
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
logTfIdf <- function(mat, prune=10, scaleTo=10000, idf=NULL) { 

  if (is(mat, "SummarizedExperiment")) {
    return(logTfIdf(assay(mat), prune=prune, scaleTo=scaleTo, idf=idf))
  }
  if(is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  message("Computing TF (term frequency) for ", nrow(mat), 
          " features across ", ncol(mat), " cells...")
  mat <- sweepSparse(mat, colSums(mat)) * scaleTo 

  if (!is.null(idf)) {
    message("Using precomputed IDF (inverse document frequency) table...")
    idfnames <- attr(idf, 'names')
    if (!all(idfnames == rownames(mat))) {
      stop("rownames(mat) != attr(idf, 'names'). Subset your matrix first.")
    }
  } else { 
    rowSm <- rowSums(mat > 0)
    toPrune <- rowSm < prune
    if (sum(toPrune) > 0) {
      message("Minimum observed feature frequency: ", min(rowSm[rowSm>0]))
      message("Minimum allowable feature frequency: ", prune)
      message("Zeroing out ",sum(toPrune)," of ",length(toPrune)," features.")
      keep <- as(Diagonal(x=!toPrune), "sparseMatrix")
      rownames(keep) <- rownames(mat) 
      mat <- keep %*% mat
      rowSm <- rowSums(mat > 0)
      message("Minimum feature frequency is ", min(rowSm[rowSm > 0]),
              " (", sum(keep), "/", length(toPrune), " features retained).")
    }
    message("Computing IDF (inverse document frequency) table...")
    idf <- as((ncol(mat) + 1)/(rowSums(mat > 0) + 1), "sparseVector")
    attr(idf, 'names') <- rownames(mat)
  }
  logtfidf <- tfidf <- (as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat)
  message("Computing log1p(TF-IDF)...", appendLF=FALSE)
  logtfidf@x <- log1p(tfidf@x)
  attr(logtfidf, 'idf') <- idf
  rownames(logtfidf) <- rownames(mat)
  message("Done.")
  return(logtfidf)

}
