#' add an NMF decomposition to the reducedDims and rowRanges of an SCE
#'
#' By default, the weights for each factor are copied to mcols(rowRanges(x)).
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param k             rank for the NMF decomposition on TFIDF (30) 
#' @param colDat        add NMF scores to colData(x)? (TRUE)
#' @param rowDat        add NMF weights to mcols(rowRanges(x))? (TRUE)
#' @param nmf_fit       an existing nmf fit to project the data (NULL) 
#' @param ...           additional arguments to pass to RcppML::nmf()
#'
#' @return              SingleCellExperiment with reducedDim(x, "NMF")
#'
#' @details             the mcols(rowRanges(x)) weights are for igvR plotting
#'                      it would be a good idea to flag a "sensible" rank here
#'
#' @seealso             archRtoSCE
#' @seealso             iSEEarchR
#' @seealso             addIgvTrack
#'
#' @import RcppML
#' @import scuttle
#' @import SingleCellExperiment
#'
#' @export
#'
addNMF <- function(x, k=30, colDat=TRUE, rowDat=TRUE, nmf_fit=NULL, ...) { 

  orig <- options()[["RcppML.verbose"]]
  if (!"TFIDF" %in% assayNames(x)) x <- addTfIdf(x)

  if (is.null(nmf_fit)) {
    message("Fitting rank-", k, " NMF model on assay(x, 'TFIDF')...")
    options("RcppML.verbose" = TRUE)
    nmf_fit <- RcppML::nmf(assay(x, 'TFIDF'), k=k, ...)
  } else { 
    nmf_fit@h <- RcppML::predict(nmf_fit, data=assay(x, 'TFIDF'), ...)
  }
  message("Saving model to metadata(x)$NMF...")
  metadata(x)$NMF <- nmf_fit
  message("Copying NMF hat matrix to reducedDim(x, 'NMF')...")
  reducedDim(x, "NMF") <- t(metadata(x)$NMF@h)
  NMFdims <- ncol(reducedDim(x, "NMF"))
  names(reducedDim(x, "NMF")) <- paste0("NMF", seq_len(NMFdims))
  if (colDat) { 
    message("Copying NMF hat columns to colData() for iSEE visualization...")
    for (i in colnames(reducedDim(x, "NMF"))) {
      message("Adding ", i, " as colData(x)$", toupper(i))
      colData(x)[, toupper(i)] <- reducedDim(x, "NMF")[, i]
    }
  }
  if (rowDat) {
    message("Copying NMF weights to mcols(rowRanges(x)) for igvR plotting...")
    for (i in colnames(metadata(x)$NMF@w)) { 
      message("Adding NMF@w[, ", i, "] as mcols(rowRanges(x))$", toupper(i))
      mcols(x)[, toupper(i)] <- metadata(x)$NMF@w[, i]
    }
  }
  options("RcppML.verbose" = orig)
  return(x) 

}
