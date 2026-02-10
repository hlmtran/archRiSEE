#' add ArchR-style log1p(TF-IDF) transformed counts matrix to assays(x, 'TFIDF')
#' 
#' @param x         a SummarizedExperiment
#' @param idf       a pre-trained idf or logTfIdf result (NULL; compute idf)
#' @param ...       additional arguments passed to logTfIdf()
#' 
#' @return          x, but with a TFIDF matrix added to assays(x)
#'
#' @details         the 'idf' parameter for logTfIdf() is of interest, as
#'                  it will be propagated into metadata(x)$idf and can be
#'                  passed in to addTfIdf as well (e.g. if projecting data)
#'
#' @seealso         logTfIdf
#'
#' @import          SummarizedExperiment
#'
#' @export
#'
addTfIdf <- function(x, idf=NULL, ...) { 

  assay(x, "TFIDF") <- logTfIdf(x, idf=idf, ...)
  message("Adding idf to metadata(x)$idf...")
  metadata(x)$idf <- attr(assay(x, "TFIDF"), 'idf')
  return(x)

}
