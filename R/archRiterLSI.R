#' retrieve relevant IterativeLSI bits from an ArchR project
#' 
#' @param   proj    the ArchRProject
#' @param   LSIdim  name of the IterativeLSI reducedDims entry ("IterativeLSI")
#' @param   keep    these ("matSVD","outliers","LSIFeatures","u","d","v","idf")
#' @param   full    just pull the entire LSI object? (FALSE)
#'
#' @return  a list with embeddings, outliers to exclude, features GRanges, etc.
#'
#' @details archRiterLSI(proj, full=TRUE) is useful for projection testing.
#'
#' @export
#'
archRiterLSI <- function(proj, LSIdim="IterativeLSI", keep=NULL, full=FALSE) {

  stopifnot(require(ArchR))
  LSI <- getReducedDims(proj, LSIdim, returnMatrix=FALSE)
  width <- archRtileSize(proj)        
  LSIF <- LSI[["LSIFeatures"]]
  LSIGR <- GRanges(seqnames=LSIF$seqnames, 
                   ranges=IRanges(start=LSIF$start, 
                                  width=rep(width, nrow(LSIF))),
                   rowSums=LSIF$rowSums)
  names(LSIGR) <- as.character(LSIGR)
  genome(LSIGR) <- archRgenome(proj)
  LSI$idf <- LSI$nCol/LSI$rowSm
  names(LSI$idf) <- rownames(LSI)
  if (full) return(LSI) 
  
  res <- list() 
  if (is.null(keep)) {
    keep <- c("matSVD", "outliers", "LSIFeatures", "u", "d", "v", "idf")
  }
  oddities <- c("LSIFeatures", "u", "d", "v", "idf")
  for (k in setdiff(keep, oddities)) res[[k]] <- LSI[[k]]
  res$LSIFeatures <- sort(sortSeqlevels(LSIGR))

  # extract right and left singular vectors and the diagonal
  if (any(c("u","d","v") %in% keep)) { 
    for (k in c("u", "d", "v")) res[[k]] <- LSI[["svd"]][[k]]
  }

  # idf is special, recompute it on the fly:
  if ("idf" %in% keep) res[["idf"]] <- LSI$nCol/LSI$rowSm


  # done
  return(res)

}
