#' retrieve relevant IterativeLSI bits from an ArchR project
#' 
#' @param   proj    the ArchRProject
#' @param   LSIdim  name of the IterativeLSI reducedDims entry ("IterativeLSI")
#' @param   keep    these ("matSVD","outliers","LSIFeatures","u","d","v","idf")
#'
#' @return  a list with embeddings, outliers to exclude, features GRanges, etc.
#'
#' @export
#'
archRiterLSI <- function(proj, LSIdim="IterativeLSI", keep=NULL) {

  if (is.null(keep)) {
    keep <- c("matSVD", "outliers", "LSIFeatures", "u", "d", "v", "idf")
  }
  oddities <- c("LSIFeatures", "u", "d", "v", "idf")

  res <- list() 
  LSI <- .LSI(proj, LSIdim)
  for (k in setdiff(keep, oddities)) res[[k]] <- LSI[[k]]

  # convert LSIFeatures to a GRanges:
  if ("LSIFeatures" %in% keep) {
    width <- archRtileSize(proj)        
    LSIF <- LSI[["LSIFeatures"]]
    LSIGR <- GRanges(seqnames=LSIF$seqnames, 
                     ranges=IRanges(start=LSIF$start, 
                                    width=rep(width, nrow(LSIF))),
                     rowSums=LSIF$rowSums)
    names(LSIGR) <- as.character(LSIGR)
    genome(LSIGR) <- archRgenome(proj)
    res$LSIFeatures <- sort(sortSeqlevels(LSIGR))
  }

  # extract right and left singular vectors and the diagonal
  if (any(c("u","d","v") %in% keep)) { 
    for (k in c("u", "d", "v")) res[[k]] <- LSI[["svd"]][[k]]
  }

  # idf is special, recompute it on the fly:
  if ("idf" %in% keep) res[["idf"]] <- LSI$nCol/LSI$rowSm
 
  # done
  return(res)

}



# helper fn
.LSI <- function(proj, LSIdim="IterativeLSI") slot(proj,"reducedDims")[[LSIdim]]
