#' retrieve relevant IterativeLSI bits from an ArchR project
#' 
#' @param   proj    the ArchRProject
#' @param   LSIdim  name of the IterativeLSI reducedDims entry ("IterativeLSI")
#' @param   keep    these ("matSVD","outliers","LSIFeatures","u","d","v","idf")
#'
#' @return  a list with embeddings, outliers to exclude, features GRanges, etc.
#'
#' @import Seqinfo
#'
#' @export
#'
archRiterLSI <- function(proj, LSIdim="IterativeLSI", keep=NULL) {

  if (is.null(keep)) {
    keep <- c("matSVD", "outliers", "LSIFeatures", "u", "d", "v", "idf")
  }
  oddities <- c("LSIFeatures", "u", "d", "v", "idf")

  res <- list() 
  for (k in setdiff(keep, oddities)) res[[k]] <- .LSI(proj, LSIdim)[[k]]

  # convert LSIFeatures to a GRanges:
  if ("LSIFeatures" %in% keep) {
    LSIF <- .LSI(proj, LSIdim)[["LSIFeatures"]]
    # in order to make this usable, we need to find out the tile size:
    width <- archRtileSize(proj)        
    LSIGR <- with(LSIF, 
                  GRanges(seqnames=seqnames, 
                          ranges=IRanges(start=start, width=width)))
    names(LSIGR) <- as.character(LSIGR)
    LSIGR$rowSums <- LSIF$rowSums
    genome(LSIGR) <- archRgenome(proj)
    si <- Seqinfo(genome=unique(genome(LSIGR)))
    suppressWarnings(seqinfo(LSIGR) <- si[seqlevels(LSIGR)])
    res$LSIFeatures <- trim(sort(sortSeqlevels(LSIGR)))
  }

  # extract right and left singular vectors and the diagonal
  if (any(c("u","d","v") %in% keep)) { 
    for (k in c("u", "d", "v")) res[[k]] <- .LSI(proj, LSIdim)[["svd"]][[k]]
  }

  # idf is special, recompute it on the fly:
  if ("idf" %in% keep) res[["idf"]] <- with(.LSI(proj, LSIdim), nCol/rowSm)
 
  # done
  return(res)

}



# helper fn
.LSI <- function(proj, LSIdim="IterativeLSI") slot(proj,"reducedDims")[[LSIdim]]
