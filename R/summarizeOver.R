#' given an SCE and a set of features, summarize counts over those features
#' 
#' @param x       a SummarizedExperiment, SingleCellExperiment, or similar
#' @param gr      a GRanges of features, or something coercible to GRanges
#' @param asy     which assay to summarize? (whichever one is first) 
#' @param lognorm log-normalize the summarized counts? (TRUE) 
#'
#' @return        a SummarizedExperiment holding the summarized counts
#'
#' @details       the result is an SE because you'll want to make it an altExp.
#'
#' @seealso       scuttle::logNormCounts
#' @seealso       Matrix::fac2sparse
#' 
#' @examples
#' 
#' \dontrun{
#'   library(rtracklayer) 
#'   K27feats <- import("NBMK27_50kbp_LSIfeats.hg38.bed")
#'   MLLK27 <- altExp(stacked, "MLLK27.FragmentCounts")
#'   rownames(MLLK27) <- as.character(granges(MLLK27))
#'   (MLLK27_50kb <- summarizeOver(MLLK27, K27feats)))
#'   MLLK27_50kb <- logNormCounts(MLLK27_50kb)
#'   altExp(MLLK27, "K27_50kb_forNBM") <- MLLK27_50kb
#' }
#' 
#' @import        Matrix 
#'
#' @export 
#'
summarizeOver <- function(x, gr, asy=NULL, lognorm=TRUE, ...) {

  # annoying, but mandatory 
  x <- .checkRowRanges(x)
  gr <- .checkFeatureRanges(gr) 
  if (is.null(asy)) asy <- assayNames(x)[1]

  x <- subsetByOverlaps(x, gr)
  ol <- findOverlaps(rowRanges(x), gr)
  mult <- fac2sparse(names(gr)[subjectHits(ol)]) 
  message("Summarizing assay ", asy, " over ", length(gr), " features...")
  summarized <- (mult %*% assay(x, asy)[queryHits(ol), ])[names(gr), ]
  stopifnot(identical(rownames(summarized), names(gr)))
  message("Done.")
  res <- x[seq_along(gr), ]
  rowRanges(res) <- as(rownames(summarized), "GRanges")
  rownames(res) <- as(rowRanges(res), "character")
  assays(res) <- list(summarized)
  assayNames(res) <- asy

  if (lognorm) { 
    message("Log-normalizing summarized ", asy, "...")
    res <- logNormCounts(res, assay.type=asy)
    message("Done.")
  }

  return(res)

}


# helper fn
.checkRowRanges <- function(x) { 

  stopifnot(is(x, "SummarizedExperiment"))
  if (!is(x, "RangedSummarizedExperiment")) {
    xx <- as(rownames(x), "GRanges")
    names(xx) <- as.character(xx)
    if (is(xx, "GRanges")) {
      rowRanges(x) <- xx
    } else { 
      stop("Cannot determine the coordinates of rows in `x`")
    }
  } else {
    # if there's a silly GRangesList in there, collapse it out
    if (is(rowRanges(x), "GRangesList")) { 
      if (length(unlist(rowRanges(x))) == nrow(x)) { 
        rowRanges(x) <- unlist(rowRanges(x)) # rare
      } else if (!is.null(rownames(x))) {
        rowRanges(x) <- as(rownames(x), "GRanges") 
      } else { 
        stop("Cannot summarize without a usable rowRanges!")
      }
    }
  }
  stopifnot(is(rowRanges(x), "GRanges"))
  rownames(x) <- as(rowRanges(x), "character")
  return(x) 

}


# helper fn
.checkFeatureRanges <- function(gr) {

  stopifnot(is(gr, "GRanges"))
  if (is.null(names(gr))) names(gr) <- as.character(gr)
  return(gr) 

}
