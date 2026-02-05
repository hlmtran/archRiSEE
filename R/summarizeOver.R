#' given an SCE and a set of features, summarize counts over those features
#' 
#' @param x       a SummarizedExperiment, SingleCellExperiment, or similar
#' @param gr      a GRanges of features, or something coercible to GRanges
#' @param asy     which assay to summarize? (whichever one is first) 
#' @param BPPARAM BiocParallel parameters (default is MulticoreParam) 
#'
#' @return        a SummarizedExperiment holding the summarized counts
#'
#' @details       the result is an SE because you'll want to make it an altExp.
#'                you'll probably also want to run logNormCounts on the result.
#'
#' @seealso       scuttle::logNormCounts
#' @seealso       SSBtools::fac2sparse
#' 
#' @examples
#' 
#' \dontrun{
#'   library(rtracklayer) 
#'   K27feats <- import("NBMK27_50kbp_LSIfeats.hg38.bed")
#'   (MLLK27_50kb <- summarizeOver(MLLK27, K27feats))
#'   MLLK27_50kb <- logNormCounts(MLLK27_50kb)
#'   altExp(MLLK27, "K27_50kb_forNBM") <- MLLK27_50kb
#' }
#' 
#' @import        SSBtools
#'
#' @export 
#'
summarizeOver <- function(x, gr, asy=NULL) {

  how <- match.arg(how)
  if (is.null(asy)) asy <- assayNames(x)[1]
  stopifnot(is(x, "RangedSummarizedExperiment"))
  stopifnot(is(gr, "GRanges"))
  if (is.null(names(gr))) {
    names(gr) <- as.character(gr)
  }
  ol <- findOverlaps(rowRanges(x), gr)
  rn <- rownames(x)[queryHits(ol)]
  grn <- names(gr)[subjectHits(ol)]
  
  stop("Still not quite done, need to test other functions though")

  message("Summarizing assay ", asy, " over ", length(gr), " features...")
  summarized <- tapply(assay(x, asy)[rn,], grn, match.fun(how))
  message("Done.")
  res <- x
  assays(res) <- list()
  rowRanges(res) <- as(rownames(summarized, "GRanges"))
  assays(res, asy) <- summarized
  assayNames(res) <- asy
  return(res)

}
