#' reformat ArchR GeneScoreMatrix rowdata into a GRanges
#'
#' @param   se    the SummarizedExperiment holding the GeneScoreMatrix
#' 
#' @return        the same data as a RangedSummarizedExperiment
#'
#' @export
#'
archRgeneScoresWithRanges <- function(se, genome=NULL) {

  rd <- rowData(se)
  rownames(rd) <- rd$name
  rd$strand <- c('+','-')[rd$strand]
  starts <- pmin(as.integer(rd$start), as.integer(rd$end))
  ends <- pmax(as.integer(rd$start), as.integer(rd$end))
  rd$start <- starts
  rd$end <- ends
  rowRanges(se) <- makeGRangesFromDataFrame(rd)
  if (!is.null(genome)) genome(se) <- genome
  return(se) 

}
