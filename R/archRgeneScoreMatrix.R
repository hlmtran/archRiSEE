#' reformat ArchR GeneScoreMatrix rowdata into a GRanges
#'
#' @param   proj  the ArchRproject holding the GeneScoreMatrix
#' 
#' @return        the gene score data as a RangedSummarizedExperiment
#'
#' @export
#'
archRgeneScoreMatrix <- function(proj, genome=NULL) {

  if (!"GeneScoreMatrix" %in% getAvailableMatrices(proj)) return(NULL)

  se <- getMatrixFromProject(proj, "GeneScoreMatrix") 
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
