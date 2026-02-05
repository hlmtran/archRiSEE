#' add ArchR peakSet ranges to a PeakMatrix-exported SE
#'
#' @param   proj    the ArchR project, which MUST HAVE A peakSet (!)
#' @param   genome  the genome to annotate the peakSet from (NULL)
#' 
#' @return          a RangedSummarizedExperiment-housed PeakMatrix, or NULL
#'
#' @export
#'
archRpeakMatrix <- function(proj, genome=NULL) {

  if (is.null(getPeakSet(proj))) return(NULL)
  if (!"PeakMatrix" %in% getAvailableMatrices(proj)) return(NULL)

  se <- getMatrixFromProject(proj, "PeakMatrix") 
  rownames(se) <- as.character(rowRanges(se))
  rowRanges(se)$idx <- NULL
  if (!is.null(genome)) genome(se) <- genome
  return(se) 

}
