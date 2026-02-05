#' figure out which genome was used by ArchR
#'
#' @param proj          an ArchRproject 
#' @param human         human genome? (TRUE; see details)
#'
#' @return              a character string indicating the genome
#'
#' @details             this function assumes BSgenome.Hsapiens.UCSC.* was used
#'
#' @export
#'
archRgenome <- function(proj, human=TRUE) { 

  require(ArchR)
  if (human) sub("BSgenome.Hsapiens.UCSC.", "", getGenome(proj))
  else stop("Currently we only handle human genome annotations. Send a PR!") 

}
