#' figure out which genome was used by ArchR
#'
#' @param proj          an ArchRproject 
#'
#' @return              a character string indicating the genome
#'
archRgenome <- function(proj) { 

  frags <- strsplit(proj@genomeAnnotation$genome, "\\.")[[1]]
  frags[length(frags)] # usually hg38 or hg19

}
