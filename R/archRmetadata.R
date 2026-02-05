#' compile ArchR metadata to pack into an SCE
#'
#' @param proj          an ArchRproject 
#' @param lst           slots to mine for metadata (see details) 
#'
#' @return              a list of metadata items
#'
#' @details             defaults to mining 'projectMetadata',
#'                      'sampleMetadata', and 'cellMetadata'
#'
#' @export
#'
archRmetadata <- function(proj, lst=NULL) { 

  require(ArchR)
  def <- c('project','sample','cell')
  if (is.null(lst)) lst <- setNames(paste0(def, "Metadata"), def)
  as.list(unlist(lapply(lst, function(x) 
                        do.call(c, lapply(slot(proj, x), unlist)))))

}
