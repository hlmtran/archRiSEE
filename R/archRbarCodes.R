#' extract unique barcodes from cellColData in an archR project
#' 
#' @param proj    the project
#' @param splt    the split string ('#')
#' @param pos     which piece to keep (2)
#'
#' @return        a character vector
#'
#' @examples
#'
#' inDEM <- archRbarCodes(MLLDEM)
#' inK27 <- archRbarCodes(MLLK27)
#' table(inDEM %in% inK27)
#' table(inK27 %in% inDEM) 
#'
#' @export
#'
archRbarCodes <- function(proj, splt="#", pos=2)  { 

  if (any(!grepl(splt, rownames(proj@cellColData)))) {
    return(rownames(proj@cellColData)) 
  } else { 
    sapply(strsplit(rownames(proj@cellColData), splt), "[", pos)
  }

}
