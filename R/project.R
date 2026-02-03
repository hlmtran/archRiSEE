#' project a (subset of an) SCE into a space defined by another
#'
#' @param x   SCE to be projected
#' @param y   SCE with projection features and saved UMAP model 
#' @param z   reducedDim name to store the result in object x
#' 
#' @return    x, with a new reducedDim
#'
#' @details   it's more fiddly than that, usually.  See examples.
#'
#' @import    uwot 
#' @import    RcppML
#'
#' @export
#'
project <- function(x, y, z="UMAP_projected") { 

  stop("Not quite done yet")

}
