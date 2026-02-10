#' fast column sweep for sparse matrices
#' 
#' @param mat       a (sparse) matrix
#' @param stats     what to sweep out
#' @param fun       how to sweep it out ('/')
#' 
#' @return          mat, with colSums swept out (divided through, by default)
#'
#' @import          Matrix
#'
#' @export
#'
sweepSparse <- function(mat, stats, fun="/") {

  mat <- t(mat) 
  f <- match.fun(fun)
  mat@x <- f(mat@x, stats[mat@i + 1])
  return(t(mat))

}
