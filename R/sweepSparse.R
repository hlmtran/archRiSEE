#' fast row/column sweep for sparse matrices
#' 
#' @param x         a (sparse) matrix
#' @param MARGIN    margin to be swept
#' @param STATS     what to sweep out
#' @param FUN       how to sweep it out ('/')
#' 
#' @return          x, with stats swept out (divided through by default)
#'
#' @import          MatrixGenerics
#' @import          Matrix
#'
#' @export
#'
sweepSparse <- function(x, MARGIN, STATS, FUN="/") {

  f <- match.fun(FUN)
  if (MARGIN == 2) x <- t(x) 
  x@x <- f(x@x, STATS[x@i + 1])
  if (MARGIN == 2) x <- t(x) 
  return(x)

}
