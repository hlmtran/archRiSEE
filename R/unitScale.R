#' scale factors (e.g. NMF factors) to unit distance
#' 
#' @param x         a (sparse) matrix
#' @param MARGIN    1: rows, 2: columns (1; unnecessary for row/colUnitScale()) 
#' 
#' @return          x, scaled 0-1
#'
#' @export
#'
unitScale <- function(x, MARGIN=1) {

  stopifnot(MARGIN %in% 1:2)
  STAT1 <- ifelse(MARGIN == 1, rowMins(x, na.rm=TRUE), colMins(x, na.rm=TRUE))
  STAT2 <- ifelse(MARGIN == 1, rowMaxs(x, na.rm=TRUE), colMaxs(x, na.rm=TRUE))
  FN <- ifelse(is(x, "sparseMatrix"), archRiSEE::sweepSparse, base::sweep)
  x <- FN(x, MARGIN, STAT1, "-")
  x <- FN(x, MARGIN, STAT2, "/")
  return(x)

}


#' @rdname unitScale 
#'
#' @export
#'
rowUnitScale <- function(x, ...) unitScale(x, MARGIN=1)


#' @rdname zScore
#'
#' @export
#'
colUnitScale <- function(x, ...) unitScale(x, MARGIN=2)
