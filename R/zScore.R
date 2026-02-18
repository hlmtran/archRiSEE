#' scale factors (e.g. PCA) to Z-scores
#' 
#' @param x         a (sparse) matrix
#' @param MARGIN    1: rows, 2: columns (1; unnecessary for row/colZscores()) 
#' 
#' @return          x, Z-scored
#'
#' @export
#'
zScore <- function(x, MARGIN=1) {

  stopifnot(MARGIN %in% 1:2)
  STAT1 <- ifelse(MARGIN == 1, rowMeans(x, na.rm=TRUE), colMeans(x, na.rm=TRUE))
  STAT2 <- ifelse(MARGIN == 1, rowSds(x, na.rm=TRUE), colSds(x, na.rm=TRUE))
  FN <- ifelse(is(x, "sparseMatrix"), archRiSEE::sweepSparse, base::sweep)
  x <- FN(x, MARGIN, STAT1, "-")
  x <- FN(x, MARGIN, STAT2, "/")
  return(x)

}


#' @rdname zScore
#'
#' @export
#'
rowZscores <- function(x, ...) zScore(x, MARGIN=1)


#' @rdname zScore
#'
#' @export
#'
colZscores <- function(x, ...) zScore(x, MARGIN=2)
