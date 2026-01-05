#' block-randomized sampling to ensure somewhat representative comparisons 
#'
#' @param x       a rectangular object that inherits from SummarizedExperiment
#' @param g       grouping factor for each row/column, or its name in colData
#' @param maxN    maximum size per-subset? (Inf; use the smallest subset size)
#' @param INDEX   block-downsample on columns (2) or rows (1)? (2) 
#' @param quietly suppress warnings about unsampled metadata? (FALSE)
#' @param ...     additional parameters to pass to sample(), e.g. replace=TRUE
#' 
#' @details       use set.seed beforehand if you want reproducibility.
#'
#' @examples
#'
#' ncells <- 100
#' u <- matrix(rpois(20000, 5), ncol=ncells)
#' v <- log2(u + 1)
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(assays=list(counts=u, logcounts=v))
#' rownames(se) <- paste0("feature", seq_len(nrow(se)))
#' colnames(se) <- paste0("subject", seq_len(ncol(se)))
#' clusts <- factor(1:5) 
#' se$Clusters <- sample(clusts, ncol(se), replace=TRUE)
#' show(se) 
#'
#' # without replacement:
#' set.seed(1234)
#' (downsample(se, "Clusters"))
#'
#' # with replacement:
#' set.seed(1234)
#' (downsample(se, "Clusters", replace=TRUE))
#' 
#' @return        a subset of x, downsampled to match the smallest N(g)
#' 
#' @import        SummarizedExperiment
#'
#' @export
#'
downsample <- function(x, g, maxN=Inf, INDEX=2, quietly=FALSE, ...) {

  if (length(g) != dim(x)[INDEX]) {
    if (is(x, "SummarizedExperiment") & INDEX == 2) {
      if (g %in% names(colData(x))) g <- colData(x)[, g]
    } else {
      stop("dim(x)[",INDEX,"] == ", dim(x)[INDEX], " != length of ", g)
    }
  }

  gs <- levels(factor(g))
  names(gs) <- gs
  n <- min(maxN, min(table(g)))
  idx <- sort(do.call(c, lapply(gs, function(h) sample(which(g == h), n, ...))))
  xx <- switch(INDEX, x[idx, ], x[, idx])

  if (is(x, "SummarizedExperiment") & INDEX == 2) {
    mdn <- names(metadata(x))
    if (length(mdn) > 0 & !quietly) { 
      mdns <- paste(mdn, collapse=", ")
      warning("Metadata (", mdns, ") discarded. May be recoverable, e.g.")
      warning('xx <- downsample(x, "', g, '", INDEX=2)')
      warning('metadata(xx) <- lapply(metadata(x), "[", j=colnames(xx))')
      warning("This is not done automatically, as metadata is unrestricted.")
    }
  }

  return(xx)

}
