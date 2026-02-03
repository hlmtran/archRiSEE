#' stack (rbind) two SCE-like objects for the same cells to integrate
#' 
#' @param x     the first object (its rows will be on top)
#' @param y     the second object (its rows will be below y)
#' @param NMF   run NMF on the merged feature set? (TRUE) 
#' @param UMAP  run UMAP on the NMF results? (TRUE; see details)
#' @param k     rank of the NMF decomposition to shoot for (default: 30)
#' @param ...   currently unused parameters for convenience functions
#'
#' @param       an object of the same type as x and y, if all goes well
#'
#' @details     rbind()'ing objects is easy. Using the result isn't always so
#'              easy. We assume that people want to use NMF and UMAP to get a
#'              somewhat interpretable look at the combined data, but these are
#'              just two of many, many possible ways to approach the task. Most
#'              of the time, users will want to run archRtoSCE on a couple of 
#'              archR projects with how='feats' and feats=commonFeaturesGRanges
#'              to ensure that the desired features from both modalities are 
#'              present. After that the process tends to be iterative, which is
#'              the whole point of using iSEE to interact with the results.
#'
#' @import      SingleCellExperiment
#' @import      RcppML
#'
#' @export
#'
stack <- function(x, y, NMF=TRUE, UMAP=TRUE, k=30, ...) { 

  stop("stack() is not fully implemented yet") 

}
