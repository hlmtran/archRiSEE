#' flexibly add clusterings with a variety of methods and benchmarks
#' 
#' @param   x           a SingleCellExperiment object or something like it 
#' @param   BLUSPARAM   a BlusterParam object (defaults are often workable)
#' @param   rdName      which reducedDim to use? (the first one) 
#' 
#' @return              a factor of length(nrow(reducedDim
#' 
#' @details this function delegates to clusterRows(), which as the name
#'          suggests, works on a matrix where the observations are rows and 
#'          the measurements are columns (i.e., a reducedDim). If no BLUSPARAM
#'          is specified, the function defaults to TwoStepParam().
#'
#' @seealso bluster::clusterRows
#' @seealso bluster::clusterSweep
#' @seealso bluster::approxSilhouette
#' @seealso bluster::bootstrapStability
#' @seealso bluster::compareClusterings
#' @seealso bluster::linkClusters
#' @seealso bluster::TwoStepParam-class
#' @seealso bluster::NNGraphParam-class
#' @seealso bluster::DbscanParam-class
#' 
#' @import bluster
#' 
#' @export
#'
addBlusters <- function(x, BLUSPARAM=NULL, rdName=NULL) { 

  if (is.null(BLUSPARAM)) BLUSPARAM <- TwoStepParam()
  if (is.null(rdName)) rdName <- reducedDimNames(x)[1]

  message("Running bluster::clusterRows on ", rdName, 
          " with ", class(BLUSPARAM)[[1]])
  clusterRows(reducedDim(x, rdName), BLUSPARAM=BLUSPARAM)

}
