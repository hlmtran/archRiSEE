#' convenience wrapper for iSEE on an ArchR-derived singlecellexperiment
#'
#' Note that this is a quick and dirty affair so don't expect much. 
#'
#' @param x             a SingleCellExperiment with UMAP and LSI reducedDims
#' @param colorColumn   the colData column to color plots ("mitoCluster")
#'
#' @details presumes NMF and UMAP. you will almost certainly want to tweak this.
#'
#' @import iSEE
#'
#' @export
#'
iSEEarchR <- function(x, colorColumn = "Clusters") { 

  # reason for this will become obvious 
  stopifnot(c("UMAP","LSI") %in% reducedDimNames(x))

  # if reducedDim(x) columns aren't already in colData(x), add them
  for (rdim in reducedDimNames(x)) x <- .reducedDimsAsColData(x, rdim=rdim)

  # launch
  iSEE(x,
     initial = 
       list(UMAP = new("ReducedDimensionPlot",
                       FontSize = 1.5,
                       PointSize = 4,
                       VisualBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       Type = "UMAP"
                       ),
            FACTORS = new("ColumnDataPlot",
                       YAxis = "TSSEnrichment",
                       XAxis = "Column data",
                       XAxisColumnData = "Sample",
                       FontSize = 1.5,
                       PointSize = 4,
                       DataBoxOpen = TRUE,
                       SelectionBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       ),
            COLDAT = new("ColumnDataTable",
                       SelectionBoxOpen = TRUE,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       )
            )
     )
}


# helper fn
.reducedDimsAsColData <- function(x, rdim="NMF") {

  stopifnot(rdim %in% reducedDimNames(x))
  cdnames <- names(colData(x))
  rdimnames <- colnames(reducedDim(x, rdim))
  if (!all(toupper(rdimnames) %in% cdnames)) {
    for (i in rdimnames) { 
      message("Adding ", i, " as colData(SCE)$", toupper(i))
      colData(x)[, toupper(i)] <- reducedDim(x, rdim)[, i]
    }
  }
  return(x) 

}
