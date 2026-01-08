#' turn an ArchR project into a SingleCellExperiment 
#'
#' The default tile size is 500bp, just like ArchR's default. 
#'
#' @param proj          an ArchRproject 
#' @param tile          tile size (500)
#' @param addNMF        add an NMF decomposition of logCounts? (FALSE)
#' @param k             rank for the NMF decomposition on logCounts (30) 
#' @param colDat        add LSI and/or NMF scores to colData(SCE)? (FALSE)
#' @param ...           additional arguments to pass to ArchR::addTileMatrix
#'
#' @return              a SingleCellExperiment 
#'
#' @details presumes LSI and UMAP. you will almost certainly want to tweak this.
#'
#' @import RcppML
#' @import scuttle
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, tile=500, addNMF=FALSE, k=30, colDat=FALSE, ...) { 

  if (!require(ArchR)) stop("This function won't work without an ArchR install")
  message("Adding TileMatrix to ArchR project (this may take a while)...")
  proj <- addTileMatrix(proj, force=TRUE, binarize=FALSE, tileSize=tile, ...)
  message("Converting to SingleCellExperiment...")
  SCE <- as(getMatrixFromProject(proj, "TileMatrix"), "SingleCellExperiment") 
  rowData(SCE)$assay<- "FragmentCounts" # for binding to any other altExps
  message("Log-normalizing fragment counts...") 
  SCE <- logNormCounts(SCE, assay.type="TileMatrix")
  mainExpName(SCE) <- "FragmentCounts"
  rowData(SCE)$end <- rowData(SCE)$start + (tile - 1) 
  rowRanges(SCE) <- as(rowData(SCE), "GRanges")
  rownames(SCE) <- as.character(rowRanges(SCE))
  colData(SCE) <- proj@cellColData
  reducedDim(SCE, "UMAP") <- proj@embeddings$UMAP$df
  names(reducedDim(SCE, "UMAP")) <- c("UMAP1", "UMAP2")

  message("Copying LSI scores to reducedDim(SCE, 'LSI')...")
  reducedDim(SCE, "LSI") <- proj@reducedDims$IterativeLSI$matSVD
  LSIdims <- ncol(reducedDim(SCE, "LSI"))
  names(reducedDim(SCE, "LSI")) <- paste0("LSI", seq_len(LSIdims))
  if (colDat) { 
    message("Copying LSI dimensions to colData(SCE) for iSEE visualization...")
    for (i in colnames(reducedDim(SCE, "LSI"))) {
      message("Adding ", i, " as colData(SCE)$", i)
      colData(SCE)[, i] <- reducedDim(SCE, "LSI")[, i]
    }
  }

  if (addNMF) SCE <- addNMF(SCE, k=k, colDat=colDat, rowDat=TRUE)

  message("Done.")
  return(SCE)

}
