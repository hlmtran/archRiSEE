#' turn an ArchR project into a SingleCellExperiment 
#'
#' @param proj          an ArchRproject 
#' @param how           where to get counts from ('tiles','feats','LSI')
#' @param feats         GRanges of features to extract (ideally 500bp; NULL)
#' @param addNMF        add an NMF decomposition of logCounts? (FALSE)
#' @param colDat        add LSI and/or NMF scores to colData(SCE)? (FALSE)
#' @param LSIdim        name of IterativeLSI reducedDim ("IterativeLSI") 
#' @param tileSize      tileSize (default is whatever is in use)
#' @param ...           additional arguments
#'
#' @return              a SingleCellExperiment 
#'
#' @details Presumes LSI and UMAP. you will almost certainly want to tweak this.
#'          If !is.null(feats), (GRanges, please!) we assume that how == "feats"
#'          If you are working with single-cell histone or transcription factor
#'          data, you will most likely want to use larger (5kb-50kb) tiles (see 
#'          Janssens et al, Nature Protocols, 2024 for sciCUT&Tag benchmarks).
#'          This function has changed recently (2026) and the assumption is
#'          that most users will want a subset of features. Pay attention to 
#'          any warning messages -- they can save you a lot of time and RAM!
#'
#' @import RcppML
#' @import scuttle
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, how=c("tiles","feats","LSI"), feats=NULL, addNMF=FALSE, colDat=FALSE, LSIdim="IterativeLSI", tileSize=NULL, ...) { 

  if (!require(ArchR)) stop("This function won't work without an ArchR install")
  how <- match.arg(how) 
  if (!is.null(feats)) how <- "feats"
  LSI <- archRiterLSI(proj)
  tile <- archRtileSize(proj)
  if (is.null(tileSize)) tileSize <- tile

  # key: which ones?
  if (how == "feats") { 

    # {{{ grab user-provided features 
    message("Adding 'feats' matrix to project (usually fast)...")
    proj <- addFeatureMatrix(proj, features=feats, matrixName="feats",
                             binarize=FALSE, force=TRUE)
    message("Converting to a SingleCellExperiment...") 
    SCE <- as(getMatrixFromProject(proj, "feats"), "SingleCellExperiment") 
    rownames(SCE) <- as.character(rowRanges(SCE))
    assayNames(SCE) <- "counts"
    # }}}

  } else if (how == "LSI") {
    
    if (tile == 500) {
      warning("IterativeLSI run with tile size of ", tile, ", be sure it's OK!")
    }

    # {{{ grab LSI-defined features
    message("Adding 'LSI' matrix to project (usually fast)...")
    proj <- addFeatureMatrix(proj, features=LSI$LSIFeatures, 
                             matrixName="LSI", binarize=FALSE,
                             force=TRUE)
    SCE <- as(getMatrixFromProject(proj, "LSI"), "SingleCellExperiment") 
    rownames(SCE) <- as.character(rowRanges(SCE))
    assayNames(SCE) <- "counts"
    # }}}

  } else { 

    if (tileSize == 500) {
      warning("tileSize set to 500 (millions of tiles). Be sure you want this!")
    }
        
    # {{{ grab everything (after possibly warning the user about the above) 
    if (tile != tileSize | !"TileMatrix" %in% getAvailableMatrices(proj)) {
      message("Adding TileMatrix to ArchR project (this may take a while)...")
      proj <- addTileMatrix(proj, force=TRUE, binarize=FALSE, tileSize=tileSize,
                            ...)
    } else { 
      message("Found existing TileMatrix with desired size, using that...")
    }
    SCE <- as(getMatrixFromProject(proj, "TileMatrix"), "SingleCellExperiment") 
    assayNames(SCE) <- "counts"
    rowData(SCE)$end <- rowData(SCE)$start + (tile - 1)
    rowRanges(SCE) <- as(rowData(SCE), "GRanges")
    
    # flag features used for LSI, if found
    LSIFeats <- archRiterLSI(proj)$LSIFeatures
    if (is(LSIFeats, "GRanges")) {
      message("Flagging features which were used for iterative LSI...")
      rowRanges(SCE)$usedForLSI <- rowRanges(SCE) %in% LSIFeats
    }

    genome(SCE) <- archRgenome(proj)
    SCE <- sort(sortSeqlevels(SCE))
    rownames(SCE) <- as.character(rowRanges(SCE))
    # }}}

  }

  # since it's feasible to stack experiments (e.g. DEM + H3K4me + H3K27me)
  rowData(SCE)$assay <- "FragmentCounts" # for binding to any other altExps
  message("Log-normalizing fragment counts...") 
  SCE <- logNormCounts(SCE, assay.type="counts")
  mainExpName(SCE) <- "FragmentCounts"
  message("You may want to update mcols(SCE)$assay to be more specific.")
  message("(For example, 'DEM' or 'H3K27me3' or 'H3K4me3' or 'ATAC'...)")

  colData(SCE) <- proj@cellColData
  
  message("Copying UMAP to reducedDim(SCE, 'UMAP')...")
  reducedDim(SCE, "UMAP") <- proj@embeddings$UMAP$df
  names(reducedDim(SCE, "UMAP")) <- c("UMAP1", "UMAP2")
  message("Copying UMAP parameters to metadata(SCE)$UMAP$params...")
  metadata(SCE)$UMAP_params <- list(params=proj@embeddings$UMAP$params)

  message("Copying LSI scores to reducedDim(SCE, 'LSI')...")
  reducedDim(SCE, "LSI") <- .LSI(proj)$matSVD
  colnames(reducedDim(SCE,"LSI")) <- 
    paste0("LSI", seq_len(ncol(reducedDim(SCE, "LSI"))))
  
  if (colDat) { 
    message("Copying LSI dimensions to colData(SCE) for iSEE visualization...")
    for (i in colnames(reducedDim(SCE, "LSI"))) {
      message("Adding ", i, " as colData(SCE)$", i)
      colData(SCE)[, i] <- reducedDim(SCE, "LSI")[, i]
    }
  }

  message("Copying metadata...")
  if (addNMF) SCE <- addNMF(SCE, k=k, colDat=colDat, rowDat=TRUE)
  md <- archRmetadata(proj)
  for (i in names(md)) metadata(SCE)[[i]] <- md[[i]]
  message("Done.")
  return(SCE)

}
