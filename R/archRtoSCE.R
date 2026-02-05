#' turn an ArchR project into a SingleCellExperiment 
#'
#' @param proj          an ArchRproject 
#' @param how           where to get counts from ('tiles','feats','LSI')
#' @param feats         GRanges of features to extract (ideally 500bp; NULL)
#' @param addNMF        add an NMF decomposition of logCounts? (FALSE)
#' @param colDat        add LSI and/or NMF scores to colData(SCE)? (FALSE)
#' @param LSIdim        name of IterativeLSI reducedDim ("IterativeLSI") 
#' @param tileSize      tileSize (default is whatever is in use)
#' @param keepbinary    keep tiles binarized if they are already? (FALSE) 
#' @param ...           additional arguments for ArchR::getMatrixFromProject()
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
#'          This function is way too big nowadays and should be refactored. 
#'
#' @import RcppML
#' @import scuttle
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, how=c("tiles","feats","LSI"), feats=NULL, addNMF=FALSE, colDat=FALSE, LSIdim="IterativeLSI", tileSize=NULL, keepbinary=TRUE, ...) {

  if (!require(ArchR)) stop("This function won't work without an ArchR install")
  how <- match.arg(how) 
  if (!is.null(feats)) how <- "feats"
  LSI <- archRiterLSI(proj)
  tile <- archRtileSize(proj)
  if (is.null(tileSize)) tileSize <- tile
  g <- archRgenome(proj)
  # "binarize beforehand"
  bb <- FALSE

  # key: which ones?
  if (how == "feats") { 

    # grab user-provided features 
    message("Adding 'feats' matrix to project (usually fast)...")
    proj <- addFeatureMatrix(proj, features=feats, matrixName="feats",
                             binarize=bb, force=TRUE)
    message("Converting to a SingleCellExperiment...") 
    SCE <- as(getMatrixFromProject(proj, "feats"), "SingleCellExperiment") 
    rownames(SCE) <- as.character(rowRanges(SCE))
    assayNames(SCE) <- "counts"

  } else if (how == "LSI") {
    
    if (tile == 500) {
      warning("IterativeLSI run with tile size of ", tile, ", be sure it's OK!")
    }

    # grab LSI-defined features
    message("Adding 'LSI' matrix to project (usually fast)...")
    proj <- addFeatureMatrix(proj, features=LSI$LSIFeatures, 
                             matrixName="LSI", binarize=bb,
                             force=TRUE)
    SCE <- as(getMatrixFromProject(proj, "LSI"), "SingleCellExperiment") 
    rownames(SCE) <- as.character(rowRanges(SCE))
    rowData(SCE)$usedForLSI <- TRUE
    assayNames(SCE) <- "counts"

  } else { 

    if (tileSize == 500) {
      warning("tileSize set to 500 (millions of tiles). Be sure you want this!")
    }
        
    # grab everything (after possibly warning the user about the above) 
    if (tile == tileSize & "TileMatrix" %in% getAvailableMatrices(proj)) {
      b <- unique(sapply(getArrowFiles(proj), h5read, "TileMatrix/Info/Class"))
      bb <- (b == "Sparse.Binary.Matrix")
      if (keepbinary) {
        message("Found existing TileMatrix with desired size, using that...")
        message("(Existing matrix is of type ", b, ", if that matters!)")
      } else { 
        if (bb) {
          message("Existing TileMatrix is binarized, replacing with counts...")
          proj <- addTileMatrix(proj, force=TRUE, binarize=FALSE, tileSize=tile,
                                )
        }
      }
    } else { 
      message("Adding TileMatrix to ArchR project (this may take a while)...")
      bb <- FALSE
      proj <- addTileMatrix(proj,force=TRUE,binarize=bb,tileSize=tileSize, ...)
    }
    
    SCE <- as(getMatrixFromProject(proj, "TileMatrix", binarize=bb, ...), 
              "SingleCellExperiment")
    assayNames(SCE) <- "counts"
    rowData(SCE)$end <- rowData(SCE)$start + (tile - 1)
    rowRanges(SCE) <- as(rowData(SCE), "GRanges")
    rownames(SCE) <- as(rowRanges(SCE), "character")
    SCE <- sort(sortSeqlevels(SCE))
    genome(SCE) <- g
    
    # flag features used for LSI, if found
    if (!"usedForLSI" %in% names(rowData(SCE))) {
      LSIFeats <- archRiterLSI(proj)$LSIFeatures
      if (is(LSIFeats, "GRanges")) { 
        message("Flagging features which were used for iterative LSI...")
        ol <- findOverlaps(rowRanges(SCE), LSIFeats)
        rowRanges(SCE)$usedForLSI <- FALSE
        rowRanges(SCE)$usedForLSI[queryHits(ol)] <- TRUE
      }
    }

  }

  # since it's feasible to stack experiments (e.g. DEM + H3K4me + H3K27me)
  mainExpName(SCE) <- "FragmentCounts"
  rowRanges(SCE)$idx <- NULL # irrelevant to us and gets in the way
  rowData(SCE)$assay <- "FragmentCounts" # for binding to any other altExps
  message("You may want to update mcols(SCE)$assay to be more specific.")
  message("(For example, 'DEM' or 'H3K27me3' or 'H3K4me3' or 'ATAC'...)")
  if (max(assay(SCE)) > 1) { 
    message("Log-normalizing fragment counts...") 
    SCE <- logNormCounts(SCE, assay.type="counts")
  } else {
    message("Binarized data, log-normalization makes no sense...")
  }
  colData(SCE) <- getCellColData(proj)
  
  message("Copying UMAP to reducedDim(SCE, 'UMAP')...")
  reducedDim(SCE, "UMAP") <- getEmbedding(proj, "UMAP")[colnames(SCE), ]
  names(reducedDim(SCE, "UMAP")) <- c("UMAP1", "UMAP2")
  message("Copying UMAP parameters to metadata(SCE)$UMAP...")
  metadata(SCE)$UMAP <- getEmbedding(proj, "UMAP", returnDF=FALSE)$params
  message("Copying LSI scores to reducedDim(SCE, 'LSI')...")
  reducedDim(SCE, "LSI") <- getReducedDims(proj)[colnames(SCE), ]

  if (colDat) { 
    message("Copying LSI dimensions to colData(SCE) for iSEE visualization...")
    for (i in colnames(reducedDim(SCE, "LSI"))) {
      message("Adding ", i, " as colData(SCE)$", i)
      colData(SCE)[, i] <- reducedDim(SCE, "LSI")[colnames(SCE), i]
    }
  }

  # could also add others if it makes sense here 
  for (mat in c("GeneScoreMatrix", "PeakMatrix")) { 
    if (mat %in% getAvailableMatrices(proj)) {
      message("Attempting to add ", mat, " to altExp(SCE, '", mat, "')...")
      if (mat == "GeneScoreMatrix") { 
        altExp(SCE, mat) <- archRgeneScoreMatrix(proj, genome=g)
      } else if (mat == "PeakMatrix" & !is.null(getPeakSet(proj))) { 
        altExp(SCE, mat) <- archRpeakMatrix(proj, genome=g)
      } else { 
        altExp(SCE, mat) <- getMatrixFromProject(proj, mat)
      }
    }
  }

  if (addNMF) {
    if (max(assay(SCE)) > 1) {
      SCE <- addNMF(SCE, k=k, colDat=colDat, rowDat=TRUE)
    } else { 
      warning("Your fragment counts are binarized, NMF will not work properly")
    }
  }

  message("Copying metadata...")
  md <- archRmetadata(proj)
  for (i in names(md)) metadata(SCE)[[i]] <- md[[i]]
  message("Done.")
  return(SCE)

}
