library(archRiSEE)

# test project
library(ArchR)
proj <- getTestProject()
proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)
#
# Checking Inputs...
# Detected less than 500 Cells.
# `filterBias` disabled.
# `outlierQuantiles` disabled
# `sampleCellsPre` disabled
# `testBias` in `addClusters` disabled
#
proj <- addTileMatrix(proj) 
SEsmall <- ArchR::getMatrixFromProject(proj, "TileMatrix", binarize=TRUE)

# convert to SCE so that we can use ReducedDims
library(SingleCellExperiment) 
SCEsmall <- as(SEsmall, "SingleCellExperiment") 
metadata(SCEsmall)$LSI <- proj@reducedDims$IterativeLSI # just a List 

# test whether by naming the features we can robustify .projectLSI
rownames(SCEsmall) 
# NULL 
metadata(SCEsmall)$tileSize <- archRtileSize(proj)
stopifnot(metadata(SCEsmall)$tileSize == metadata(SCEsmall)$LSI$tileSize)
rowData(SCEsmall)$end <- rowData(SCEsmall)$start + metadata(SCEsmall)$tileSize
rowData(SCEsmall)$start <- rowData(SCEsmall)$start + 1
rowRanges(SCEsmall) <- as(rowData(SCEsmall), "GRanges") 
rownames(SCEsmall) <- as.character(rowRanges(SCEsmall))
head(rownames(SCEsmall))

# now do the same with LSI$LSIFeatures 
LSIfeats <- metadata(SCEsmall)$LSI$LSIFeatures
LSIfeats$end <- LSIfeats$start + metadata(SCEsmall)$tileSize
LSIfeats$start <- LSIfeats$start + 1 
LSIgr <- as(LSIfeats, "GRanges") 
names(LSIgr) <- as(LSIgr, "character") 

# are these the features we are looking for?
identical(unname(rowSums(assay(subsetByOverlaps(SCEsmall, LSIgr)))), 
          LSIgr$rowSums)
# [1] TRUE

ol <- findOverlaps(SCEsmall, LSIgr)
rowData(SCEsmall)$usedForLSI <- FALSE
rowData(SCEsmall)$usedForLSI[queryHits(ol)] <- TRUE 
table(rowData(SCEsmall)$usedForLSI)
# 
# FALSE  TRUE 
# 30593  1000 

# need to mask 0-sum rows 
rowData(SCEsmall)$rowSm <- rowSums(assay(SCEsmall))

# test projection and LSI mapping to archRiSEE
# 'TFIDF' not in names(assays(<RangedSummarizedExperiment>))
SCEsmall <- addTfIdf(SCEsmall, prune=1)

# try again
res <- try(addLSI(SCEsmall))

if (!inherits(res, "try-error")) {
  subsample <- sample(colnames(SCEsmall), 100)
  toProject <- assay(SCEsmall)[, subsample]
  LSI <- metadata(SCEsmall)$LSI
  projectedMatSVD <- archRiSEE::projectLSI(toProject, LSI)
  # Subsetting TF-IDF matrix...
  # Running SVD...
  # Warning in irlba::irlba(mat, nDimensions + 5, nDimensions + 5) :
  #   convergence criterion below machine epsilon
  # Warning in irlba::irlba(mat, nDimensions + 5, nDimensions + 5) :
  #   did not converge--results might be invalid!; try increasing work or maxit
  # Scaling matSVD...
  # Checking for depth-correlated columns...
  # Kept 0% of columns.
  # Error in matSVD[, toKeep][, seq_len(nDimensions)] : 
  #   subscript out of bounds
  testError <- projectedMatSVD - LSI$matSVD[subsample, ] 
}

# test stacking with LSI mapping 


# test projecting on stacked data
