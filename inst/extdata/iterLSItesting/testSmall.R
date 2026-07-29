library(archRiSEE)

# test project
library(ArchR)
proj <- getTestProject()
proj <- addTileMatrix(proj,,tileSize=20000,binarize=TRUE,force=TRUE)
proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)
#
# Checking Inputs...
# Detected less than 500 Cells.
# `filterBias` disabled.
# `outlierQuantiles` disabled
# `sampleCellsPre` disabled
# `testBias` in `addClusters` disabled
#
# proj <- addTileMatrix(proj) 
SEsmall <- ArchR::getMatrixFromProject(proj, "TileMatrix", binarize=TRUE)
metadata(SEsmall)$LSI <- proj@reducedDims$IterativeLSI # just a List 
LSIfeats <- match(metadata(SEsmall)$LSI$idx, mcols(SEsmall)$idx)
mcols(SEsmall)$usedForLSI <- (mcols(SEsmall)$idx %in% metadata(SEsmall)$LSI$idx)
all(mcols(SEsmall)$usedForLSI[LSIfeats])

# test whether by naming the features we can robustify .projectLSI
rownames(SEsmall) 
# NULL 
metadata(SEsmall)$tileSize <- archRtileSize(proj)
rowData(SEsmall)$end <- rowData(SEsmall)$start + metadata(SEsmall)$tileSize
rowData(SEsmall)$start <- rowData(SEsmall)$start + 1
rowRanges(SEsmall) <- as(rowData(SEsmall), "GRanges") 
rownames(SEsmall) <- as.character(rowRanges(SEsmall))
head(rownames(SEsmall))

# test projection and LSI mapping to archRiSEE
res <- try(addLSI(SEsmall))
# 'TFIDF' not in names(assays(<RangedSummarizedExperiment>))
SEsmall <- addTfIdf(SEsmall)

# try again
res <- try(addLSI(SEsmall))

if (!inherits(res, "try-error")) {
  subsample <- sample(colnames(SEsmall), 100)
  toProject <- assay(SEsmall)[, subsample]
  LSI <- metadata(SEsmall)$LSI
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


projected_test <- archRiSEE::projectLSI(mat_iter1_20k, LSI_iter1_20k)
identical(projected_test,matSVD_iter1_20k)
