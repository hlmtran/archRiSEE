library(archRiSEE)

# test project
library(ArchR)
# Checking Inputs...
# Detected less than 500 Cells.
# `filterBias` disabled.
# `outlierQuantiles` disabled
# `sampleCellsPre` disabled
# `testBias` in `addClusters` disabled
#
proj <- addTileMatrix(proj, binarize=FALSE, force=TRUE)
SEsmall <- ArchR::getMatrixFromProject(proj, "TileMatrix")
metadata(SEsmall)$LSI <- proj@reducedDims$IterativeLSI # just a List 

# test whether by naming the features we can robustify .projectLSI
rownames(SEsmall) 
# NULL 
metadata(SEsmall)$tileSize <- archRtileSize(proj)
rowData(SEsmall)$end <- rowData(SEsmall)$start + metadata(SEsmall)$tileSize
rowData(SEsmall)$start <- rowData(SEsmall)$start + 1
rowRanges(SEsmall) <- as(rowData(SEsmall), "GRanges") 
rownames(SEsmall) <- as.character(rowRanges(SEsmall))
head(rownames(SEsmall))



LSIFeats <- archRiterLSI(proj)$LSIFeatures
# if (is(LSIFeats, "GRanges")) {
#   message("Flagging features which were used for iterative LSI...")
  ol <- findOverlaps(rowRanges(SEsmall), LSIFeats)
  rowRanges(SEsmall)$usedForLSI <- FALSE
  rowRanges(SEsmall)$usedForLSI[queryHits(ol)] <- TRUE
# }

# test projection and LSI mapping to archRiSEE
res <- try(addLSI(SEsmall,useMatrix="TileMatrix"))
if (!inherits(res, "try-error")) {
  subsample <- sample(colnames(SEsmall), 100)
  toProject <- assay(SEsmall)[, subsample]
  LSI <- metadata(SEsmall)$LSI
  projectedMatSVD <- archRiSEE::projectLSI(toProject, LSI)
  testError <- projectedMatSVD - LSI$matSVD[subsample, ] 
}

# test stacking with LSI mapping 


# test projecting on stacked data
