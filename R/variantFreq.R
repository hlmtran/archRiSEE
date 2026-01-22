#' tabulate (usually) offtarget vs. ontarget variants
#'
#' @param SE        SingleCellExperiment or a SummarizedExperiment (see Details)
#' @param grouping  names of factors in colData(SE) to group by (or NULL) 
#' @param loci      rownames to include (or NULL) 
#' @param altexp    'mtVariants', by default; ignored unless SE is an SCE
#' @param asy       'allele_frequency', by default; can also be NULL (assay 1) 
#' @param cutoff    0.05 (five percent), by default; if 0, just return AFs 
#' 
#' @return          a vector or table of cells passing the cutoff at each locus
#'
#' @details         this function is not very clever, it just sets up HET calls
#'                  for tabulating by (eg) cellLine and edit status
#'
#' @examples
#'
#' \dontrun{
#'
#'   # this takes time, even though it IS subsetted
#'   system.time(ewsSub <- readRDS("ewsSub.rds"))
#'   #    user  system elapsed 
#'   # 200.518   5.787 206.487 
#'   
#'   grouping <- c("cellLine", "edited")
#'   loci <- c("11922G>A", "11925C>T", "11946C>T")
#'   res <- variantFreq(ewsSub, grouping=grouping, loci=loci)
#'   names(res) 
#'   # [1] "11922G>A" "11925C>T" "11946C>T" "cellLine" "edited"  
#'   
#'   (names(res)[1:3] <- paste0("chrM_", sub('>', '_', names(res)[1:3])))
#'   # [1] "chrM_11922G_A" "chrM_11925C_T" "chrM_11946C_T"
#'   
#'   with(res, table(chrM_11922G_A, edited, cellLine))  # on-target
#'   with(res, table(chrM_11925C_T, edited, cellLine))  # off-target 
#' 
#' \}
#'
#'
#' @seealso archRtoSCE 
#'
#' @import SingleCellExperiment
#'
#' @export 
#'
variantFreq   <- function(SE, 
                          grouping=NULL, 
                          loci=NULL,
                          altexp="mtVariants", 
                          assay="allele_frequency", 
                          cutoff=0.05) {

  if (is(SE, "SingleCellExperiment")) {
    stopifnot(altexp %in% c(mainExpName(SE), altExpNames(SE)))
    if (altexp %in% altExpNames(SE)) SE <- altExp(SE, altexp)
  }
  if (!is.null(loci)) SE <- SE[loci, ] 
  if (!is.null(grouping)) stopifnot(all(grouping %in% names(colData(SE))))

  res <- assay(SE, assay) 
  if (cutoff > 0) res <- (res > cutoff)
  res <- as.data.frame(data.matrix(t(res))) # often a dgCMatrix or lgCMatrix
  if (!is.null(grouping)) for (g in grouping) res[, g] <- colData(SE)[, g]
  return(res)

}
