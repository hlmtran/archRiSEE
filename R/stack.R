#' stack (rbind) two SCE-like objects for the same cells to integrate
#' 
#' @param lst   a list of SCE-like objects with identical colnames()
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
#'              Currently, NMF and UMAP on the combined features is disabled
#'              because the default binary output makes no sense with NMF.
#'              If you plan to use raw or lognormalized FragmentCounts again, 
#'              make sure that the object you want to use these on is lst[[1]],
#'              else you'll have to 
#'
#' @import      SingleCellExperiment
#'
#' @export
#'
stack <- function(lst, ...) {

  if (is.null(names(lst))) {
    stop("Your list of experiments must have names to disambiguate them.")
  }
  if (!Reduce(identical, lapply(lst, colnames))) {
    stop("Error: you must ensure identical column names across assays.")
  }
  for (i in names(lst)) { 
    names(colData(lst[[i]])) <- paste0(i, ".", names(colData(lst[[i]])))
    names(rowData(lst[[i]])) <- paste0(i, ".", names(rowData(lst[[i]])))
    mainExpName(lst[[i]]) <- paste0(i, ".", mainExpName(lst[[i]]))
    assayNames(lst[[i]]) <- paste0(i, ".", assayNames(lst[[i]]))
    altExpNames(lst[[i]]) <- paste0(i, ".", altExpNames(lst[[i]]))
    reducedDimNames(lst[[i]]) <- paste0(i, ".", reducedDimNames(lst[[i]]))
  }
  
  # yes it's ugly, no I don't care
  stacked <- lst[[1]]
  lst <- lst[-1]
  for (j in names(lst)) {
    message("Adding data from ", j, "...") 
    colData(stacked)[, names(colData(lst[[j]]))] <- colData(lst[[j]])
    for (k in reducedDimNames(lst[[j]])) {
      if (k %in% reducedDimNames(stacked)) {
        warning("reducedDim(stacked, '", k, "') exists and will be replaced")
      }
      message("Adding ", k, " to reducedDims(stacked)...")
      reducedDim(stacked, k) <- reducedDim(lst[[j]], k)
    }
    for (l in altExpNames(lst[[j]])) {
      if (l %in% c(altExpNames(stacked), mainExpName(stacked))) {
        warning("altExp(stacked, '", l, "') exists and will be replaced")
      }
      message("Adding ", l, " to altExps(stacked)...")
      altExp(stacked, l) <- altExp(lst[[j]], l)
    }
    m <- mainExpName(lst[[j]])
    if (m %in% c(altExpNames(stacked), mainExpName(stacked))) {
      warning("altExp(stacked, '", m, "') exists and will be replaced")
    }
    altExp(stacked, m) <- as(lst[[j]], "SummarizedExperiment")
    
    # flag features of interest in each constituent dataset
    if (identical(rownames(stacked), rownames(lst[[j]]))) {
      n <- grep(paste0(j, ".usedForLSI"), names(rowData(lst[[j]])), val=TRUE)
      rowData(stacked)[, n] <- rowData(lst[[j]])[, n]
    }
  }

  return(stacked) 

}
