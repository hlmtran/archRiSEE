#' compute sparsity for a Matrix, SummarizedExperiment, or SingleCellExperiment
#' 
#' @param object  the object whose sparsity we would like to evaluate
#' 
#' @return        1 - (length(object@x)/prod(dim(object))), i.e. (1 - density)
#' 
#' @details       just adding a method for Matrix objects to RcppML's generic.
#'                if a SummarizedExperiment is passed, the first assay is used.
#'
#' @import        RcppML
#' @import        Matrix
#' @import        SummarizedExperiment
#'
#' @rdname        sparsity 
#'
#' @export
#'
setMethod("sparsity", "Matrix", 
          function(object) 1 - (length(object@x) / prod(dim(object))))
setMethod("sparsity", "SummarizedExperiment", 
          function(object) 1 - (length(assay(object)@x) / prod(dim(object))))
