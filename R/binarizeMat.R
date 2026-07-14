#' pretty much what it sounds like: binarize a sparse matrix
#'
#' @param mat the matrix
#' 
#' @return    a binary version of mat
#'
#' @import    Matrix
#'
#' @export
#'
binarizeMat <- function(mat) { 

  if (is(mat, "sparseMatrix")) { 
    mat@x[mat@x > 0] <- 1
  } else {
    ifelse(mat > 0, 1, 0)
  }
  return(mat)

}
