#' figure out what tile size was used for an ArchR iterative LSI run
#' 
#' @param   proj    the ArchRProject
#'
#' @return  unique(h5read(ArrowFile, "TileMatrix/Info/Params")$tileSize)
#'
#' @details the above is an oversimplification (!), see code for grim details
#'
#' @export
#'
archRtileSize <- function(proj) { 
  
  sampleTileSizes <- 
    sapply(getSampleColData(proj)[, "ArrowFiles"],
           function(x) unique(h5read(x, "TileMatrix/Info/Params/")$tileSize))
  names(sampleTileSizes) <- sapply(names(sampleTileSizes), basename)
  res <- unique(sampleTileSizes)
  if (length(res) > 1) warning("You have inconsistent tile sizes. This is bad.")
  return(res) 

}
