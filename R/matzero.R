
#' Title matzero
#'
#' @param mat matrix of spectra
#' @param spks number of the peak selected as the starting point
#' @param x variable to be reset according to the position of the selected peak
#' @param colname name of the selected column
#' @param myby  varialbe selected for subsetting the matrix
#'
#' @description search in the checked matrix by checkmat, the starting of the selected peak for each heating/cooling cycles and zeroes the time (in seconds) 
#'
matzero <- function(mat, spks=1, x=mat$time.seconds.zero, colname="v.check",myby="id_cycle") {
  pks <- NULL
  mat.zero <- mat[mat[, pks == spks]][, c(colname) := (x - min(x)), by = myby][]
  return(mat.zero)
}
