#' Dataframe of TARSs library.
#'
#' @param comb composed list returned from train() function
#' @return Dataframe containing the significant target-associated receptor sequences.
#' @export
getLib <- function(comb) {
  lib <- comb$l
  l <- lib[, unique(colnames(lib)), with=FALSE]
  colnames(l) <- c("Sequence",	"Positive Present",	"Negative Present",	"Total",	"Negative Absent",	"Positive Absent",	"P-Value")
  #l <- l[, c(1, 2, 6, 3, 5, 7), with=FALSE]
  l <- l[, c(1, 7), with=FALSE]
  return(l)
}
