#' Read input training data.
#'
#' @param list Vector of training samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @param posOrNeg String indicating whether the data is negative or positive.
#' @return Dataframe containing unique sequences and their frequencies in the samples.
#' @export
readTrn <- function(list, field, posOrNeg) {
  if (length(list) == 0) {
    return(NULL)
  }

  fs <- strsplit(field, ' ')[[1]]
  #combine files and count repeats

  final <- lapply(list, processFiles, field=field)
  final <- data.table::rbindlist(final)
  x <- rep("", length(final[[1]]))
  for (i in 1:length(fs))
    x <- paste(x, final[[i]])
  x <- substring(x, 2)
  final <-
    data.table::data.table(x)
  colnames(final) <- "names"
  final <- final[, .N, by = names(final)]
  secCol <- paste(posOrNeg, "amounts", sep='')
  colnames(final) <- c("names", secCol)
  return(final)
}
