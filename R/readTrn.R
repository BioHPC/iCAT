#' Read input training data.
#'
#' @param list Vector of training samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @param count String containing the column name for colontype counts.
#' @param copyrange Integer of the minimum copy of a sequence, within a sample, to be considered.  
#' @param posOrNeg String indicating whether the data is negative or positive.
#' @return Dataframe containing unique sequences and their frequencies in the samples.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' COUNT <- "count (templates)"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' COPY_RANGE <- "1 99"
#' 
#' listPos <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, COUNT, COPY_RANGE, "naive")
#' vaccs <- readTrn(listPos, FIELD, COUNT, COPY_RANGE, "vacc")
readTrn <- function(list, field, count, copyrange, posOrNeg) {
  if (length(list) == 0) {
    return(NULL)
  }

  fs <- strsplit(field, ' ')[[1]]
  copyrange <- as.integer(strsplit(copyrange, ' ')[[1]])
  #combine files and count repeats

  final <- lapply(list, processFiles, field=field, count=count)
  final <- data.table::rbindlist(final)
  final <- final[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
  final <- final[, names := do.call(paste,.SD), .SDcols=!count]
  final <- final[, c(fs, count) := NULL]
  final <- final[, .N, by = names(final)]
  
  
  # final <- final[,  lapply(.SD, sum) , by = fs, .SDcols = count]
  # final <- final[, names := do.call(paste,.SD), .SDcols=!count]
  # final <- final[, (fs) := NULL]
  #setcolorder(final, c("names", count))

  secCol <- paste(posOrNeg, "amounts", sep='')
  colnames(final) <- c("names", secCol)
  return(final)
}