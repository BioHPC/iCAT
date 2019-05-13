#' Internal function for transforming input files.
#'
#' @param x Input file.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @export
#' @examples
#' lapply(list, processFiles)
processFiles <- function(x, field) {
  fs <- strsplit(field, ' ')[[1]]
  # get sequenceStatus column to use as filter, then remove it
  #dat <- readLines(x)
  #dat <- dat[grep('^[A-Za-z]', dat)]
  #dat <- read.table(textConnection(dat), sep="\t", header=T)
  #dat <- select(dat, fs, "sequenceStatus")
  dat <- data.table::fread(x, select = c(fs, "sequenceStatus"))
  dat <- dat[dat$sequenceStatus == "In", , drop=F]

  dat <- dat[, 1:length(fs), drop=F]

  dat <- dat[!duplicated(dat), ,drop=F]
  dat <- dat[order(dat[,1:length(fs)]), ,drop=F]
  #print(dat)
}
