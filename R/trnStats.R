#' Statistics about training data.
#'
#' @param listPre Vector of negative training samples.
#' @param listPost Vector of positive training samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @param count String containing the column name for colontype counts.
#' @param copyrange Integer Vector of the min and max copy of a sequence, within a sample, to be considered.  
#' @return Dataframe containing number of samples, colonotypes, and unique seqs in training data.
#' @export
#' @examples 
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' COUNT <- "copy"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' COPY_RANGE <- "1 99"
#' 
#' listPos <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' 
#' trnStats(listPos, listNeg, FIELD, COUNT, COPY_RANGE)
trnStats <- function(listPre, listPost, field, count, copyrange) {
  fs <- strsplit(field, ' ')[[1]]
  copyrange <- as.integer(strsplit(copyrange, ' ')[[1]])
  
  pre <- rbindlist(lapply(listPre, function(x) {
    dat <- fread(x, select = c(fs, count))
  }))
  
  DTpre <- as.data.table(pre)
  DTpre <- DTpre[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
  pre <- DTpre[, lapply(.SD, sum), by = fs, .SDcols = count]
  pref <- c(length(listPre), length(pre[[count]]), sum(pre[[count]]))
  
  post <- rbindlist(lapply(listPost, function(x) {
    dat <- fread(x, select = c(fs, count))
  }))
  
  DTpost <- as.data.table(post)
  DTpost <- DTpost[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
  post <- DTpost[, lapply(.SD, sum), by = fs, .SDcols = count]
  postf <- c(length(listPost), length(post[[count]]), sum(post[[count]]))
  
  both <- rbind(pref, postf)
  colnames(both) <- c("# Samples", "# Unique Sequences", "# Clonotypes")
  rownames(both) <- c("Negative", "Positive")
  
  return(both)
}
