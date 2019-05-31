#' Statistics about training data.
#'
#' @param listPre Vector of negative training samples.
#' @param listPost Vector of positive training samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @return Dataframe containing number of samples, colonotypes, and unique seqs in training data.
#' @export
trnStats <- function(listPre, listPost, field) {
  fs <- strsplit(field, ' ')[[1]]
  
  pre <- rbindlist(lapply(listPre, function(x) {
    dat <- fread(x, select = c(fs, "count (templates)", "copy"))
  }))
  
  DTpre <- as.data.table(pre)
  pre <- DTpre[, lapply(.SD, sum), by = fs, .SDcols = "count (templates)"]
  pref <- c(length(listPre), length(pre$`count (templates)`), sum(pre$`count (templates)`))
  
  post <- rbindlist(lapply(listPost, function(x) {
    dat <- fread(x, select = c(fs, "count (templates)"))
  }))
  
  DTpost <- as.data.table(post)
  post <- DTpost[, lapply(.SD, sum), by = fs, .SDcols = "count (templates)"]
  postf <- c(length(listPost), length(post$`count (templates)`), sum(post$`count (templates)`))
  
  both <- rbind(pref, postf)
  colnames(both) <- c("# Samples", "# Clonotypes", "# Unique Sequences")
  rownames(both) <- c("Negative", "Positive")
  
  return(both)
}
