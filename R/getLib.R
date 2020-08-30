#' Dataframe of TARSs library.
#'
#' @param comb composed list returned from train() function
#' @return Dataframe containing the significant target-associated receptor sequences.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' COUNT <- "copy"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' COPY_RANGE <- "1 99"
#' 
#' listPos <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, COUNT, COPY_RANGE, "naive")
#' vaccs <- readTrn(listPos, FIELD, COUNT, COPY_RANGE, "vacc")  
#' 
#' mod <- train(naive, vaccs, listNeg, listPos, FIELD, COUNT, COPY_RANGE, P_CUTOFF, MIN_PUBLIC, NULL)
#' getLib(mod)
getLib <- function(comb) {
  lib <- comb$l
  l <- lib[, unique(colnames(lib)), with=FALSE]
  colnames(l) <- c("Sequence",	"Positive Present",	"Negative Present",	"Total",	"Negative Absent",	"Positive Absent",	"P-Value")
  l <- l[, c(1, 2, 6, 3, 5, 7), with=FALSE]
  #l <- l[, c(1, 7), with=FALSE]
  return(l)
}
