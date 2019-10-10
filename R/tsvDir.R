#' List TSV files in a directory.
#'
#' @param dir String of directory.
#' @return List of *.tsv files in the directory.
#' @export
#'  FIELD <- "vGeneName aminoAcid jGeneName"
#' listPos <- tsvDir("path/to/positve/samples/")
#' listNeg <- tsvDir("path/to/negative/samples/")
#' 
#' naive <- readTrn(listNeg, FIELD, "naive")
#' vaccs <- readTrn(listPos, FIELD, "vacc")  
tsvDir <- function (dir) {
  list <-
    list.files(
      path = dir,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "*.tsv"
    )
}
