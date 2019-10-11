#' List TSV files in a directory.
#'
#' @param dir String of directory.
#' @return List of *.tsv files in the directory.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' listPos <- tsvDir("iCAT/extdata/Post/")
#' listNeg <- tsvDir("iCAT/extdata/Pre/")
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
