#' List TSV files in a directory.
#'
#' @param dir String of directory.
#' @return List of *.tsv files in the directory.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' listPos <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, "naive")
#' vaccs <- readTrn(listPos, FIELD, "vacc")  
tsvDir <- function (dir) {
    list.files(
      path = dir,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "*.tsv"
    )
}
