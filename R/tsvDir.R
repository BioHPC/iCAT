#' List TSV files in a directory.
#'
#' @param dir String of directory.
#' @return List of *.tsv files in the directory.
#' @export
tsvDir <- function (dir) {
  list <-
    list.files(
      path = dir,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "*.tsv"
    )
}
