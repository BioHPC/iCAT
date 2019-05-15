

#' Launch interactive interface
#'
#' @examples
#' if(interactive()) {iCATinteractive()}
#' @export
iCATinteractive <- function() {

  dir <- system.file("app", package = "iCAT")
  shiny::runApp(dir, quiet=TRUE, launch.browser=TRUE)
}


