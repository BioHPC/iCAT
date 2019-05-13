

#' Launch interactive interface
#'
#' @examples
#' if(interactive()) {iCATinteractive()}
#' @export
iCATinteractive <- function() {

  print(getwd())
  shiny::runApp('app', quiet=TRUE, launch.browser=TRUE)
}


