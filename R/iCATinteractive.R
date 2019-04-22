

#' Launch interactive interface
#'
#'
#' @examples
#' if(interactive()) {iCATinteractive()}
#' @export
iCATinteractive <- function() {

print(getwd())
runApp('iCAT/R/app', quiet=TRUE, launch.browser=TRUE)

}


