#' Save plot to pdf.
#'
#' @param file name to give plot file
#' @param plotIn plot object
#' @return ggplot2::ggsave return
#' @export
#' @examples 
#' savePlot("plot.pdf", plotObj)
savePlot <- function(file, plotIn) {
  ggplot2::ggsave(
    device = "pdf",
    plot = plotIn,
    width = 5,
    height = 5,
    units = "in",
    dpi = 200
  )
}
