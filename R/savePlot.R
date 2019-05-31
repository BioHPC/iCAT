#' Save plot to pdf.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @export
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
