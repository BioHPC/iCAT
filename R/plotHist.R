#' Plot histogram of clonotype frequencies in negative and positive samples.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @export
#' @example
#' comb <- train(...)
#' plotHist(comb)
plotHist <- function(comb) {
  comb_ <- list(n=comb$n, v=comb$v)
  navvac_hist <-
    ggplot2::ggplot(melt(comb_), ggplot2::aes(x = value, fill = L1)) +
    ggplot2::geom_histogram(color = 1) +
    ggplot2::labs(
      y = "Frequency\n",
      x = "\nPercentages",
      fill = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values=c("royalblue1", "tomato1"), labels = c("Negative", "Positive")) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size=18), panel.grid.major.y = ggplot2::element_line(color = 1, linetype = 'dashed')) +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  navvac_hist # this line outputs the plot if using a console/RStudio
}
