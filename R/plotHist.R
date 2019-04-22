#' Plot histogram of clonotype frequencies in negative and positive samples.
#' 
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @export
#' @example
#' comb <- train(...)
#' plotHist(comb)
plotHist <- function(comb) {
  comb_ <- list(n=comb$n, v=comb$v)
  print(comb)
  print(comb_)
  navvac_hist <-
    ggplot(melt(comb_), aes(x = value, fill = L1)) +
    geom_histogram(color = 1) +
    labs(
      y = "Frequency\n",
      x = "\nPercentages",
      fill = element_blank()
    ) +
    scale_fill_manual(values=c("royalblue1", "tomato1"), labels = c("Negative", "Positive")) +
    theme_classic() +
    theme(text = element_text(size=18), panel.grid.major.y = element_line(color = 1, linetype = 'dashed')) +
    scale_y_continuous(expand = c(0, 0))
  
  navvac_hist # this line outputs the plot if using a console/RStudio
}
