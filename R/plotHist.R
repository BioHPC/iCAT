#' Plot histogram of clonotype frequencies in negative and positive samples.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @return ggplot histogram of clonotype frequencies in negative and positive samples
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' 
#' listPos <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, "naive")
#' vaccs <- readTrn(listPos, FIELD, "vacc")  
#' 
#' mod <- train(naive, vaccs, listNeg, listPos, FIELD, P_CUTOFF, MIN_PUBLIC, NULL)
#' plotHist(mod)
utils::globalVariables(c("Pvalue","Pid"))
plotHist <- function(comb) {
  comb_ <- list(n=comb$n, v=comb$v)
  m <- melt(comb_)
  colnames(m) <- c('Pvalue', 'Pid')
  navvac_hist <-
    ggplot2::ggplot(m, ggplot2::aes(x = Pvalue, fill = Pid)) +
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
