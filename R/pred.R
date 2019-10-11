#' Predict the exposure status of an independent sample.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @param indpt Vector of independent samples file paths.
#' @param names Vector of labels for independent samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @return Matrix with \% correct predictions from training data.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' 
#' listPos <- tsvDir("iCAT/extdata/Post/")
#' listNeg <- tsvDir("iCAT/extdata/Pre/")
#' 
#' naive <- readTrn(listNeg, FIELD, "naive")
#' vaccs <- readTrn(listPos, FIELD, "vacc")  
#' 
#' mod <- train(naive, vaccs, listNeg, listPos, FIELD, P_CUTOFF, MIN_PUBLIC, NULL)
#' pred(mod, "iCAT/extdata/Post/KJW122_M83_D08.tsv", "unknown-sample-label", FIELD)
pred <- function(comb, indpt, names, field) {
  fs <- strsplit(field, ' ')[[1]]
  nums <- length(indpt)

  greplistppost <- lapply(indpt, function(x) {
    dat <- fread(x, select = fs)
    x <- rep("", length(dat[[1]]))
    for (i in 1:length(fs))
      x <- paste(x, dat[[i]])
    x <- substring(x, 2)
    dat <-
      data.table::data.table(x)
    dat <- dat[c(grep("^[A-Z*]", dat$x)), , drop=FALSE]
    dat <- unique(dat[,freq := .N, by = x], drop=FALSE) # similar to sort | uniq -c
    h <- hash::hash(dat$x, dat$freq) # build hash of counts
    return(h)
  })
  
  lib <- comb$l
  idcounts <- list()
  for (i in 1:nums)
    idcounts[i] <-
    sum(unlist(
      lapply(lib$names, function(x)
        greplistppost[[i]][[x]])
    ))

  idtotals <- lapply(greplistppost, function(x) sum(hash::values(x)))

  idpercs <- (as.numeric(idcounts) / as.numeric(idtotals)) * 100

  navpercs <- comb$n
  vacpercs <- comb$v

  navmean <- mean(navpercs)
  navmean <- navmean * 3 # accounting for overfitting

  navsd <- sd(navpercs)
  navsd <- navsd * 3 # accounting for overfitting

  vacmean <- mean(vacpercs)

  vacsd <- sd(vacpercs)

  idnavdnorm <- dnorm(idpercs,
                       mean = navmean,
                       sd = navsd,
                       log = FALSE)

  idvacdnorm <- dnorm(idpercs,
                       mean = vacmean,
                       sd = vacsd,
                       log = FALSE)

  class <- ifelse(idnavdnorm > idvacdnorm, "Negative", "Positive")


  df <- cbind(basename(names), class, idpercs)
  colnames(df) <- c("Sample", "Prediction", "% TARS")

  return(df)
}
