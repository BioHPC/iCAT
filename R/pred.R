#' Predict the exposure status of an independent sample.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @param indpt Vector of independent samples file paths.
#' @param names Vector of labels for independent samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @param count String containing the column name for colontype counts.
#' @param copyrange Integer Vector of the min and max copy of a sequence, within a sample, to be considered.
#' @return Matrix with \% correct predictions from training data.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' COUNT <- "copy"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' COPY_RANGE <- "1 99"
#' 
#' 
#' listPos <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, COUNT, COPY_RANGE, "naive")
#' vaccs <- readTrn(listPos, FIELD, COUNT, COPY_RANGE, "vacc")  
#' 
#' mod <- train(naive, vaccs, listNeg, listPos, FIELD, COUNT, COPY_RANGE, P_CUTOFF, MIN_PUBLIC, NULL)
#' pred(mod, system.file("extdata", "Pre", "KJW100_HLA-A2_10_PRE.tsv", package="iCAT"), "unknown-sample-label", FIELD, COUNT, COPY_RANGE)
pred <- function(comb, indpt, names, field, count, copyrange) {
  fs <- strsplit(field, ' ')[[1]]
  copyrange <- as.integer(strsplit(copyrange, ' ')[[1]])
  
  nums <- length(indpt)

  greplistppost <- lapply(indpt, function(x) {
    dat <- data.table(fread(x, select = c(fs,count)))
    dat <- dat[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
    dat <- dat[, names := do.call(paste,.SD), .SDcols=!count]
    dat <- dat[, c(fs, count) := NULL]
    
    # x <- rep("", length(dat[[1]]))
    # for (i in 1:length(fs))
    #   x <- paste(x, dat[[i]])
    # x <- substring(x, 2)
    # dat <-
    #   data.table::data.table(x)
    dat <- dat[c(grep("^[A-Z*]", dat$names)), , drop=FALSE]
    freq <- NULL
    dat <- unique(dat[,freq := .N, by = names], drop=FALSE) # similar to sort | uniq -c
    h <- hash::hash(dat$names, dat$freq) # build hash of counts
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
  
  e <- exp(1)
  x <- ifelse(idnavdnorm > idvacdnorm, idnavdnorm - idvacdnorm, idvacdnorm - idnavdnorm)

  df <- cbind(basename(names), class, idpercs, 100*(1-(e^x/((e^x)+1))))
              
  colnames(df) <- c("Sample", "Prediction", "% TARS", "% Uncertainty")

  return(df)
}
