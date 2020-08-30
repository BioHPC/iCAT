#' Main training function.
#' @param negatives Dataframe of sequence frequencies in negative samples.
#' @param positives Datafram of sequence frequencies in positive samples.
#' @param prelist Vector of negative training samples.
#' @param postlist Vector of positive training samples.
#' @param field String containing the column or columns (space-delimited) of interest.
#' @param count String containing the column name for colontype counts.
#' @param copyrange Integer Vector of the min and max copy of a sequence, within a sample, to be considered.
#' @param pcut P-value threshold for fisher-exact test.
#' @param minpublic Sequence frequency threshold to be considered.
#' @param updateProgress Function for updating a progress bar in a Shiny interface.
#' @return List containing both negtive (n) and positive (v) clonotype percentages.
#' @export
#' @examples
#' FIELD <- "vGeneName aminoAcid jGeneName"
#' COUNT <- "copy"
#' P_CUTOFF <- 0.1
#' MIN_PUBLIC <- 2
#' COPY_RANGE <- "1 99"
#' 
#' listPos <- tsvDir(system.file("extdata", "Post", package="iCAT"))
#' listNeg <- tsvDir(system.file("extdata", "Pre", package="iCAT"))
#' 
#' naive <- readTrn(listNeg, FIELD, COUNT, COPY_RANGE, "naive")
#' vaccs <- readTrn(listPos, FIELD, COUNT, COPY_RANGE, "vacc")  
#' 
#' mod <- train(naive, vaccs, listNeg, listPos, FIELD, COUNT, COPY_RANGE, P_CUTOFF, MIN_PUBLIC, NULL)
train <- function(negatives, positives, prelist, postlist, field, count, copyrange, pcut, minpublic, updateProgress=NULL) {

  fs <- strsplit(field, ' ')[[1]]
  copyrange <- as.integer(strsplit(copyrange, ' ')[[1]])

  numpre <- length(prelist)
  numpost <- length(postlist)

  colnames(negatives)[2] <- "naiveamounts"
  colnames(positives)[2] <- "vaccamounts"
  #merge vaccinated and unvaccinated samples
  all <- merge(positives, negatives, all.x = TRUE)
  all[is.na(all)] <- 0


  #calculate and eliminate samples that are upregulated in naive
  totals <- (all$vaccamounts * (numpre / numpost) - all$naiveamounts)

  all <- data.table::data.table(cbind(all, totals))
  all <- all[order(totals),]
  all <- subset(all, totals > 0)
  #calculate number of samples that each amino acid is absent from
  naiveabsent <- numpre - all$`naiveamounts`
  vaccabsent <- numpost - all$`vaccamounts`
  #bind those numbers to the rest of the table
  all <- cbind(all, naiveabsent, vaccabsent)

  #pull out fisher-test relevant values and run test
  fishervalues <- all[, c(2, 3, 6, 5)]
  if (!is.null(updateProgress))
    updateProgress(detail = "Calculating p-value")
  # slow step...
  # instead of repeating fisher tests only do the unique ones then merge back
  uniquefishervalues <- unique(fishervalues[, 1:4])
  uniquefishervalues$pvals <-
    apply(uniquefishervalues, 1, function(x)
      fisher.test(matrix(x, nrow = 2), alternative = "greater")$p.value)
  # AR: **kyle sumcv idea**
  pvalfishies <-
    merge(
      uniquefishervalues,
      fishervalues,
      by = c("vaccamounts", "naiveamounts", "vaccabsent",   "naiveabsent"),
      sort = FALSE
    )


  uniquefishervalues$sumcv <-
    lapply(uniquefishervalues$pvals, function (x) sum(pvalfishies$vaccamounts[pvalfishies$pvals <= x]))
  uniquefishervalues$covvac <-
    as.numeric(uniquefishervalues$sumcv) / numpost
  uniquefishervalues$sumcn <-
    lapply(uniquefishervalues$pvals, function (x) sum(pvalfishies$naiveamounts[pvalfishies$pvals <= x]))
  uniquefishervalues$covnav <- as.numeric(uniquefishervalues$sumcn) / numpre



  #bind p values to table
  all <- cbind(all, pvalfishies)

  #all <- all[order(pvalfishies)]
  #grep for only viable sequences
  all <- all[c(grep("^[A-Z*]", all$names)),]

  #all <- transform(all, Cv = ave(vaccamounts, FUN = cumsum))

  #all <- transform(all, Cn = ave(naiveamounts, FUN = cumsum))

  uniquefishervalues$ratio <- ifelse(uniquefishervalues$covnav < 1,
                                     uniquefishervalues$covvac,
                                     uniquefishervalues$covvac/uniquefishervalues$covnav)


  testvalue <- subset(uniquefishervalues, uniquefishervalues$pvals < as.double(pcut))


  pval <- testvalue$pvals[which.max(testvalue$ratio)]

  #get final list
  if (!is.null(updateProgress))
    updateProgress(detail = "Separating sequence library")
  finally <- all[all$pvals <= pval,]
  finally <- finally[finally$vaccamounts >= as.numeric(minpublic),]

  lib <- finally


  greplistpre <- lapply(prelist, function(x) {
    dat <- data.table::fread(x, select = c(fs, count))
    dat <- dat[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
    dat <- dat[, names := do.call(paste,.SD), .SDcols=!count]
    dat <- dat[, c(fs, count) := NULL]

    dat <- dat[c(grep("^[A-Z*]", dat$names)), ]
    freq <- NULL
    dat <- unique(dat[,freq := .N, by = names]) # similar to sort | uniq -c
    h <- hash::hash(dat$names, dat$freq)
    return(h)
  })


  greplistppost <- lapply(postlist, function(x) {
    dat <- data.table(fread(x, select = c(fs,count)))
    dat <- dat[get(count) >= copyrange[1] & get(count) <= copyrange[2]]
    dat <- dat[, names := do.call(paste,.SD), .SDcols=!count]
    dat <- dat[, c(fs, count) := NULL]
    
    dat <- dat[c(grep("^[A-Z*]", dat$names)), ]
    freq <- NULL
    dat <- unique(dat[,freq := .N, by = names]) # similar to sort | uniq -c
    h <- hash::hash(dat$names, dat$freq) # build hash of counts
    return(h)
  })



  #find percentage of significant clonotypes
  navcounts <- list()
  for (i in 1:numpre) {
    if (!is.null(updateProgress))
      updateProgress(detail = paste0("Finding % of significant clonotypes [Negative Sample #", i, "]"))
    navcounts[i] <-
    sum(unlist(
      lapply(finally$names, function(x)
        greplistpre[[i]][[x]])
    ))
  }

  navtotals <- sapply(greplistpre, function(x) sum(hash::values(x)))

  navpercs <- (as.numeric(navcounts) / navtotals) * 100

  vaccounts <- list()
  for (i in 1:numpost) {
    if (!is.null(updateProgress))
      updateProgress(detail = paste0("Finding % of significant clonotypes [Positive Sample #", i, "]"))
    vaccounts[i] <-
    sum(unlist(
      lapply(finally$names, function(x)
        greplistppost[[i]][[x]])
    ))
  }

  vactotals <- lapply(greplistppost, function(x) sum(hash::values(x)))

  vacpercs <- (as.numeric(vaccounts) / as.numeric(vactotals)) * 100

  return(list(n=navpercs, v=vacpercs, l=finally))
}
