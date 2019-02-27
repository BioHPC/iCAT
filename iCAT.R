library(data.table)
library(dplyr)
library(plyr)
library(hash)
library(ggplot2)

start_time <- Sys.time()
r <- c()
lib <- c()
n <- c()
p <- c()

trnStats <- function(listPre, listPost, field) {
  fs <- strsplit(field, ' ')[[1]]
  
  pre <- rbindlist(lapply(listPre, function(x) {
    dat <- fread(x, select = c(fs, "count (templates)", "copy"))
  }))
  
  DTpre <- as.data.table(pre)
  pre <- DTpre[, lapply(.SD, sum), by = fs, .SDcols = "count (templates)"]
  pref <- c(length(listPre), length(pre$`count (templates)`), sum(pre$`count (templates)`))
  
  post <- rbindlist(lapply(listPost, function(x) {
    dat <- fread(x, select = c(fs, "count (templates)"))
  }))
  
  DTpost <- as.data.table(post)
  post <- DTpost[, lapply(.SD, sum), by = fs, .SDcols = "count (templates)"]
  postf <- c(length(listPost), length(post$`count (templates)`), sum(post$`count (templates)`))
  
  both <- rbind(pref, postf)
  colnames(both) <- c("# Samples", "# Clonotypes", "# Unique Sequences")
  rownames(both) <- c("Negative", "Positive")
  
  return(both)
}

#process files
processFiles <- function(x, field) {
  fs <- strsplit(field, ' ')[[1]]
  # get sequenceStatus column to use as filter, then remove it
  #dat <- readLines(x)
  #dat <- dat[grep('^[A-Za-z]', dat)]
  #dat <- read.table(textConnection(dat), sep="\t", header=T)
  #dat <- select(dat, fs, "sequenceStatus")
  dat <- fread(x, select = c(fs, "sequenceStatus"))
  dat <- dat[dat$sequenceStatus == "In",]
  dat <- dat[, 1:length(fs)]
  dat <- dat[!duplicated(dat),]
  dat <- dat[order(dat[,1:length(fs)]),]
  #print(dat)
}

listDir <- function (dir) {
  list <-
    list.files(
      path = dir,
      recursive = TRUE,
      full.names = TRUE,
      pattern = "*.tsv"
    )
}

readPre <- function(list, field) {
  if (length(list) == 0) {
    return(NULL)
  }
  
  fs <- strsplit(field, ' ')[[1]]
  #combine files and count repeats
  
  final <- lapply(list, processFiles, field=field)
  final <- rbindlist(final)
  x <- rep("", length(final[[1]]))
  for (i in 1:length(fs))
    x <- paste(x, final[[i]])
  x <- substring(x, 2)
  final <-
    data.table(x)
  colnames(final) <- "names"
  final <- final[, .N, by = names(final)]
  colnames(final) <- c("names", "naiveamounts")
  return(final)
}

readPost <- function(list, field) {
  if (length(list) == 0) {
    return(NULL)
  }
  
  #combine files and count repeats
  fs <- strsplit(field, ' ')[[1]]
  final <- lapply(list, processFiles, field=field)
  final <- rbindlist(final)
  x <- rep("", length(final[[1]]))
  for (i in 1:length(fs))
    x <- paste(x, final[[i]])
  x <- substring(x, 2)
  final <-
    data.table(x)
  colnames(final) <- "names"
  final <- final[, .N, by = names(final)]
  colnames(final) <- c("names", "vaccamounts")
  return(final)
}


analyse <- function(naive, vaccs, prelist, postlist, field, pcut, minpublic, updateProgress) {
  
  fs <- strsplit(field, ' ')[[1]]
  
  numpre <- length(prelist)
  numpost <- length(postlist)
  
  #merge vaccinated and unvaccinated samples
  all <- merge(vaccs, naive, all.x = T)
  all[is.na(all)] <- 0
  
  
  
  #calculate and eliminate samples that are upregulated in naive
  totals <-
    (all$vaccamounts * (numpre / numpost) - all$naiveamounts)
  all <- data.table(cbind(all, totals))
  all <- all[order(totals),]
  all <- subset(all, totals > 0)
  #calculate number of samples that each amino acid is absent from
  naiveabsent <- numpre - all$`naiveamounts`
  vaccabsent <- numpost - all$`vaccamounts`
  #bind those numbers to the rest of the table
  all <- cbind(all, naiveabsent, vaccabsent)

  
  #pull out fisher-test relevant values and run test
  fishervalues <- all[, c(2, 3, 6, 5)]
  #updateProgress(detail = "Calculating p-value")
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
  #updateProgress(detail = "Separating sequence library")
  finally <- all[all$pvals <= pval,]
  finally <- finally[finally$vaccamounts >= as.numeric(minpublic),]

  lib <<- finally
  

  greplistpre <- lapply(prelist, function(x) {
    dat <- fread(x, select = c(fs))
    x <- rep("", length(dat[[1]]))
    for (i in 1:length(fs))
      x <- paste(x, dat[[i]])
    x <- substring(x, 2)
    dat <-
      data.table(x)
    dat <- dat[c(grep("^[A-Z*]", dat$x)), ]
    dat <- unique(dat[,freq := .N, by = x]) # similar to sort | uniq -c
    h <- hash(dat$x, dat$freq)
    return(h)
  })
  
  
  greplistppost <- lapply(postlist, function(x) {
    dat <- fread(x, select = fs)
    x <- rep("", length(dat[[1]]))
    for (i in 1:length(fs))
      x <-paste(x, dat[[i]])
    x <- substring(x, 2)
    dat <-
      data.table(x)
    dat <- dat[c(grep("^[A-Z*]", dat$x)), ]
    dat <- unique(dat[,freq := .N, by = x]) # similar to sort | uniq -c
    h <- hash(dat$x, dat$freq) # build hash of counts
    return(h)
  })
  
  
  
  #find percentage of significant clonotypes
  navcounts <- list()
  for (i in 1:numpre) {
    #updateProgress(detail = paste0("Finding % of significant clonotypes [Negative Sample #", i, "]"))
    navcounts[i] <-
    print(sum(unlist(
      lapply(finally$names, function(x)
        greplistpre[[i]][[x]])
    )))
  }

  navtotals <- sapply(greplistpre, function(x) sum(values(x)))
  
  navpercs <- (as.numeric(navcounts) / navtotals) * 100
  
  vaccounts <- list()
  for (i in 1:numpost) {
    #updateProgress(detail = paste0("Finding % of significant clonotypes [Positive Sample #", i, "]"))
    vaccounts[i] <-
    print(sum(unlist(
      lapply(finally$names, function(x)
        greplistppost[[i]][[x]])
    )))
  }
  
  vactotals <- lapply(greplistppost, function(x) sum(values(x)))
  
  vacpercs <- (as.numeric(vaccounts) / as.numeric(vactotals)) * 100
  r <<- list(n=navpercs, v=vacpercs )
  return(list(n=navpercs, v=vacpercs))
}

# swithed to ggplot for flexibility
plotHist <- function(comb) {
  navvac_hist <-
    ggplot(melt(comb), aes(x = value, fill = L1)) +
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

savePlot <- function(file, plotIn) {
  ggsave(
    device = "pdf",
    plot = plotIn,
    width = 5,
    height = 5,
    units = "in",
    dpi = 200
  )
}

# Kyle's Method (modified to compile)
classMat <- function(comb) {
  navpercs <- comb$n
  vacpercs <- comb$v
  

  navmean <- mean(navpercs)
  navmean <- navmean * 3 # accounting for overfitting
  
  navsd <- sd(navpercs)
  navsd <- navsd * 3 # accounting for overfitting
  
  vacmean <- mean(vacpercs)
  
  vacsd <- sd(vacpercs)
  
  navnavdnorm <- dnorm(navpercs,
                       mean = navmean,
                       sd = navsd,
                       log = FALSE)
  
  navvacdnorm <- dnorm(navpercs,
                       mean = vacmean,
                       sd = vacsd,
                       log = FALSE)
  
  vacvacdnorm <- dnorm(vacpercs,
                       mean = vacmean,
                       sd = vacsd,
                       log = FALSE)
  
  vacnavdnorm <- dnorm(vacpercs,
                       mean = navmean,
                       sd = navsd,
                       log = FALSE)
  
  navclass <- ifelse(navnavdnorm > navvacdnorm, "NAIVE", "INFECTED")
  
  vacclass <- ifelse(vacvacdnorm > vacnavdnorm, "INFECTED", "NAIVE")
  
  navnavdnorm <- round(navnavdnorm, digits = 4)
  
  navvacdnorm <- round(navvacdnorm, digits = 4)
  
  vacvacdnorm <- round(vacvacdnorm, digits = 4)
  
  vacnavdnorm <- round(vacnavdnorm, digits = 4)
  
  navdata <- cbind(navpercs, navnavdnorm, navvacdnorm, navclass)
  
  vacdata <- cbind(vacpercs, vacvacdnorm, vacnavdnorm, vacclass)
  
  classmatrix <-
    matrix(
      c(
        sum(navdata[, 4] == "NAIVE"),
        sum(navdata[, 4] == "INFECTED"),
        sum(vacdata[, 4] == "NAIVE"),
        sum(vacdata[, 4] == "INFECTED")
      ),
      nrow = 2,
      ncol = 2,
      byrow = TRUE
    )
  
  navclasscorrect <- (classmatrix[1, 1] / sum(classmatrix[1, ])) * 100
  
  vacclasscorrect <- (classmatrix[2, 2] / sum(classmatrix[2, ])) * 100
  
  classmatrix <-
    matrix(
      c(classmatrix, navclasscorrect, vacclasscorrect),
      nrow = 2,
      ncol = 3
    )
  
  rownames(classmatrix) <- c("Negative", "Positive")
  
  colnames(classmatrix) <- c("Unexposed", "Exposed", "% Correct")
  
  classmatrix
}

getLib <- function() {
  l <- lib[, unique(colnames(lib)), with=FALSE]
  colnames(l) <- c("Sequence",	"Positive Present",	"Negative Present",	"Total",	"Negative Absent",	"Positive Absent",	"P-Value")
  l <- l[, c(1, 2, 6, 3, 5, 7), with=FALSE]
  return(l)
}

pred <- function(comb, iccs, indpt, names, field) {
  fs <- strsplit(field, ' ')[[1]]
  nums <- length(indpt)
  
  greplistppost <- lapply(indpt, function(x) {
    dat <- fread(x, select = fs)
    x <- rep("", length(dat[[1]]))
    for (i in 1:length(fs))
      x <- paste(x, dat[[i]])
    x <- substring(x, 2)
    dat <-
      data.table(x)
    dat <- dat[c(grep("^[A-Z*]", dat$x)), ]
    dat <- unique(dat[,freq := .N, by = x]) # similar to sort | uniq -c
    h <- hash(dat$x, dat$freq) # build hash of counts
    return(h)
  })
  
print(greplistppost)
  idcounts <- list()
  for (i in 1:nums)
    idcounts[i] <-
    print(sum(unlist(
      lapply(lib$names, function(x)
        greplistppost[[i]][[x]])
    )))
  
  idtotals <- lapply(greplistppost, function(x) sum(values(x)))
  
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
