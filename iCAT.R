library(data.table)
library(dplyr)
library(plyr)
library(hash)
library(ggplot2)

NaiveFiles <- "Pre"
VaccFiles <- "Post"

# iCAT <- function(NaiveFiles, VaccFiles) {

start_time <- Sys.time()

#process files
processFiles <- function(x) {
  # get sequenceStatus column to use as filter, then remove it
  dat <-
    fread(x,
          select = c("aminoAcid", "vGeneName", "jGeneName", "sequenceStatus"))
  dat <- dat[dat$sequenceStatus == "In",]
  dat <- dat[, 1:3]
  dat <- dat[!duplicated(dat),]
  dat <- dat[order(dat$aminoAcid),]
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

readPre <- function(list) {
  if (length(list) == 0) {
    message("Error in locating Naive Files")
    return()
    
  }
  
  #combine files and count repeats
  
  final <- lapply(list, processFiles)
  final <- rbindlist(final)
  final <-
    data.table(paste(final$aminoAcid, final$vGeneName, final$jGeneName, sep = " "))
  colnames(final) <- "names"
  final <- final[, .N, by = names(final)]
  colnames(final) <- c("names", "naiveamounts")
  return(final)
}

readPost <- function(list) {
  if (length(list) == 0) {
    message("Error in locating Naive Files")
    return()
    
  }
  
  #combine files and count repeats
  
  final <- lapply(list, processFiles)
  final <- rbindlist(final)
  final <-
    data.table(paste(final$aminoAcid, final$vGeneName, final$jGeneName, sep = " "))
  colnames(final) <- "names"
  final <- final[, .N, by = names(final)]
  colnames(final) <- c("names", "vaccamounts")
  return(final)
}


# processtime <- Sys.time()
#
# elapsed_time = round(difftime(processtime, start_time, units = "secs"), digits = 2)
# message("Process time: ", elapsed_time, " seconds")


analyse <- function(naive, vaccs, prelist, postlist) {
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
  
  # slow step...
  # instead of repeating fisher tests only do the unique ones then merge back
  uniquefishervalues <- unique(fishervalues[, 1:4])
  uniquefishervalues$pvals <-
    apply(uniquefishervalues, 1, function(x)
      fisher.test(matrix(x, nrow = 2), alternative = "greater")$p.value)
  
  pvalfishies <-
    merge(
      uniquefishervalues,
      fishervalues,
      by = c("vaccamounts", "naiveamounts", "vaccabsent",   "naiveabsent"),
      sort = FALSE
    )
  
  #bind p values to table
  all <- cbind(all, pvalfishies)
  #all <- all[order(pvalfishies)]
  #grep for only viable sequences
  all <- all[c(grep("^[A-Z*]", all$names)),]
  all <- transform(all, Cv = ave(vaccamounts, FUN = cumsum))
  
  all <- transform(all, Cn = ave(naiveamounts, FUN = cumsum))
  
  covvac <- all$Cv / numpost
  covnav <- all$Cn / numpre
  all <- cbind(all, covvac, covnav)
  all$covnav[all$covnav < 1] <- 1
  ratio <- all$covvac / all$covnav
  all <- cbind(all, ratio)
  all <- all[all$pvals < .1,]
  #all <- subset(all, pvalfishies < .1)
  maxrat <- max(all$ratio)
  pval <- all$pvalfishies[ratio == maxrat]
  #get final list
  finally <- all[all$pvals <= .1,]
  #  finally <<- subset(all, pvalfishies <= .1)
  
  # greplistpre <- lapply(prelist, function(x) {
  #   f  <- file(x, open = "r")
  #   h <- hash()
  #   while (length(oneLine <- readLines(f, n = 1, warn = FALSE)) > 0) {
  #     # aminoAcid, vGeneName, jGeneName
  #     id <- paste((strsplit(oneLine, "\t", fixed=TRUE))[[1]][c(2,8,22)], collapse=" ") 
  #     if (!startsWith(id, " ")) {
  #       
  #       if (has.key(id, h)) {
  #         h[[id]] <- h[[id]] + 1
  #       } else {
  #         h[[id]] <- 1
  #       }
  #     }
  #   }
  #   close(f)
  # 
  #   return(h)
  # })
  greplistpre <- lapply(prelist, function(x) {
    dat <- fread(x, select = c("aminoAcid", "vGeneName", "jGeneName"))
    dat <-
      data.table(paste(dat$aminoAcid, dat$vGeneName, dat$jGeneName, sep = " "))
    dat <- dat[c(grep("^[A-Z*]", dat$V1)), ]
    dat <- unique(dat[,freq := .N, by = V1]) # similar to sort | uniq -c
    h <- hash(dat$V1, dat$freq)
    return(h)
  })
  
  
  greplistppost <- lapply(postlist, function(x) {
    dat <- fread(x, select = c("aminoAcid", "vGeneName", "jGeneName"))
    dat <-
      data.table(paste(dat$aminoAcid, dat$vGeneName, dat$jGeneName, sep = " "))
    dat <- dat[c(grep("^[A-Z*]", dat$V1)), ]
    dat <- unique(dat[,freq := .N, by = V1]) # similar to sort | uniq -c
    h <- hash(dat$V1, dat$freq) # build hash of counts
    return(h)
  })
  
  
  
  #find percentage of significant clonotypes
  navcounts <- list()
  for (i in 1:numpre)
    navcounts[i] <-
    print(sum(unlist(
      lapply(finally$names, function(x)
        greplistpre[[i]][[x]])
    )))

  navtotals <- sapply(greplistpre, function(x) sum(values(x)))
  
  navpercs <- (as.numeric(navcounts) / navtotals) * 100
  
  vaccounts <- list()
  for (i in 1:numpost)
    vaccounts[i] <-
    print(sum(unlist(
      lapply(finally$names, function(x)
        greplistppost[[i]][[x]])
    )))
  
  vactotals <- lapply(greplistppost, function(x) sum(values(x)))
  
  vacpercs <- (as.numeric(vaccounts) / as.numeric(vactotals)) * 100
  return(cbind(navpercs, vacpercs))
}

# swithed to ggplot for flexibility
plotHist <- function(comb) {
  navvac_hist <-
    ggplot(melt(comb), aes(x = value, fill = Var2)) +
    geom_histogram(binwidth = 0.01, color = 1) +
    labs(
      title = "Percentage of Clonotypes",
      y = "Frequency",
      x = "Percentages",
      fill = element_blank()
    ) +
    scale_fill_discrete(labels = c("Naive", "Vaccinated")) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = 1, linetype = 'dashed')) +
    scale_y_continuous(expand = c(0, 0))
  
  navvac_hist # this line outputs the plot if using a console/RStudio
  
}

savePlot <- function(plotIn) {
  ggsave(
    filename = "HistClonotypes.pdf",
    plot = plotIn,
    width = 5,
    height = 5,
    units = "in",
    dpi = 200
  )
}

# ** start of Kyle's Method (modified to compile) **
# navmean <- mean(navpercs)
#
# navsd <- sd(navpercs)
#
# vacmean <- mean(vacpercs)
#
# vacsd <- sd(vacpercs)
#
# navnavdnorm <- dnorm(navpercs, mean=navmean, sd=navsd, log=FALSE)
#
# navvacdnorm <- dnorm(navpercs, mean=vacmean, sd=vacsd, log=FALSE)
#
# vacvacdnorm <- dnorm(vacpercs, mean=vacmean, sd=vacsd, log=FALSE)
#
# vacnavdnorm <- dnorm(vacpercs, mean=navmean, sd=navsd, log=FALSE)
#
# navclass <- ifelse(navnavdnorm>navvacdnorm, "NAÏVE", "INFECTED")
#
# vacclass <- ifelse(vacvacdnorm>vacnavdnorm, "INFECTED", "NAÏVE")
#
# navnavdnorm <- round(navnavdnorm, digits = 4)
#
# navvacdnorm <- round(navvacdnorm, digits = 4)
#
# vacvacdnorm <- round(vacvacdnorm, digits = 4)
#
# vacnavdnorm <- round(vacnavdnorm, digits = 4)
#
# navdata <- cbind(navpercs, navnavdnorm, navvacdnorm, navclass)
#
# vacdata <- cbind(vacpercs, vacvacdnorm, vacnavdnorm, vacclass)
#
# classmatrix <-
#   matrix(
#     c(
#       sum(navdata[,4] == "NAÏVE"), # assuming this is what intended
#       sum(navdata[,4] == "INFECTED"),
#       sum(vacdata[,4] == "NAÏVE"),
#       sum(vacdata[,4] == "INFECTED")
#     ),
#     nrow = 2,
#     ncol = 2,
#     byrow = TRUE
#   )
#
# navclasscorrect <- (classmatrix[1,1]/sum(classmatrix[1,]))*100
#
# vacclasscorrect <- (classmatrix[2,2]/sum(classmatrix[2,]))*100
#
# classmatrix <- matrix(c(classmatrix, navclasscorrect, vacclasscorrect), nrow = 2, ncol = 3)
#
# rownames(classmatrix) <- c("NAÏVE", "INFECTED")
#
# colnames(classmatrix) <- c("UNEXPOSED", "EXPOSED", "% CORRECT")

# ** end of Kyle's method **

# end_time <- Sys.time()
# elpased_time = end_time - start_time
# cat("Total time: ", elpased_time)
