library(data.table)
library(dplyr)
library(plyr)

NaiveFiles <- "Pre"
VaccFiles <- "Post"

# iCAT <- function(NaiveFiles, VaccFiles) {

  start_time <- Sys.time()
  prelist <- list.files(path = NaiveFiles, recursive = TRUE, full.names = TRUE, pattern = "*.tsv")
  if ( length(prelist) == 0 ) {
    message("Error in locating Naive Files")
    return();
  }

  #process files
  processFiles <- function(x){
    dat <- fread(x, select = c("aminoAcid", "vGeneName", "jGeneName"))
    dat <- dat[!duplicated(dat), ]
    dat <- dat[order(dat$aminoAcid), ]
  }

  #combine files and count repeats
  naive <- lapply(prelist, processFiles)
  naive <- rbindlist(naive)
  naive <- data.table(paste(naive$aminoAcid, naive$vGeneName, naive$jGeneName, sep = " "))
  colnames(naive) <- "names"
  naive <- naive[,.N, by = names(naive)]
  colnames(naive) <- c("names", "naiveamounts")

  postlist = list.files(path = VaccFiles, recursive = TRUE, full.names = TRUE,pattern = '*.tsv')
  if ( length(postlist) == 0 ) {
    message("Error in locating Vacc Files")
    return();
  }

  #combine files and condense/count repeats
  vaccs <- lapply(postlist, processFiles)
  vaccs <- rbindlist(vaccs)
  vaccs <- data.table(paste(vaccs$aminoAcid, vaccs$vGeneName, vaccs$jGeneName, sep = " "))
  colnames(vaccs) <- "names"
  vaccs <- vaccs[,.N, by = names(vaccs)]
  colnames(vaccs) <- c("names", "vaccamounts")

  numpre <- length(prelist)
  numpost <- length(postlist)


  processtime <- Sys.time()
  elapsed_time = round(difftime(processtime, start_time, units = "secs"), digits = 2)
  message("Process time: ", elapsed_time, " seconds")

  #merge vaccinated and unvaccinated samples
  all <- merge(vaccs,naive, all.x = T)
  all[is.na(all)] <- 0



  #calculate and eliminate samples that are upregulated in naive
  totals <- (all$`vaccamounts`*(numpre/numpost) - all$`naiveamounts`)
  all <- data.table(cbind(all, totals))
  all <- all[order(totals),]
  all <- subset(all, totals > 0)
  #calculate number of samples that each amino acid is absent from
  naiveabsent <- numpre - all$`naiveamounts`
  vaccabsent <- numpost - all$`vaccamounts`
  #bind those numbers to the rest of the table
  all <- cbind(all, naiveabsent, vaccabsent)



  #pull out fisher-test relevant values and run test
  fishervalues <- all[,c(2,3,6,5)]

  # slow step...
  # instead of repeating fisher tests only do the unique ones then merge back
  uniquefishervalues <- unique(fishervalues[,1:4])
  uniquefishervalues$pvals <- apply(uniquefishervalues,1, function(x) fisher.test(matrix(x,nrow = 2), alternative = "greater")$p.value)

  pvalfishies <- merge(uniquefishervalues, fishervalues, by = c("vaccamounts", "naiveamounts", "vaccabsent",   "naiveabsent"), sort = FALSE)

  #bind p values to table
  all <- cbind(all, pvalfishies)
  #all <- all[order(pvalfishies)]
  #grep for only viable sequences
  all <- all[c(grep("^[A-Z*]", all$names)),]
  all <- transform(all, Cv = ave(vaccamounts, FUN = cumsum))

  all <- transform(all, Cn = ave(naiveamounts, FUN = cumsum))

  covvac <- all$Cv/numpost
  covnav <- all$Cn/numpre
  all <- cbind(all, covvac, covnav)
  all$covnav[all$covnav < 1] <- 1
  ratio <- all$covvac/all$covnav
  all <- cbind(all,ratio)
  all <- all[all$pvals < .1,]
  #all <- subset(all, pvalfishies < .1)
  maxrat <- max(all$ratio)
  pval <- all$pvalfishies[ratio == maxrat]
  #get final list
  finally <- all[all$pvals <= .1,]
#  finally <<- subset(all, pvalfishies <= .1)

  VATStime <- Sys.time()
  elapsed_time = round(difftime(VATStime, start_time, units = "secs"), digits = 2)
  message("VATS time: ", elapsed_time, " seconds")

  greplistpre <- lapply(prelist, function(x){
    dat <- fread(x, select = c("aminoAcid", "vGeneName", "jGeneName"))
    dat <- data.table(paste(dat$aminoAcid, dat$vGeneName, dat$jGeneName, sep = " "))
    dat <- dat[c(grep("^[A-Z*]", dat$V1)),]
  })  
  
  greplistppost <- lapply(postlist, function(x){
    dat <- fread(x, select = c("aminoAcid", "vGeneName", "jGeneName"))
    dat <- data.table(paste(dat$aminoAcid, dat$vGeneName, dat$jGeneName, sep = " "))
    dat <- dat[c(grep("^[A-Z*]", dat$V1)),]
  })
  
  #find percentage of significant clonotypes
  navcounts<- list()
  for (i in 1:numpre) navcounts[i] <- print(sum(unlist(lapply(finally$names, function(x) length(grep(x, greplistpre[i]))))))
  
  navtotals <- list()
  navtotals <- sapply(greplistpre,nrow)
  
  navpercs <- (as.numeric(navcounts)/navtotals)*100
  
  vaccounts <- list()
  for (i in 1:numpost) vaccounts[i] <- print(sum(unlist(lapply(finally$names, function(x) length(grep(x, greplistppost[i]))))))
  
  vactotals <- lapply(greplistppost,nrow)
  
  vacpercs <- (as.numeric(vaccounts)/as.numeric(vactotals))*100
  
  hist(navpercs, xlim = c(0,.2),ylim = c(0,25), breaks = 2, main = "Percentages of clonotypes", xlab = "Percentages", col = 'pink')
  par(new=TRUE)
  hist(vacpercs, xlim = c(0,.2),ylim = c(0,25), breaks = 10,main = "Percentages of clonotypes", xlab = "Percentages", col = 'blue')
  
  end_time <- Sys.time()
  elpased_time = end_time - start_time
  cat("Total time: ", elpased_time)
  
