#' Classification matrix estimating accuracy of model.
#'
#' @param comb List containing both negtive (n) and positive (v) clonotype percentages.
#' @return Matrix with % correct predictions from training data.
#' @export
#' @examples
#' comb <- train(...)
#' classMat(comb)
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
