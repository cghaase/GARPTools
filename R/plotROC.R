#' Plots the Receiver Operating Charactersitic (ROC) Curve
#'
#' \code{plotROC} Plots the ROC curve of the best model subset output from GARP, relating 1-model specificity to model sensitivity
#'
#' @param n a numeric value specifying the number of models in the best subset output from GARP
#' @param x a raster object of the summated raster of the best models output from GARP
#' @param points a spatial object of presence data to use for testing locations
#'
#' @return A plot of the ROC curve.
#'
#' @details For discrete cutpoints (\code{n}), represents the number of models that agree on a predicted presence location.
#'
#' The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'
#' The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}.
#'
#' @seealso \code{\link{aucGARP}}
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   plotROC(n = 10, x = r, points = SpatialPoints(hs[,1:2]))
#'
#' @import raster
#' @import sp
#'
#' @export

plotROC <- function(n, x, points){

  #Extract values from summated grid at test point locations
  grid = x
  taxa.models <- extract(grid, points)
  if(any(is.na(taxa.models))){stop("One or more testing points do not overlap raster.")}

  #Create cutpoints dataframe
  cutpoints <- data.frame(seq(0,n,1))
  names(cutpoints) <- "cutpoint"

  #Summarize each cutpoint
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$taxa.present[i]   <- sum(taxa.models == i-1)
    cutpoints$cutpoint.area[i]  <- freq(grid, value = i-1)
    cutpoints$no.taxa.pixels[i] <- freq(grid, value = i-1) - sum(taxa.models == i-1)
    cutpoints$cum.area[i] <- ifelse(cutpoints$cutpoint[i] == 0, 0,
                                    ifelse(cutpoints$cutpoint[i] == 1,
                                           (cutpoints$no.taxa.pixels[1] + cutpoints$no.taxa.pixels[i]),
                                           (cutpoints$no.taxa.pixels[i] + cutpoints$cum.area[i-1])))
  }

  #Calculate confusion matrix
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$a[i] <- ifelse(cutpoints[i,1] == 0, (length(taxa.models) - cutpoints$taxa.present[i]), (cutpoints$a[i-1] - cutpoints$taxa.present[i]))
    cutpoints$b[i] <- ifelse(cutpoints[i,1] == 0, cutpoints$no.taxa.pixels[i], (cutpoints$b[i-1] + cutpoints$no.taxa.pixels[i]))
    cutpoints$c[i] <- ifelse(cutpoints[i,1] == 0, cutpoints$taxa.present[i], (cutpoints$c[i-1] + cutpoints$taxa.present[i]))
    cutpoints$d[i] <- ifelse(cutpoints[i,1] == 0, (sum(cutpoints$no.taxa.pixels) - cutpoints$no.taxa.pixels[i]), (cutpoints$d[i-1] - cutpoints$no.taxa.pixels[i]))
  }

  #Calculate sensitivity (true positives) and 1-specificity (false positives)
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$sensitivity[i]     <- cutpoints$a[i]/(cutpoints$a[i] + cutpoints$c[i])
    cutpoints$one.specificity[i] <- 1-(cutpoints$b[i]/(cutpoints$b[i] + cutpoints$d[i]))
  }

  #Calculate AUC value for each cutpoint
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$AUC[i] <- ifelse(cutpoints[i,1] == 0, (((1+cutpoints$sensitivity[i])/2) * (1-cutpoints$one.specificity[i])),
                               (((cutpoints$sensitivity[i] + cutpoints$sensitivity[i-1])/2) * (cutpoints$one.specificity[i-1] - cutpoints$one.specificity[i])))
  }

  #Plot ROC
  plot(cutpoints$one.specificity,cutpoints$sensitivity,
       main = paste("ROC Curve of ", n, " Best Models", sep = ""),
       type = "l", col = "blue", lwd = 2, tck = -0.02,
       xlab = "1-Specificity",
       ylab = "Sensitivity",
       ylim = c(0,1.0),
       xlim = c(0,1.0))
  lines(seq(max(cutpoints$one.specificity), 1, length.out = 2),
        c(max(cutpoints$sensitivity),max(cutpoints$sensitivity)), col = "blue", lwd = 2)
  lines(seq(0, 1, 0.1), seq(0, 1, 0.1),lty = 2, lwd = 1)
  points(cutpoints$one.specificity, cutpoints$sensitivity, pch = 16)
  legend("bottomright", legend = (c("Models", "Reference")), cex = 0.85,
         bg = "", bty = "n", col = c("blue", "black"), lty = c(1,2), lwd = c(2,2))
}
