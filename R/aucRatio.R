#' Ratio of Area Under the Receiver Operating Characteristic (ROC) curve (AUC)
#'
#' \code{aucRatio} Calculates the ratio of the area under the ROC curve (AUC) as suggested by Peterson et al. (2008). AUC represents the accuracy of the subset of the best models produced by GARP. AUC evaluates the predictive performance of the
#' best models by relating the model sensitivity (true positive rate) to 1-specificity (true negative rate), and can be described as the probability that any given cell is correctly predicted
#' as present or absent. An AUC of 0.5 is that predicted at random, while an AUC of 1 represents a perfect prediction
#'
#' @param n a numeric value specifying the number of models in the best subset outputted by GARP
#' @param x a raster object of the summated raster of the best models output by GARP
#' @param points a spatial object of presence data to use for testing locations
#' @param E the amount of error admissible along the true-positive axis (less than 1.0)
#'
#' @return Plots the modified Receiver Operating Characteristic (ROC) curve and returns a list containing:
#'  \itemize{
#'   \item \code{Modified.AUC} the total area under the modified ROC curve
#'   \item \code{Modified.Reference} the area under the modified reference curve
#'   \item \code{AUC.ratio} the ratio of the modified AUC to the modified reference value
#'  }
#'
#' @details The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'   The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}.
#'
#' @seealso \code{\link{aucGARP}}, \code{\link{zAUC}}, \code{\link{seAUC}}, \code{\link{wGARP}}, \code{\link{plotROC}}
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   aucRatio(n = 10, x = r, points = SpatialPoints(hs[,1:2]), E = 0.20)
#'
#' @references Peterson, A.T., M. Papes, J. Soberon (2007) Rethinking receiver operating characteristic analysis applications in ecological niche modeling. \emph{Ecological Modelling}. \bold{213(1):63-72}.
#'
#' @import raster
#' @import sp
#' @export

aucRatio <- function(n, x, points, E){

  #Extract values from summated grid at test point locations
  grid = x
  taxa.models <- raster::extract(grid, points)
  if(any(is.na(taxa.models))){stop("One or more testing points do not overlap raster.")}

  #Create cutpoints dataframe
  cutpoints <- data.frame(seq(0, n, 1))
  names(cutpoints) <- "cutpoint"

  #Summarize each cutpoint
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$taxa.present[i]   <- sum(taxa.models == i-1)
    cutpoints$cutpoint.area[i]  <- freq(grid, value = i-1)
    cutpoints$no.taxa.pixels[i] <- freq(grid, value = i-1) - sum(taxa.models == i-1)
    cutpoints$cum.area[i] <- ifelse(cutpoints$cutpoint[i] == 0,0,
                                    ifelse(cutpoints$cutpoint[i] == 1,
                                           (cutpoints$no.taxa.pixels[1] + cutpoints$no.taxa.pixels[i]),
                                           (cutpoints$no.taxa.pixels[i] + cutpoints$cum.area[i-1])))
  }

  #Calculate confusion matrix
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$a[i] <- ifelse(cutpoints[i,1] == 0, (length(taxa.models) - cutpoints$taxa.present[i]),(cutpoints$a[i-1] - cutpoints$taxa.present[i]))
    cutpoints$b[i] <- ifelse(cutpoints[i,1] == 0, cutpoints$no.taxa.pixels[i], (cutpoints$b[i-1] + cutpoints$no.taxa.pixels[i]))
    cutpoints$c[i] <- ifelse(cutpoints[i,1] == 0, cutpoints$taxa.present[i], (cutpoints$c[i-1] + cutpoints$taxa.present[i]))
    cutpoints$d[i] <- ifelse(cutpoints[i,1] == 0, (sum(cutpoints$no.taxa.pixels) - cutpoints$no.taxa.pixels[i]), (cutpoints$d[i-1] - cutpoints$no.taxa.pixels[i]))
  }

  #Calculate sensitivity and 1-specificity
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$sensitivity[i]     <- cutpoints$a[i]/(cutpoints$a[i] + cutpoints$c[i])
    cutpoints$one.specificity[i] <- 1-(cutpoints$b[i]/(cutpoints$b[i] + cutpoints$d[i]))
  }

  #Calculate 1-E cutpoints
  if((1-E)>max(cutpoints$sensitivity)){stop("Maximum sensitivity less than error threshold (E). Select a higher value for E.")}
  cutpoints2 <- cutpoints[cutpoints$sensitivity >= (1-E),]

  for(i in 1:dim(cutpoints2)[1]){
    cutpoints2$sensitivity[i]     <- cutpoints2$a[i]/(cutpoints2$a[i] + cutpoints2$c[i])
    cutpoints2$one.specificity[i] <- 1-(cutpoints2$b[i]/(cutpoints2$b[i] + cutpoints2$d[i]))
  }

  #Calculate the area under the null curve given E
  null <- ((min(cutpoints2$one.specificity)+1)/2)*(1-min(cutpoints2$one.specificity))

  #Calculate AUC value for each cutpoint
  for(i in 1:dim(cutpoints2)[1]){
    cutpoints2$AUC[i] <- ifelse(cutpoints2[i,1] == 0, (((1+cutpoints2$sensitivity[i])/2) * (1-cutpoints2$one.specificity[i])),
                               (((cutpoints2$sensitivity[i] + cutpoints2$sensitivity[i-1])/2) * (cutpoints2$one.specificity[i-1] - cutpoints$one.specificity[i])))
  }

  #Calculate modified AUC ratio
  AUC.ratio <- sum(cutpoints2$AUC)/null

  #Plot modified ROC
  plot(cutpoints2$one.specificity,cutpoints2$sensitivity,
       main = paste("Modified ROC Curve of ", n ," Best Models", sep = ""),
       type = "l", col = "blue", lwd = 2, tck = -0.02,
       xlab = "1-Specificity",
       ylab = "Sensitivity",
       ylim = c(0,1.0),
       xlim = c(0,1.0))
  lines(seq(max(cutpoints2$one.specificity), 1, length.out = 2),
        c(max(cutpoints2$sensitivity), max(cutpoints2$sensitivity)), col = "blue", lwd = 2)
  lines(seq(0, 1, 0.1), seq(0, 1, 0.1),lwd = 1)
  segments(0,min(cutpoints2$sensitivity), min(cutpoints2$one.specificity, min(cutpoints2$sensitivity)), lty = 2)
  segments(min(cutpoints2$one.specificity),0, min(cutpoints2$one.specificity), min(cutpoints2$sensitivity), lty = 2)

  legend(x=0.5, y = 0.5, legend=(c("Models", "Reference", "1-E")),
         bg = "", bty = "n", col = c("blue", "black", "black"), lty = c(1,1,2), lwd = c(2,2,2))

  return(list(Modified.AUC= sum(cutpoints2$AUC), Modified.Reference = null, AUC.Ratio = AUC.ratio))
}

