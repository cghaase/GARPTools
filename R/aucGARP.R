#' Area Under the Receiver Operating Characteristic (ROC) curve (AUC)
#'
#' \code{aucGARP} Calculates the area under the ROC curve (AUC), which represents the accuracy of the subset of the best models produced by GARP. AUC evaluates the predictive performance of the
#' best models by relating the model sensitivity (true positive rate) to 1-specificity (true negative rate), and can be described as the probability that any given cell is correctly predicted
#' as present or absent. An AUC of 0.5 is that predicted at random, while an AUC of 1 represents a perfect prediction
#'
#' @param n a numeric value specifying the number of models in the best subset outputted by DesktopGARP
#' @param x a raster object of the summated raster of the best models output by DesktopGARP
#' @param points a spatial object of presence data to use for testing locations
#'
#' @return Plots the Receiver Operating Characteristic (ROC) curve and returns a list containing:
#'  \itemize{
#'   \item \code{AUC} the total area under the ROC curve
#'   \item \code{Wilcoxon} the Wilcoxon test statistic; essentially equal to the AUC
#'   \item \code{Standard.Error} the standard error of the AUC value
#'   \item \code{Z.score} the z-score associated with the AUC value
#'  }
#'
#' @details When calculating AUC for discrete cutpoints (\code{n}), represented by the number of models that agree on a predicted presence location, the AUC of the entire curve
#'   is essentially equal to the Wilcoxon test statistic.
#'   The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'   The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}.
#'
#' @seealso \code{\link{zAUC}}, \code{\link{seAUC}}, \code{\link{wGARP}}, \code{\link{plotROC}}
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   aucGARP(n = 10, x = r, points = SpatialPoints(hs[,1:2]))
#'
#' @references Cortes, C. and Mohri, M. (2004) Confidence Intervals for the Area Under the ROC curve. \emph{Advances in Neural Information Processing Systems}. \bold{6}.
#'
#' @import raster
#' @import sp
#' @export

aucGARP <- function(n, x, points){

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

  #Calculate AUC value for each cutpoint
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$AUC[i] <- ifelse(cutpoints[i,1] == 0, (((1+cutpoints$sensitivity[i])/2) * (1-cutpoints$one.specificity[i])),
                               (((cutpoints$sensitivity[i] + cutpoints$sensitivity[i-1])/2) * (cutpoints$one.specificity[i-1] - cutpoints$one.specificity[i])))
  }
  AUC.sum <- sum(cutpoints$AUC)

  #Plot ROC
  plot(cutpoints$one.specificity,cutpoints$sensitivity,
       main = paste("ROC Curve of ", n ," Best Models", sep = ""),
       type = "l", col = "blue", lwd = 2, tck = -0.02,
       xlab = "1-Specificity",
       ylab = "Sensitivity",
       ylim = c(0,1.0),
       xlim = c(0,1.0))
  lines(seq(max(cutpoints$one.specificity), 1, length.out = 2),
        c(max(cutpoints$sensitivity), max(cutpoints$sensitivity)), col = "blue", lwd = 2)
  lines(seq(0, 1, 0.1), seq(0, 1, 0.1),lty = 2, lwd = 1)
  points(cutpoints$one.specificity, cutpoints$sensitivity, pch = 16)
  legend("bottomright", legend=(c("Models", "Reference")), cex = 0.85,
         bg = "", bty = "n", col = c("blue", "black"), lty = c(1,2), lwd = c(2,2))

  #Calculate the z-score
  mu.exp <- (sum(cutpoints$taxa.present) + sum(cutpoints$no.taxa.pixels))/2
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$mid.rank[i] <- ifelse(cutpoints[i,1] == 0, (0.5 * (cutpoints$cutpoint.area[i] + 1)),
                                    (cutpoints$mid.rank[i-1] + (0.5 * (cutpoints$cutpoint.area[i] + 1))))
    cutpoints$mu[i] <- cutpoints$taxa.present[i] * cutpoints$mid.rank[i]
    cutpoints$d3[i] <- ((cutpoints$taxa.present[i] + cutpoints$no.taxa.pixels[i])^3) - (cutpoints$taxa.present[i] + cutpoints$no.taxa.pixels[i])
  }
  mu.obs <- sum(cutpoints$mu)
  var <- ((sum(cutpoints$taxa.present) * sum(cutpoints$no.taxa.pixels) * (sum(cutpoints$taxa.present) + sum(cutpoints$no.taxa.pixels) + 1))/12) -
    ((sum(cutpoints$taxa.present)* sum(cutpoints$no.taxa.pixels) * sum(cutpoints$d3))/(12*(sum(cutpoints$taxa.present) +
                                                                                             sum(cutpoints$no.taxa.pixels)) * ((sum(cutpoints$taxa.present) + sum(cutpoints$no.taxa.pixels)) - 1)))
  stdev <- sqrt(var)
  z <- (mu.obs-mu.exp)/stdev

  #Calculate the Wilcoxon statistic
  for(i in 1:dim(cutpoints)[1]){
    cutpoints$W[i] <- (cutpoints$no.taxa.pixels[i] * cutpoints$a[i]) + (0.5*cutpoints$no.taxa.pixels[i] * cutpoints$taxa.present[i])
  }
  total.w <- sum(cutpoints$W)/(sum(cutpoints$no.taxa.pixels) * sum(cutpoints$taxa.present))
  total.w.sq <- total.w^2

  #Calculate the standard error
  Q1 <- total.w/(2 - total.w)
  Q2 <- (2 * (total.w.sq))/(1 + total.w)
  SE <- sqrt(((total.w * (1 - total.w)) + ((sum(cutpoints$taxa.present) - 1) *
                                             (Q1 - total.w.sq)) + ((sum(cutpoints$no.taxa.pixels) - 1) *
                                                                     (Q2 - total.w.sq)))/(sum(cutpoints$taxa.present) * sum(cutpoints$no.taxa.pixels)))
  return(list(AUC= AUC.sum, Wilcoxon = total.w, Standard.Error = SE, Z.score = z))
}

