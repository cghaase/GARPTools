#' Calculates the confusion matrix
#'
#' \code{confuseMatrix} Returns the confusion matrix containing the sensitivity and specificity of the best model subset of a GARP run.
#'
#' @param n a numeric value specifying the number of models in the best subset output from GARP
#' @param x a raster object of the summated raster of the best models output from GARP
#' @param points a spatial object of presence data to use for testing locations
#'
#' @return A data.frame object containing:
#'  \itemize{
#'   \item \code{a} True positive
#'   \item \code{b} False negative
#'   \item \code{c} False positive
#'   \item \code{d} True negative
#'   \item \code{sensitivity} the proportion of occupied pixels correctly predicted as present
#'   \item \code{specificity} the proportion of non-occupied pixels correctly predicted as absent
#'   \item \code{1-specificity} the proportion of non-occupied pixels incorrectly predicted as present
#'  }
#'
#' @details For discrete cutpoints (\code{n}), represents the number of models that agree on a predicted presence location.
#'
#' The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'
#' The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}
#'
#' @seealso \code{\link{aucGARP}}
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   confuseMatrix(n = 10, x = r, points = SpatialPoints(hs[,1:2]))
#'
#' @import raster
#' @import sp
#'
#' @export

confuseMatrix <- function(n, x, points){

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

  df <- data.frame(cutpoints$cutpoint, cutpoints$a, cutpoints$b, cutpoints$c, cutpoints$d, cutpoints$sensitivity, 1-cutpoints$one.specificity, cutpoints$one.specificity)
  colnames(df) <- c("Models", "a", "b", "c", "d", "sensitivity", "specificity", "1-specificity")
  return(df)
}
