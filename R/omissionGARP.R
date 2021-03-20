#' Calculates omission metrics of model outputs
#'
#' \code{omissionGARP} Calculates the omission, or the true occurrence points that are not predicted as present, of the best models output from GARP.
#'
#' @param n a numeric value specifying the number of models in the best subset output from GARP
#' @param x a raster object of the summated raster of the best models output from GARP
#' @param points a spatial object of presence data to use for testing locations
#'
#' @return A list object containing:
#'  \itemize{
#'   \item \code{Average.Omission} the average omission of the model subset
#'   \item \code{Total.Omission} the total omission of the model subset
#'  }
#'
#' @details The discrete cutpoints (\code{n}) represent the number of models that agree on a predicted presence location.
#'
#' The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'
#' The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}.
#'
#' @seealso \code{\link{commissionGARP}}
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   omissionGARP(n = 10, x = r, points = SpatialPoints(hs[,1:2]))
#'
#' @import raster
#' @import sp
#'
#' @export

omissionGARP <- function(n, x, points){

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
    cutpoints$cum.area[i] <- ifelse(cutpoints$cutpoint[i] == 0, 0,
                                    ifelse(cutpoints$cutpoint[i] == 1,
                                           (cutpoints$no.taxa.pixels[1] + cutpoints$no.taxa.pixels[i]),
                                           (cutpoints$no.taxa.pixels[i] + cutpoints$cum.area[i-1])))
  }

  #Calculate omission
  total.omission <- cutpoints$taxa.present[1]/sum(cutpoints$taxa.present)
  avg.omission <- (1 - ((sum(cutpoints$cutpoint * cutpoints$taxa.present) + total.omission)/(n * sum(cutpoints$taxa.present)))) * 100
  return(list(Total.Omission = total.omission, Average.Omission = avg.omission))
}
