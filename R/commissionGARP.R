#' Calculates commission metrics of model outputs
#'
#' \code{commissionGARP} Calculates the commission, or the proportion of the landscape predicted as present, of the best models output by GARP.
#'
#' @param n a numeric value specifying the number of models in the best subset output from GARP
#' @param x a raster object of the summated raster of the best models output from GARP
#' @param points a spatial object of presence data to use for testing locations
#' @param file.path file directory of folder containing best model rasters output from GARP (usually within athe "GARPruns folder")
#'
#' @return A list object containing:
#'  \itemize{
#'   \item \code{Average.Commission} the average commission of the model subset
#'   \item \code{Total.Commission} the total commission of the model subset
#'  }
#'
#' @details The discrete cutpoints (\code{n}) represent the number of models that agree on a predicted presence location.
#'
#' The raster object (\code{x}) should be a raster representing the number of models that agree on a predicted presence location per pixel and that outtput by \code{\link{sumRasters}}.
#'
#' The shapefile \code{points} should presence locations that were not used by GARP for model training and those output by \code{\link{splitData}}.
#'
#' @seealso \code{\link{omissionGARP}}
#'
#' @examples
#' \dontrun{
#' commissionGARP(n = 10, x = sum.raster, points = test_1, file.path = "C/GARPruns/BestSubsets/")
#' }
#'
#' @import raster
#' @import sp
#'
#' @export

commissionGARP <- function(n, x, points, file.path){

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

  #Calculate total commission
  total.commission <- cutpoints$cutpoint.area[n+1]/sum(cutpoints$cutpoint.area) * 100

  #Calculate average commission
  list.rasters <- list.files(path = paste(file.path), full.names=TRUE)

  #Remove info folder
  for(i in 1:length(list.rasters)){
    if("info" %in% strsplit(list.rasters[i],split="/")[[1]][length(strsplit(list.rasters[i], split = "/")[[1]])]){list.rasters <- list.rasters[-i]}
  }

  commission.df <- data.frame(rep(0, length(list.rasters)))
  names(commission.df) <- "task.name"
  for(i in 1:length(list.rasters)){
    commission.df$task.name[i] <- paste("task",i,sep="")
  }
  for(i in 1:dim(commission.df)[1]){
    task <- raster(list.rasters[i])
    commission.df$total.pixels[i] <- ncell(task) - freq(task, value=NA)
    commission.df$model.area[i]   <- freq(task, value=1)
    commission.df$idv.com[i]      <- commission.df$model.area[i]/commission.df$total.pixels[i] * 100
  }
  avg.commission <- mean(commission.df$idv.com)
  return(list(Average.Commission = avg.commission, Total.Commission = total.commission))
}
