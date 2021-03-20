#' Check that input points overlap with input rasters
#'
#' \code{checkData} Checks that the training or testing points overlap environmental layers used in GARP.
#'
#' @param points a spatial object of the presence data used for training GARP, testing output models, or both training and testing points
#' @param grid a raster object of any environmental layer that will be input into GARP
#' @param mask optional; a spatial object of the study or sampling area
#'
#' @return If any point falls outside raster, will return error to resample or recrop raster to allow overlay. If all points fall within raster, will return message that data is okay to use in GARP.
#'
#' @details The raster used for \code{grid} should either be cropped to study or sampling area, or \code{mask} should be used to specify area to crop raster.
#'
#' @examples
#'   set.seed(0)
#'   library(raster)
#'   r   <- raster(ncols = 100, nrows = 100)
#'   r[] <- rbinom(5, 10, 0.3)
#'   hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54), "Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#'   checkData(points = SpatialPoints(hs[,1:2]), grid = r)
#'
#' @import raster
#' @import sp
#'
#' @export

checkData <- function(points, grid, mask = NULL){
  #Read in mask to clip rasters
  if(is.null(mask) == FALSE){
    if(is.na(projection(grid) ==  projection(mask))){
      projection(mask) <- projection(grid)
      grid <- crop(grid, mask, snap = "out")
      grid <- mask(grid, mask = mask)
    }
    else if(projection(grid) !=  projection(mask)){
      projection(mask) <- projection(grid)
      grid <- crop(grid, mask, snap = "out")
      grid <- mask(grid, mask = mask)
    }else{
      grid <- crop(grid, mask, snap = "out")
      grid <- mask(grid, mask = mask)
    }
  }

  taxa.models <- raster::extract(grid, points)
  if(any(is.na(taxa.models))){stop("One or more testing points do not overlap raster. Suggest buffering mask.")}else{print("All points fall within cropped raster.")}

}
