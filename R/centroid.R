#' Returns centroid coordinates of raster cells containing presence locations
#'
#' \code{centroid} The desktop version of GARP records the instance of a presence location in a grid cell, and therefore does not distinguish between multiple locations within a single cell. This function
#' returns the centroid location of all grid cells that contain a presence location, which circumvents multiple locations. It is useful when input rasters cells are large or
#' when sampling was unevenly clustered. This function uses the cell size and extent of user-defined raster to generate centroid coordinates of all cells that contain presence locations.
#'
#' @param x an object of class raster that represents the spatial resolution and extent of raster inputted into GARP
#' @param points a data.frame or SpatialPointsDataFrame object containing latitude and longitude values of presence locations
#' @param xy a vector object of latitude and longitude coordinates from \code{points} when \code{points} is class data.frame. Not needed if \code{points} is a SpatialPointsDataFrame object
#' @param species the vector of species names from the \code{points} data.frame or SpatialPointsDataFrame object
#'
#' @return A SpatialPointsDataFrame object containing "Latitude", "Longitude", and "Species" vectors.
#'
#' @details \code{x} should be a raster that represents the spatial resolution and extent of all environmental layers that will be inputted into GARP. Recommended to use a raster output by \code{\link{rasterPrep}}.
#'
#' @examples
#' set.seed(0)
#' library(raster)
#' r   <- raster(ncols = 100, nrows = 100)
#' r[] <- rbinom(5, 10, 0.3)
#' hs  <- data.frame("Latitude" = c(-89, 72, 63, 42, 54),"Longitude" = c(-12, 13, 24, 26, 87), "Species" = rep("Homo_sapiens", 5))
#' centroid(x = r, points = hs, xy = cbind(hs$Latitude, hs$Longitude), species = hs$Species)
#'
#' @import raster
#' @import sp
#' @export


centroid <- function(x, points, xy = NULL, species, output.path = NULL){

  if(class(points) == "data.frame"){
    #Convert dataframe to SpatialPointsDataFrame
    points.sp <- SpatialPointsDataFrame(as.matrix(xy), data = points)
  }else{points.sp <- points}

  #Extract the cell numbers from the grid file
  points.sp$cells <- extract(x, points.sp, cellnumbers = TRUE)[,1]

  #Obtain centroid coordinates
  points.sp$centroids <- xyFromCell(x, cell = points.sp$cells)

  points.df <- data.frame(points.sp)

  species <- levels(as.factor(species))

  points.all <- data.frame()

  for(i in species){
    points.sub <- subset(points.df, species == i)

    #Pull unique centroids and subset points
    points.uni <- subset(points.sub, !duplicated(points.sub$cells))

    points.all <- rbind(points.all, points.uni)
  }

  #Convert to spatial points dataframe
  centroids <- SpatialPointsDataFrame(as.matrix(points.all$centroids), data = points.all)
  return(centroids)

}
